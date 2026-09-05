#!/usr/bin/env python3
"""Retrieve KEGG ortholog members for the animal groups analysed here.

The script resolves a gene symbol or KEGG Orthology (KO) identifier, filters
members to the target phyla, and writes their KEGG and NCBI protein identifiers
as a tab-separated table.
"""

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
import sys
from dataclasses import dataclass
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlopen


KEGG_REST = "https://rest.kegg.jp"

TARGET_PHYLA = [
    ("Vertebrata", ("Vertebrates", "Vertebrata", "Mammals", "Birds", "Reptiles", "Amphibians", "Fishes", "Cartilaginous fishes")),
    ("Tunicata", ("Tunicates", "Tunicata", "Urochordata")),
    ("Cephalochordata", ("Cephalochordates", "Cephalochordata")),
    ("Hemichordata", ("Hemichordates", "Hemichordata")),
    ("Echinodermata", ("Echinoderms", "Echinodermata")),
    ("Annelida", ("Annelids", "Annelida")),
    ("Mollusca", ("Mollusks", "Mollusca")),
    ("Arthropoda", ("Arthropods", "Arthropoda", "Insects", "Crustaceans", "Arachnids", "Myriapods", "Chelicerates")),
    ("Cnidaria", ("Cnidarians", "Cnidaria")),
    ("Placozoa", ("Placozoans", "Placozoa")),
]


@dataclass(frozen=True)
class Organism:
    code: str
    species: str
    lineage: str


def kegg_text(path: str, timeout: int = 60) -> str:
    url = f"{KEGG_REST}/{path.lstrip('/')}"
    try:
        with urlopen(url, timeout=timeout) as response:
            return response.read().decode("utf-8")
    except HTTPError as exc:
        raise RuntimeError(f"KEGG request failed ({exc.code}): {url}") from exc
    except URLError as exc:
        raise RuntimeError(f"Could not reach KEGG: {exc.reason}") from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch KEGG ortholog members for selected animal phyla and write a TSV."
    )
    parser.add_argument(
        "query",
        help="KEGG Orthology ID, e.g. K13806. A gene symbol can also be used.",
    )
    parser.add_argument("output_tsv", help="Output TSV path, e.g. dagl.tsv")
    return parser.parse_args()


def load_organisms() -> dict[str, Organism]:
    organisms: dict[str, Organism] = {}

    for line in kegg_text("list/organism").splitlines():
        columns = line.split("\t")
        if len(columns) < 4:
            continue

        _, code, species, lineage = columns[:4]
        organisms[code] = Organism(code=code, species=species, lineage=lineage)

    if not organisms:
        raise RuntimeError("KEGG returned no organisms.")

    return organisms


def normalize_species(species: str) -> str:
    scientific_name = species.split(" (", 1)[0].strip()
    return re.sub(r"\s+", "_", scientific_name)


def numeric_kegg_id(gene_id: str) -> str:
    return gene_id.split(":", 1)[1] if ":" in gene_id else gene_id


def protein_ids_for_genes(
    gene_ids: list[str], batch_size: int = 10, max_workers: int = 8
) -> dict[str, str]:
    def fetch_batch(batch: list[str]) -> dict[str, str]:
        batch_ids: dict[str, str] = {}
        text = kegg_text(f"conv/ncbi-proteinid/{'+'.join(batch)}", timeout=60)

        for line in text.splitlines():
            columns = line.split("\t")
            if len(columns) != 2:
                continue

            gene_id, protein_id = columns
            batch_ids[gene_id] = protein_id.removeprefix("ncbi-proteinid:")

        return batch_ids

    batches = [
        gene_ids[start:start + batch_size]
        for start in range(0, len(gene_ids), batch_size)
    ]
    protein_ids: dict[str, str] = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(fetch_batch, batch) for batch in batches]
        for future in as_completed(futures):
            protein_ids.update(future.result())

    return protein_ids


def lineage_parts(lineage: str) -> set[str]:
    return {part.strip().casefold() for part in lineage.split(";") if part.strip()}


def matching_phylum(lineage: str) -> str | None:
    parts = lineage_parts(lineage)

    for label, aliases in TARGET_PHYLA:
        if any(alias.casefold() in parts for alias in aliases):
            return label

    return None


def query_to_ko(query: str) -> str:
    query = query.strip()
    if re.fullmatch(r"K\d{5}", query, flags=re.IGNORECASE):
        return query.upper()

    hits = kegg_text(f"find/genes/{quote(query)}", timeout=30).splitlines()
    if not hits:
        raise RuntimeError(f"No KEGG gene hits found for {query!r}.")

    selected_gene = None
    for hit in hits:
        gene = hit.split("\t", 1)[0]
        if gene.startswith("hsa:"):
            selected_gene = gene
            break

    selected_gene = selected_gene or hits[0].split("\t", 1)[0]
    entry = kegg_text(f"get/{selected_gene}", timeout=30)

    match = re.search(r"^ORTHOLOGY\s+(K\d{5})", entry, flags=re.MULTILINE)
    if not match:
        match = re.search(r"\b(K\d{5})\b", entry)
    if not match:
        raise RuntimeError(f"Could not determine a KO for {query!r}.")

    return match.group(1)


def fetch_gene_links(ko: str) -> list[str]:
    genes: list[str] = []

    for line in kegg_text(f"link/genes/{ko}").splitlines():
        columns = line.split("\t")
        if len(columns) != 2:
            continue
        genes.append(columns[1])

    if not genes:
        raise RuntimeError(f"No KEGG gene links found for {ko}.")

    return genes


def build_rows(ko: str, organisms: dict[str, Organism]) -> list[list[str]]:
    rows: list[list[str]] = []

    for gene_id in fetch_gene_links(ko):
        organism_code = gene_id.split(":", 1)[0]
        organism = organisms.get(organism_code)
        if organism is None:
            continue

        phylum = matching_phylum(organism.lineage)
        if phylum is None:
            continue

        rows.append(
            [
                phylum,
                normalize_species(organism.species),
                numeric_kegg_id(gene_id),
                gene_id,
                organism_code,
            ]
        )

    protein_ids = protein_ids_for_genes([row[3] for row in rows])
    for row in rows:
        row.append(protein_ids.get(row[3], ""))

    phylum_order = {name: index for index, (name, _) in enumerate(TARGET_PHYLA)}
    rows.sort(key=lambda row: (phylum_order[row[0]], row[1], row[2]))
    return rows


def write_rows(output_tsv: str, rows: list[list[str]]) -> None:
    with open(output_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Phylum", "Species", "KEGG_GeneID", "NCBI_ProteinID"])
        for row in rows:
            writer.writerow([row[0], row[1], row[2], row[5]])


def main() -> int:
    args = parse_args()

    try:
        ko = query_to_ko(args.query)
        print(f"Using KO: {ko}", file=sys.stderr)

        organisms = load_organisms()
        rows = build_rows(ko, organisms)
        write_rows(args.output_tsv, rows)

    except RuntimeError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {len(rows)} rows to {args.output_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
