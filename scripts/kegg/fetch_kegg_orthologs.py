#!/usr/bin/env python3
"""Retrieve KEGG ortholog members for the animal groups analysed here.

The script resolves a gene symbol or KEGG Orthology (KO) identifier, filters
members to the target phyla, and writes their KEGG and NCBI protein identifiers
as a tab-separated table.
"""

import time
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
import socket
import sys
from dataclasses import dataclass
from threading import Lock
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


KEGG_REST = "https://rest.kegg.jp"
USER_AGENT = "Mozilla/5.0 (compatible; kegg-fetch-script/1.0)"
# KEGG asks API users to make no more than three requests per second.  This is
# shared by all worker threads, because conversion requests are concurrent.
KEGG_MIN_REQUEST_INTERVAL = 0.4
_request_lock = Lock()
_last_request_started = 0.0

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


def kegg_text(path: str, timeout: int = 120, retries: int = 5) -> str:
    global _last_request_started

    url = f"{KEGG_REST}/{path.lstrip('/')}"
    request = Request(url, headers={"User-Agent": USER_AGENT})

    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with _request_lock:
                wait = KEGG_MIN_REQUEST_INTERVAL - (time.monotonic() - _last_request_started)
                if wait > 0:
                    time.sleep(wait)
                _last_request_started = time.monotonic()
            with urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8")
        except HTTPError as exc:
            # KEGG can temporarily return 400 when an IP has been rate-limited.
            # Do not retry ordinary permanent client errors beyond this small
            # bounded retry loop.
            if exc.code in (400, 408, 429, 500, 502, 503, 504):
                last_error = exc
                if attempt < retries - 1:
                    time.sleep(2.0 ** (attempt + 1))
                    continue
                break
            raise RuntimeError(f"KEGG request failed ({exc.code}): {url}") from exc
        except (URLError, TimeoutError, socket.timeout) as exc:
            last_error = exc
            if attempt < retries - 1:
                time.sleep(2.0 ** (attempt + 1))
                continue
            reason = getattr(exc, "reason", exc)
            raise RuntimeError(
                f"Could not reach KEGG after {retries} attempts: {reason}"
            ) from exc

    hint = ""
    if isinstance(last_error, HTTPError) and last_error.code in (400, 429):
        hint = " KEGG may have temporarily rate-limited this IP; wait a few minutes and retry."
    raise RuntimeError(f"KEGG request failed after {retries} attempts: {url}.{hint}") from last_error

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch KEGG ortholog members for selected animal phyla as TSV."
    )
    parser.add_argument(
        "query",
        help="KEGG Orthology ID, e.g. K13806. A gene symbol can also be used.",
    )
    parser.add_argument(
        "output_tsv",
        nargs="?",
        help="Optional output TSV path; omit to write TSV to standard output.",
    )
    return parser.parse_args()


def load_organisms() -> dict[str, Organism]:
    organisms: dict[str, Organism] = {}

    try:
        organism_listing = kegg_text("list/organism")
    except RuntimeError:
        return load_animal_genomes()

    for line in organism_listing.splitlines():
        columns = line.split("\t")
        if len(columns) < 4:
            continue

        _, code, species, lineage = columns[:4]
        organisms[code] = Organism(code=code, species=species, lineage=lineage)

    if not organisms:
        raise RuntimeError("KEGG returned no organisms.")

    return organisms


def load_animal_genomes() -> dict[str, Organism]:
    """Load animal taxonomy from GENOME when the special organism list fails.

    KEGG documents ``list/genome/animals`` and permits up to ten genome entries
    in each ``get`` request.  Genome records contain the LINEAGE field needed
    to classify genes into the target phyla.
    """
    genome_ids = [
        line.split("\t", 1)[0]
        for line in kegg_text("list/genome/animals").splitlines()
        if line.startswith("T") and "\t" in line
    ]
    organisms: dict[str, Organism] = {}

    for start in range(0, len(genome_ids), 10):
        entries = kegg_text(f"get/{'+'.join(genome_ids[start:start + 10])}")
        for entry in entries.split("///"):
            code_match = re.search(r"^ORG_CODE\s+(\S+)", entry, flags=re.MULTILINE)
            name_match = re.search(r"^NAME\s+(.+)$", entry, flags=re.MULTILINE)
            lineage_match = re.search(
                r"^  LINEAGE\s+(.+)$", entry, flags=re.MULTILINE
            )
            group_match = re.search(r"^ORG_GROUP\s+(.+)$", entry, flags=re.MULTILINE)
            if not (code_match and name_match and lineage_match):
                continue

            code = code_match.group(1)
            organisms[code] = Organism(
                code=code,
                species=name_match.group(1).strip(),
                lineage="; ".join(
                    part.group(1).strip()
                    for part in (lineage_match, group_match)
                    if part is not None
                ),
            )

    if not organisms:
        raise RuntimeError("KEGG returned no animal genome records.")

    return organisms


def normalize_species(species: str) -> str:
    scientific_name = species.split(" (", 1)[0].strip()
    return re.sub(r"\s+", "_", scientific_name)


def numeric_kegg_id(gene_id: str) -> str:
    return gene_id.split(":", 1)[1] if ":" in gene_id else gene_id


def protein_ids_for_genes(
    gene_ids: list[str], batch_size: int = 10, max_workers: int = 2
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


def write_rows(output_tsv: str | None, rows: list[list[str]]) -> None:
    if output_tsv is None:
        writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        writer.writerow(["Phylum", "Species", "KeggID", "NCBI-ProteinID"])
        for row in rows:
            writer.writerow([row[0], row[1], row[2], row[5]])
        return

    with open(output_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Phylum", "Species", "KeggID", "NCBI-ProteinID"])
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

    destination = args.output_tsv or "standard output"
    print(f"Wrote {len(rows)} rows to {destination}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
