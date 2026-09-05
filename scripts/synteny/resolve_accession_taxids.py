#!/usr/bin/env python3
"""Resolve NCBI protein accessions to taxonomy IDs using E-utilities."""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen


EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve protein accession IDs in a .acc file to NCBI taxonomy IDs."
    )
    parser.add_argument("input_acc", type=Path, help="One protein accession per line.")
    parser.add_argument("output", type=Path, help="Output file with: accession<TAB>taxid.")
    parser.add_argument("--email", default=os.getenv("NCBI_EMAIL"), help="Contact email (or set NCBI_EMAIL).")
    parser.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"), help="NCBI API key (or set NCBI_API_KEY).")
    parser.add_argument("--batch-size", type=int, default=200, help="Accessions per request (default: 200).")
    parser.add_argument("--retries", type=int, default=5, help="Retries per failed request (default: 5).")
    return parser.parse_args()


def read_accessions(path: Path) -> list[str]:
    accessions = []
    for line in path.read_text().splitlines():
        accession = line.strip()
        if accession and not accession.startswith("#"):
            accessions.append(accession)
    if not accessions:
        raise ValueError(f"No accession IDs found in {path}.")
    return list(dict.fromkeys(accessions))


def completed_accessions(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return {line.split("\t", 1)[0] for line in path.read_text().splitlines() if line.strip()}


def request_batch(ids: list[str], email: str, api_key: str | None, retries: int) -> dict[str, str]:
    parameters = {"db": "protein", "id": ",".join(ids), "retmode": "json", "email": email, "tool": "lipidbindinggpcr-synteny"}
    if api_key:
        parameters["api_key"] = api_key
    url = f"{EUTILS_URL}?{urlencode(parameters)}"

    for attempt in range(retries):
        try:
            with urlopen(url, timeout=60) as response:
                result = json.load(response).get("result", {})
            resolved = {}
            for uid in result.get("uids", []):
                record = result.get(uid, {})
                taxid = str(record.get("taxid", "")).strip()
                if taxid and taxid != "0":
                    for accession in (record.get("accessionversion", ""), record.get("caption", "")):
                        if accession:
                            resolved[accession] = taxid
            return resolved
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as error:
            if attempt == retries - 1:
                raise RuntimeError(f"NCBI request failed for {ids[0]}: {error}") from error
            time.sleep(2 ** attempt)
    return {}


def main() -> int:
    args = arguments()
    if not args.email:
        print("Error: provide --email or set NCBI_EMAIL.", file=sys.stderr)
        return 2
    if args.batch_size < 1:
        print("Error: --batch-size must be positive.", file=sys.stderr)
        return 2

    accessions = read_accessions(args.input_acc)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    done = completed_accessions(args.output)
    pending = [accession for accession in accessions if accession not in done]
    print(f"Accessions: {len(accessions)}; already resolved: {len(done)}; pending: {len(pending)}", file=sys.stderr)

    delay = 0.11 if args.api_key else 0.34
    unresolved = []
    with args.output.open("a", encoding="utf-8") as handle:
        for start in range(0, len(pending), args.batch_size):
            batch = pending[start:start + args.batch_size]
            resolved = request_batch(batch, args.email, args.api_key, args.retries)
            for accession in batch:
                taxid = resolved.get(accession) or resolved.get(accession.split(".", 1)[0])
                if taxid:
                    handle.write(f"{accession}\t{taxid}\n")
                else:
                    unresolved.append(accession)
            handle.flush()
            time.sleep(delay)

    if unresolved:
        print(f"Warning: {len(unresolved)} accession(s) could not be resolved; rerun to retry them.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
