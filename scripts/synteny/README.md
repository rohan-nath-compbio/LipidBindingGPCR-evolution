# Gene-neighborhood pipeline

`run_synteny_pipeline.sh` converts a list of NCBI protein accessions into a taxonomically arranged gene-neighborhood report. It resolves taxonomy identifiers, obtains complete TaxonKit lineages, groups accessions by major vertebrate clades and invertebrates, and extracts local genomic neighborhoods.

## Quick start

From the repository root, run the supplied example:

```bash
bash scripts/synteny/run_synteny_pipeline.sh data/synteny/example.acc
```

The accession-resolution step may prompt for the request metadata required by NCBI. For non-interactive execution, consult `bash scripts/synteny/run_synteny_pipeline.sh --help` before running the workflow.

## Input format

The input must use the `.acc` extension and contain one NCBI protein accession per line. Empty lines and lines beginning with `#` are ignored.

```text
NP_149421.2
NP_031752.1
XP_038540064.1
```

## What the workflow does

1. Installs a project-local TaxonKit launcher when one is not available.
2. Downloads and builds the local NCBI gene-neighborhood database when required.
3. Resolves input accessions to NCBI taxonomy identifiers.
4. Retrieves complete taxonomic lineages with TaxonKit and replaces spaces with underscores.
5. Extracts neighborhoods, assigns them to taxonomic groups, and writes the arranged report.

The first run downloads public reference data and builds a local database. It can take substantial time and disk space; later runs reuse these resources.

> **Memory recommendation:** Perform the initial database build on a system with **more than 32 GB RAM**. The build keeps large NCBI accession and annotation mappings in memory while also running a parallel sort, so systems with less memory can fail during parsing or sorting. Once `gene_neighborhood.tsv` has been created, later runs reuse it and require substantially less memory.

## Requirements

The workflow is intended for Linux and macOS. It requires:

- Conda
- Bash and Python 3.10 or later
- `curl`, `tar`, `gzip`, `awk`, `grep`, `sed`, and `sort`
- Internet access for first-run setup and accession resolution

Python dependencies are limited to the standard library.

## Output

For an input named `cnr1.acc`, a successful run creates:

```text
results/synteny/cnr1/cnr1.arranged.neighbor
```

This report is arranged under taxonomic headers, including `MAMMALIA`, other vertebrate groups, and `INVERTEBRATE`. To keep result directories concise, intermediate accession, taxonomy, lineage, and neighborhood files are removed automatically after the final report is written. If report generation does not complete, intermediates are retained for inspection and rerunning.

## Reproducibility notes

The pipeline uses local copies of the required taxonomy and gene-neighborhood resources under `resources/synteny/`, and installs its TaxonKit launcher under `tools/bin/`. These derived resources and output reports are intentionally excluded from version control.

The helper scripts in this directory are invoked automatically by `run_synteny_pipeline.sh`; normally, run only the command shown above.

## Troubleshooting

- Run the command from the repository root so that its input and output paths resolve correctly.
- Confirm that Conda is available on `PATH` before the first run.
- Use `bash scripts/synteny/run_synteny_pipeline.sh --help` to check command usage.
- If a run ends before the final report is created, inspect the retained files in `results/synteny/<input-name>/` and rerun the same command.
