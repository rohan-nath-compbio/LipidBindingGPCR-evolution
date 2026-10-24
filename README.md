# Large-Scale Comparative Analysis Uncovers the Evolutionary History and Hidden Diversity of Lipid-Responsive GPCRs

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R: >=4.1](https://img.shields.io/badge/R-%3E%3D4.1-276DC3.svg?logo=r)](https://www.r-project.org/)
[![Python: >=3.10](https://img.shields.io/badge/Python-%3E%3D3.10-3776AB.svg?logo=python)](https://www.python.org/)

Reproducible analysis workflows, scripts, and curated datasets supporting the manuscript:
> *Large-Scale Comparative Analysis Uncovers the Evolutionary History and Hidden Diversity of Lipid-Responsive GPCRs*

---

## Overview

This repository provides reproducible code and curated input data supporting the comparative genomic analyses described in the manuscript. Specifically, the repository includes routines to:

- Reconstruct phyletic distribution patterns across vertebrate and invertebrate clades (taxonomic distributions, copy-number modes, sequence counts).
- Perform micro-synteny analyses to map conserved genomic neighborhoods surrounding lipid-responsive GPCR loci (e.g., *CNR1*).
- Systematically retrieve and filter animal orthology datasets from the public KEGG API.
- Execute an end-to-end command-line workflow generating taxonomically arranged gene-neighborhood reports directly from NCBI protein accessions.

---

## CNR1 Gene-Neighborhood Architecture (Synteny Demo)

The figure below illustrates the micro-syntenic environment surrounding the **CNR1** locus across representative chordate genomes, rendered dynamically via `geneviewer`:

<p align="center">
  <img src="results/synteny/synteny.svg" alt="CNR1 Gene Neighborhood Synteny" width="100%">
</p>

*Figure: Synteny map showing the target gene (CNR1, crimson), conserved paralogon orthologs (pink), surrounding genomic neighborhood (grey), strand orientations, and aligned chromosomal coordinates.*

---

## Repository Structure

```text
LipidBindingGPCR-evolution/
├── data/
│   ├── kegg/                   # Curated KEGG ortholog snapshots across animal phyla
│   ├── phyletic-distribution/  # Input datasets for taxonomic counts and copy-number modes
│   └── synteny/                # Curated CNR1 synteny dataset and example accession list
├── scripts/
│   ├── kegg/                   # KEGG REST retrieval utility (fetch_kegg_orthologs.py)
│   ├── phyletic-distribution/  # R scripts generating manuscript distribution heatmaps
│   └── synteny/                # Synteny visualization script (synteny.R) and pipeline
└── results/
    └── synteny/                # Generated publication figures (synteny.svg)
```

---

## Requirements

Run all commands from the repository root.

| Workflow | Requirements |
| :--- | :--- |
| **Phyletic-distribution figures** | R ≥ 4.1, `ggplot2`, `dplyr`, `egg`, `tidyr`, `stringr`, `scales`, `svglite` |
| **CNR1 synteny figure** | R ≥ 4.1, `geneviewer`, `readr`, `htmlwidgets`, `chromote`, and Google Chrome / Chromium |
| **KEGG retrieval** | Python ≥ 3.10 (Standard Library only; no third-party dependencies required) |
| **Gene-neighborhood pipeline** | Linux / macOS, Conda, Bash, Python ≥ 3.10, `curl`, `tar`, `gzip`, `awk`, `grep`, `sed`, `sort` |

### R Package Installation

Install required R packages once via CRAN:

```r
install.packages(c(
  "ggplot2", "dplyr", "egg", "tidyr", "stringr", "scales", "svglite",
  "readr", "geneviewer", "htmlwidgets", "chromote"
))
```

---

## Reproduce Manuscript Figures

Execute the following commands from the repository root:

```bash
# Generate phyletic distribution heatmaps
Rscript scripts/phyletic-distribution/plot_taxonomic_distribution_vertebrates.R
Rscript scripts/phyletic-distribution/plot_taxonomic_distribution_invertebrates.R
Rscript scripts/phyletic-distribution/plot_copy_number_mode_vertebrates.R
Rscript scripts/phyletic-distribution/plot_sequence_counts_vertebrates.R

# Generate CNR1 synteny figure
Rscript scripts/synteny/synteny.R
```

The scripts create output directories in `results/` as needed and generate the following SVG files:

| Analysis | Script | Output File |
| :--- | :--- | :--- |
| Vertebrate taxonomic distribution | `plot_taxonomic_distribution_vertebrates.R` | `results/phyletic-distribution/taxonomic_distribution_vertebrates.svg` |
| Invertebrate taxonomic distribution | `plot_taxonomic_distribution_invertebrates.R` | `results/phyletic-distribution/taxonomic_distribution_invertebrates.svg` |
| Vertebrate copy-number mode | `plot_copy_number_mode_vertebrates.R` | `results/phyletic-distribution/copy_number_mode_vertebrates.svg` |
| Vertebrate sequence counts | `plot_sequence_counts_vertebrates.R` | `results/phyletic-distribution/sequence_counts_vertebrates.svg` |
| CNR1 gene neighborhood | `synteny.R` | `results/synteny/synteny.svg` |

---

## Generate Gene Neighborhoods from Accessions

Provide an input `.acc` file containing one NCBI protein accession per line. To run the supplied test accession set:

```bash
bash scripts/synteny/run_synteny_pipeline.sh data/synteny/example.acc
```

> **Memory recommendation:** Run the initial database build on a system with **more than 32 GB RAM**. The build parses large NCBI annotation and accession files in memory and performs a parallel sort; lower-memory systems can run out of memory during these stages. Subsequent runs reuse the completed database and have substantially lower resource requirements.

For full pipeline options, input specifications, and troubleshooting, see the [Gene-Neighborhood Pipeline Guide](scripts/synteny/README.md).

---

## Regenerate KEGG Ortholog Tables

The deposited tables in [`data/kegg`](data/kegg) are reproducible snapshots. To query the public KEGG REST API and regenerate a table for a gene symbol or KEGG Orthology (KO) identifier:

```bash
python3 scripts/kegg/fetch_kegg_orthologs.py DAGL data/kegg/dagl.tsv
```

The script includes automated rate-limiting, exponential backoff retries, and animal taxonomy fallbacks. See [data/kegg/README.md](data/kegg/README.md) for column definitions and scope.

---

## Input Data

| Dataset | Description |
| :--- | :--- |
| `data/phyletic-distribution/receptor_taxonomic_distribution.csv` | Unique species counts by receptor and taxonomic group |
| `data/phyletic-distribution/receptor_copy_number_mode.csv` | Modal receptor copy number by receptor and taxonomic group |
| `data/phyletic-distribution/receptor_sequence_counts.csv` | Sequence counts by receptor and taxonomic group |
| `data/synteny/synteny.csv` | Curated genomic neighborhood coordinates used for synteny mapping |
| `data/synteny/example.acc` | Example accession input for the command-line gene-neighborhood workflow |
| `data/kegg/*.tsv` | KEGG-derived ortholog snapshots |

In the phyletic-distribution tables, `Phyla` identifies the taxonomic group and `Receptor` identifies a receptor or pipe-delimited receptor block. A `/LIKE` suffix denotes a receptor-like sequence.

---

## Citation

If you use these materials, please cite the associated manuscript:
> *Large-Scale Comparative Analysis Uncovers the Evolutionary History and Hidden Diversity of Lipid-Responsive GPCRs.* (Full bibliographic information and DOI will be added upon publication).

---

## License

Code and data are released under the [MIT License](LICENSE).
