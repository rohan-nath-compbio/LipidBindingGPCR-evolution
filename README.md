# Evolution of lipid-responsive GPCRs

Reproducible code and curated input data supporting the manuscript *Large-Scale Comparative Analysis Uncovers the Evolutionary History and Hidden Diversity of Lipid-Responsive GPCRs*.

The repository generates the phyletic-distribution and CNR1 gene-neighborhood figures, retrieves comparable KEGG ortholog tables, and provides a command-line workflow for producing taxonomically arranged gene neighborhoods from protein accessions.

## Contents

```text
data/
├── kegg/                   KEGG-derived ortholog snapshots
├── phyletic-distribution/  Tables underlying the phyletic-distribution figures
└── synteny/                CNR1 neighborhood data and example accession input
scripts/
├── kegg/                   KEGG ortholog retrieval utility
├── phyletic-distribution/  R scripts for manuscript figures
└── synteny/                CNR1 figure script and neighborhood pipeline
results/                    Generated locally; not versioned
```

## Requirements

Run all commands from the repository root.

| Workflow | Requirements |
| --- | --- |
| Phyletic-distribution figures | R 4.1+ and the packages listed below |
| CNR1 synteny figure | The R requirements above, plus Chrome or Chromium available to `chromote` |
| KEGG retrieval | Python 3.10+; no third-party Python packages |
| Gene-neighborhood pipeline | Linux or macOS, Conda, Bash, Python 3.10+, `curl`, `tar`, `gzip`, `awk`, `grep`, `sed`, and `sort` |

Install the R packages once:

```r
install.packages(c(
  "ggplot2", "dplyr", "egg", "tidyr", "stringr", "scales", "svglite",
  "readr", "geneviewer", "htmlwidgets", "chromote"
))
```

## Reproduce the manuscript figures

```bash
Rscript scripts/phyletic-distribution/plot_taxonomic_distribution_vertebrates.R
Rscript scripts/phyletic-distribution/plot_taxonomic_distribution_invertebrates.R
Rscript scripts/phyletic-distribution/plot_copy_number_mode_vertebrates.R
Rscript scripts/phyletic-distribution/plot_sequence_counts_vertebrates.R
Rscript scripts/synteny/synteny.R
```

The scripts create `results/` as needed and write the following SVG files:

| Analysis | Output |
| --- | --- |
| Vertebrate taxonomic distribution | `results/phyletic-distribution/taxonomic_distribution_vertebrates.svg` |
| Invertebrate taxonomic distribution | `results/phyletic-distribution/taxonomic_distribution_invertebrates.svg` |
| Vertebrate copy-number mode | `results/phyletic-distribution/copy_number_mode_vertebrates.svg` |
| Vertebrate sequence counts | `results/phyletic-distribution/sequence_counts_vertebrates.svg` |
| CNR1 gene neighborhood | `results/synteny/synteny.svg` |

## Generate gene neighborhoods from accessions

Provide a `.acc` file containing one NCBI protein accession per line. The supplied example can be run with:

```bash
bash scripts/synteny/run_synteny_pipeline.sh data/synteny/example.acc
```

On its first run, the workflow installs TaxonKit into the project, downloads the required public NCBI resources, and builds a local gene-neighborhood database. This setup can require substantial time and disk space. A successful run writes one final report to `results/synteny/<input-name>/<input-name>.arranged.neighbor`.

> **Memory recommendation:** Run the initial database build on a system with **more than 32 GB RAM**. The build parses large NCBI annotation and accession files in memory and performs a parallel sort; lower-memory systems can run out of memory during these stages. Subsequent runs reuse the completed database and have substantially lower resource requirements.

See the [gene-neighborhood pipeline guide](scripts/synteny/README.md) for input format, workflow behavior, and troubleshooting.

## Regenerate KEGG ortholog tables

The deposited tables in [`data/kegg`](data/kegg) are reproducible snapshots. To retrieve a table for a gene symbol or KEGG Orthology identifier:

```bash
python3 scripts/kegg/fetch_kegg_orthologs.py DAGL data/kegg/dagl.tsv
```

The utility queries the public [KEGG REST API](https://rest.kegg.jp/) and writes `Phylum`, `Species`, `KEGG_GeneID`, and `NCBI_ProteinID` columns. Because KEGG is updated continuously, a regenerated table can differ from the deposited snapshot. See [data/kegg/README.md](data/kegg/README.md) for column definitions and scope.

## Input data

| Data set | Purpose |
| --- | --- |
| `data/phyletic-distribution/receptor_taxonomic_distribution.csv` | Unique species counts by receptor and taxonomic group |
| `data/phyletic-distribution/receptor_copy_number_mode.csv` | Modal receptor copy number by receptor and taxonomic group |
| `data/phyletic-distribution/receptor_sequence_counts.csv` | Sequence counts by receptor and taxonomic group |
| `data/synteny/synteny.csv` | Curated neighborhoods used for the CNR1 figure |
| `data/synteny/example.acc` | Example input for the command-line gene-neighborhood workflow |
| `data/kegg/*.tsv` | KEGG-derived ortholog snapshots |

In the phyletic-distribution tables, `Phyla` identifies the taxonomic group and `Receptor` identifies a receptor or pipe-delimited receptor block. A `/LIKE` suffix denotes a receptor-like sequence.

## Citation

If you use these materials, please cite the associated manuscript. Full bibliographic information and DOI will be added upon publication.

## License

Code and data are released under the [MIT License](LICENSE).
