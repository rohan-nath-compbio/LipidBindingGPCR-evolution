# KEGG ortholog tables

This directory contains the KEGG-derived ortholog tables used in the comparative analysis. Each TSV file is named for the query gene symbol used to retrieve its KEGG Orthology (KO) group.

## Columns

| Column | Description |
| --- | --- |
| `Phylum` | Broad animal group assigned from the KEGG organism lineage. |
| `Species` | KEGG organism name, normalized with underscores. |
| `KEGG_GeneID` | Numeric component of the organism-specific KEGG gene identifier. |
| `NCBI_ProteinID` | NCBI protein accession returned by KEGG's `conv/ncbi-proteinid` endpoint; blank where no mapping was available. |

## Scope and provenance

The retrieval workflow keeps members assigned by KEGG to Vertebrata, Tunicata, Cephalochordata, Hemichordata, Echinodermata, Annelida, Mollusca, Arthropoda, Cnidaria, or Placozoa. The tables are deposited snapshots; KEGG is continuously updated, so rerunning the script may produce different results.

Regenerate a table from the repository root with, for example:

```bash
python3 scripts/kegg/fetch_kegg_orthologs.py DAGL data/kegg/dagl.tsv
```

The script uses the public [KEGG REST API](https://rest.kegg.jp/) and requires only the Python standard library.
