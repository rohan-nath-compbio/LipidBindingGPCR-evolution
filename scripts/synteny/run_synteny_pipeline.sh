#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_synteny_pipeline.sh FILE.acc

Create taxonomic lineage annotations and gene-neighborhood output from a file
containing one NCBI protein accession per line.

Arguments:
  FILE.acc                   Input accession file (one accession per line).

Environment variables (optional):
  NCBI_EMAIL                 Contact email for NCBI E-utilities. If not set,
                             the script asks for it interactively.
  NCBI_API_KEY               NCBI API key for faster requests.

Stages:
  1. Install TaxonKit and build the NCBI gene-neighborhood database if missing.
  2. Resolve protein accessions to NCBI taxonomy IDs.
  3. Obtain TaxonKit lineages and write <name>.acc.taxid.lineage.
  4. Extract, clean, and arrange gene neighborhoods.

Examples:
  scripts/synteny/run_synteny_pipeline.sh data/synteny/example.acc

  NCBI_EMAIL=you@example.org NCBI_API_KEY=... \
    scripts/synteny/run_synteny_pipeline.sh my_receptor.acc
EOF
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then usage; exit 0; fi
[[ $# -eq 1 ]] || { echo "Error: provide exactly one .acc file." >&2; usage >&2; exit 2; }
input=$1; email=${NCBI_EMAIL:-}; api_key=${NCBI_API_KEY:-}; resources_dir=""; output_dir=""; setup=true
if [[ -z "$email" ]]; then
    [[ -t 0 ]] || { echo "Error: set NCBI_EMAIL when running non-interactively." >&2; exit 2; }
    read -r -p "NCBI contact email: " email
fi
[[ -n "$email" ]] || { echo "Error: an NCBI contact email is required." >&2; exit 2; }
[[ -f "$input" ]] || { echo "Error: input file not found: $input" >&2; exit 1; }
[[ "$input" == *.acc ]] || { echo "Error: input file must use the .acc extension." >&2; exit 2; }

launcher_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$launcher_dir/../.." && pwd)
script_dir="$launcher_dir"
input=$(cd -- "$(dirname -- "$input")" && pwd)/$(basename -- "$input")
name=$(basename -- "$input" .acc)
resources_dir=${resources_dir:-"$repo_root/resources/synteny"}
output_dir=${output_dir:-"$repo_root/results/synteny/$name"}
taxonkit="$repo_root/tools/bin/taxonkit"
database="$resources_dir/gene-neighborhood/gene_neighborhood.tsv"
taxdump="$resources_dir/gene-neighborhood/taxdump"

if [[ ! -x "$taxonkit" ]]; then
    $setup || { echo "Error: TaxonKit is missing: $taxonkit" >&2; exit 1; }
    "$script_dir/install_taxonkit.sh" --bin-dir "$repo_root/tools/bin"
fi
if [[ ! -s "$database" ]]; then
    $setup || { echo "Error: database is missing: $database" >&2; exit 1; }
    echo 'Gene-neighborhood database is absent; starting its one-time NCBI build.'
    echo 'This expands and externally sorts large source files. Progress is reported per stage.'
    "$script_dir/build_gene_neighborhood_database.sh" --output-dir "$resources_dir/gene-neighborhood" --taxonkit "$taxonkit"
fi
[[ -d "$taxdump" ]] || { echo "Error: TaxonKit taxonomy data are missing: $taxdump" >&2; exit 1; }

mkdir -p "$output_dir"
work_acc="$output_dir/$name.acc"
taxid_file="$output_dir/$name.acc.taxid"
lineage_file="$output_dir/$name.acc.taxid.lineage"
lineage_map="$output_dir/$name.taxid.lineage.tsv"
cp "$input" "$work_acc"

echo '[1/3] Resolving protein accessions to taxonomy IDs...'
resolver=(python3 "$script_dir/resolve_accession_taxids.py" "$work_acc" "$taxid_file" --email "$email")
[[ -n "$api_key" ]] && resolver+=(--api-key "$api_key")
"${resolver[@]}"
[[ -s "$taxid_file" ]] || { echo "Error: no taxonomy IDs were resolved." >&2; exit 1; }

echo '[2/3] Fetching TaxonKit lineages...'
awk '{print $2}' "$taxid_file" | "$taxonkit" lineage --data-dir "$taxdump" | sed 's/ /_/g' > "$lineage_map"
awk -F'\t' 'BEGIN {OFS="\t"} NR == FNR {lineage[$1] = $2; next} $2 in lineage {print $1, $2, lineage[$2]}' \
    "$lineage_map" "$taxid_file" > "$lineage_file"
[[ -s "$lineage_file" ]] || { echo "Error: no lineage annotations were produced." >&2; exit 1; }

echo '[3/3] Extracting gene neighborhoods...'
(
    cd "$output_dir"
    GENE_NEIGHBORHOOD_DB="$database" "$script_dir/extract_gene_neighborhoods.sh"
)

arranged_file="$output_dir/$name.arranged.neighbor"
if [[ -s "$arranged_file" ]]; then
    find "$output_dir" -maxdepth 1 -type f -name "$name.*" ! -name "$name.arranged.neighbor" -delete
    echo "Completed. Arranged neighborhoods: $arranged_file"
else
    echo "Warning: no arranged-neighborhood report was created; retaining intermediate files in $output_dir." >&2
fi
exit 0
