#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: build_gene_neighborhood_database.sh --output-dir DIRECTORY --taxonkit PATH

Download NCBI gene-neighborhood source files and create gene_neighborhood.tsv.
This is a large, one-time operation. Download archives are retained in
DIRECTORY/downloads for reproducibility; intermediate processing files are
removed after a successful build.

EOF
}

output_dir=""; taxonkit=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir) output_dir=${2:?"--output-dir requires a directory"}; shift 2 ;;
        --taxonkit) taxonkit=${2:?"--taxonkit requires a path"}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Error: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n "$output_dir" && -n "$taxonkit" ]] || { usage >&2; exit 2; }
[[ -x "$taxonkit" ]] || { echo "Error: TaxonKit is not executable: $taxonkit" >&2; exit 1; }
for command in curl tar awk grep sort sed gzip; do command -v "$command" >/dev/null || { echo "Error: $command is required." >&2; exit 1; }; done

mkdir -p "$(dirname -- "$output_dir")"
output_dir=$(cd -- "$(dirname -- "$output_dir")" && pwd)/$(basename -- "$output_dir")
database="$output_dir/gene_neighborhood.tsv"
[[ -s "$database" ]] && { echo "Database already exists: $database"; exit 0; }
mkdir -p "$output_dir/downloads" "$output_dir/taxdump"
downloads="$output_dir/downloads"; taxdump="$output_dir/taxdump"
threads=$(nproc)

validate_tar() { tar -tzf "$1" >/dev/null; }
validate_gzip() { gzip -t "$1"; }

download() {
    local url=$1 destination=$2 validator=$3 temporary="${2}.part"

    if [[ -s "$destination" ]] && "$validator" "$destination" >/dev/null 2>&1; then
        echo "Verified existing download: $(basename -- "$destination")"
        return
    fi
    [[ -e "$destination" ]] && echo "Removing incomplete or corrupt download: $(basename -- "$destination")"
    rm -f "$destination" "$temporary"
    echo "Downloading: $(basename -- "$destination") (large NCBI source file)..."
    curl --fail --location --retry 5 --retry-all-errors --output "$temporary" "$url"
    "$validator" "$temporary" >/dev/null
    mv "$temporary" "$destination"
}

echo 'Checking one-time NCBI source downloads...'
download "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz" "$downloads/new_taxdump.tar.gz" validate_tar
download "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_neighbors.gz" "$downloads/gene_neighbors.gz" validate_gzip
download "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz" "$downloads/gene_info.gz" validate_gzip
download "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz" "$downloads/gene2accession.gz" validate_gzip

echo 'Preparing NCBI taxonomy files...'
tar -xzf "$downloads/new_taxdump.tar.gz" -C "$taxdump"
work=$(mktemp -d "$output_dir/.build.XXXXXX")
trap 'rm -rf "$work"' EXIT

if command -v pigz >/dev/null; then
    decompressor=(pigz -dk)
else
    decompressor=(gunzip -k)
fi
for archive in gene_neighbors gene_info gene2accession; do
    if [[ ! -s "$downloads/$archive" ]]; then
        echo "Extracting: $archive.gz"
        "${decompressor[@]}" "$downloads/$archive.gz"
    else
        echo "Verified existing extraction: $archive"
    fi
done

echo 'Preprocessing taxonomy, gene annotations, and protein accessions...'
awk '{print $1}' "$taxdump/taxidlineage.dmp" | "$taxonkit" filter --data-dir "$taxdump" -E Species -E Subspecies | \
    "$taxonkit" lineage --data-dir "$taxdump" -n -L | sed 's/\./_/g; s/ /_/g' > "$work/taxid2species" &
grep -E '#tax_id|protein-coding' "$downloads/gene_info" | sed 's/ /_/g' > "$work/gene_info_protein" &
awk '$6 != "-" { key=$1" "$2; values[key]=values[key]"|"$6 } END { for (key in values) print key, substr(values[key], 2) }' \
    "$downloads/gene2accession" \
    > "$work/gene2accession_merged" &
wait

echo 'Mapping species names to gene neighborhoods...'
awk 'NR==FNR {species[$1]=$2; next} {print $1, species[$1], $2, $5, $6, $7, $8, $3}' \
    "$work/taxid2species" "$downloads/gene_neighbors" > "$work/neighbors_v2"

echo 'Enriching neighborhoods with gene symbols and descriptions...'
awk 'NR==FNR {metadata[$2]=$3"\t"$9; next} {print $1, $2, $3, metadata[$3], $4, $5, $6, metadata[$9], $7, $8}' \
    "$work/gene_info_protein" "$work/neighbors_v2" > "$work/neighbors_v3"

echo 'Merging accessions and sorting the final database...'
{
    printf '#taxID\tspecies\tgeneID\tsymbol\tdescription\tstart\tend\torientation\tchromosome\tgenomic_accession.version\tprotein_accession.version\n'
    awk 'NR==FNR {accessions[$2]=$3; next} {print $0, accessions[$3]}' \
        "$work/gene2accession_merged" "$work/neighbors_v3" | \
        awk 'NF>10' | tail -n +2 | sort --parallel="$threads" -k10,10 -k6,6n | sed 's/ /\t/g'
} > "$database"

echo "Created database: $database"
