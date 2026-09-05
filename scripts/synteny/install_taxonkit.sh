#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: install_taxonkit.sh [--bin-dir DIRECTORY]

Create (when needed) the Conda environment named "taxonkit" with:
  conda create --name taxonkit -c bioconda taxonkit

A project-local launcher is then written to DIRECTORY/taxonkit so the main
pipeline can invoke TaxonKit consistently.
EOF
}

bin_dir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bin-dir) bin_dir=${2:?"--bin-dir requires a directory"}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Error: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

command -v conda >/dev/null || { echo "Error: Conda is required to install TaxonKit." >&2; exit 1; }
if ! conda env list | awk 'NR > 2 {print $1}' | grep -Fxq taxonkit; then
    echo "Creating Conda environment: taxonkit"
    conda create --name taxonkit -c bioconda taxonkit -y
fi
conda run --no-capture-output --name taxonkit taxonkit version

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
bin_dir=${bin_dir:-"$repo_root/tools/bin"}
mkdir -p "$bin_dir"
cat > "$bin_dir/taxonkit" <<'EOF'
#!/usr/bin/env bash
exec conda run --no-capture-output --name taxonkit taxonkit "$@"
EOF
chmod 0755 "$bin_dir/taxonkit"
echo "TaxonKit launcher created: $bin_dir/taxonkit"
