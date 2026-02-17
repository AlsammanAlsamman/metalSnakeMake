#!/usr/bin/env bash
set -euo pipefail

# Split each SNPdb build VCF into per-chromosome VCFs in parallel.
# Default input:  resources/SNPdb/*.vcf.gz
# Default output: resources/SNPdb/by_chr/<BUILD>/<CHROM>.vcf.gz
#
# Usage:
#   bash utility_scripts/split_snpdb_by_chr.sh
#   bash utility_scripts/split_snpdb_by_chr.sh --workers 16 --force
#   bash utility_scripts/split_snpdb_by_chr.sh --input-dir resources/SNPdb --output-root resources/SNPdb/by_chr

INPUT_DIR="resources/SNPdb"
OUTPUT_ROOT="resources/SNPdb/by_chr"
WORKERS="$(nproc 2>/dev/null || echo 8)"
FORCE=0

print_help() {
    cat <<'EOF'
Split SNPdb VCFs by chromosome.

Options:
  --input-dir <path>     Directory containing build VCFs (default: resources/SNPdb)
  --output-root <path>   Root output dir (default: resources/SNPdb/by_chr)
  --workers <int>        Parallel chromosomes per build (default: nproc or 8)
  --force                Recreate existing chromosome outputs
  -h, --help             Show this help

Expected input files:
  <INPUT_DIR>/GRCh37.vcf.gz
  <INPUT_DIR>/GRCh38.vcf.gz
  (or any *.vcf.gz build files)

Output layout:
  <OUTPUT_ROOT>/<BUILD>/<CHROM>.vcf.gz
  <OUTPUT_ROOT>/<BUILD>/<CHROM>.vcf.gz.tbi
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output-root)
            OUTPUT_ROOT="$2"
            shift 2
            ;;
        --workers)
            WORKERS="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        -h|--help)
            print_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            print_help
            exit 1
            ;;
    esac
done

if ! command -v bcftools >/dev/null 2>&1; then
    echo "Error: bcftools is required but not found in PATH." >&2
    exit 1
fi

if ! command -v tabix >/dev/null 2>&1; then
    echo "Error: tabix is required but not found in PATH." >&2
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: input directory not found: $INPUT_DIR" >&2
    exit 1
fi

if ! [[ "$WORKERS" =~ ^[0-9]+$ ]] || [[ "$WORKERS" -lt 1 ]]; then
    echo "Error: --workers must be a positive integer." >&2
    exit 1
fi

mkdir -p "$OUTPUT_ROOT"

split_one_chrom() {
    local vcf="$1"
    local chrom="$2"
    local out_dir="$3"
    local force="$4"

    local out_vcf="$out_dir/${chrom}.vcf.gz"

    if [[ "$force" -eq 0 && -s "$out_vcf" && -s "${out_vcf}.tbi" ]]; then
        echo "  [skip] ${chrom} already exists"
        return 0
    fi

    bcftools view -r "$chrom" -Oz -o "$out_vcf" "$vcf"
    tabix -f -p vcf "$out_vcf"
    echo "  [done] ${chrom}"
}

process_build_vcf() {
    local vcf="$1"
    local build
    build="$(basename "$vcf" .vcf.gz)"

    local build_out_dir="$OUTPUT_ROOT/$build"
    mkdir -p "$build_out_dir"

    echo "Processing build: $build"
    echo "  input:  $vcf"
    echo "  output: $build_out_dir"

    if [[ ! -s "${vcf}.tbi" ]]; then
        echo "  Index missing; creating: ${vcf}.tbi"
        tabix -f -p vcf "$vcf"
    fi

    local chroms=()

    # Preferred: use index summary (fast), but this can fail on older index metadata.
    mapfile -t chroms < <(bcftools index -s "$vcf" 2>/dev/null | awk '{print $1}')

    # Fallback 1: read sequence names from index.
    if [[ ${#chroms[@]} -eq 0 ]]; then
        mapfile -t chroms < <(tabix -l "$vcf" 2>/dev/null)
    fi

    # Fallback 2: parse contigs from VCF header.
    if [[ ${#chroms[@]} -eq 0 ]]; then
        mapfile -t chroms < <(
            bcftools view -h "$vcf" \
            | awk -F'[=,>]' '/^##contig=<ID=/{print $3}'
        )
    fi

    if [[ ${#chroms[@]} -eq 0 ]]; then
        echo "Error: could not detect chromosomes for $vcf (index + header methods failed)" >&2
        return 1
    fi

    echo "  chromosomes detected: ${#chroms[@]}"

    local pids=()
    local active=0

    for chrom in "${chroms[@]}"; do
        split_one_chrom "$vcf" "$chrom" "$build_out_dir" "$FORCE" &
        pids+=("$!")
        ((active+=1))

        if [[ "$active" -ge "$WORKERS" ]]; then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
            ((active-=1))
        fi
    done

    for pid in "${pids[@]}"; do
        wait "$pid"
    done

    echo "Finished build: $build"
    echo
}

shopt -s nullglob
build_vcfs=("$INPUT_DIR"/*.vcf.gz)
shopt -u nullglob

if [[ ${#build_vcfs[@]} -eq 0 ]]; then
    echo "Error: no .vcf.gz files found in $INPUT_DIR" >&2
    exit 1
fi

echo "Found ${#build_vcfs[@]} build VCF(s) in $INPUT_DIR"
echo "Workers per build: $WORKERS"
echo

for vcf in "${build_vcfs[@]}"; do
    process_build_vcf "$vcf"
done

echo "All builds completed successfully."
