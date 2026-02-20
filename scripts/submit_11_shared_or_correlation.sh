#!/bin/bash
# Submit Step 11: Shared SNP OR Correlation and Harmonization

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR" || exit 1

DRY_RUN=""
CORES=1
JOBS=10

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN="--dry-run"
            shift
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        --jobs)
            JOBS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--dry-run] [--cores N] [--jobs N]"
            exit 1
            ;;
    esac
done

./submit.sh \
    --snakefile rules/11_shared_or_correlation.smk \
    --cores "$CORES" \
    --jobs "$JOBS" \
    $DRY_RUN
