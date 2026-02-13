#!/bin/bash
# Submit Step 7: Enrich Meta-Analysis with Study-Specific EAF
# Process all METAL combinations in parallel

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR" || exit 1

# Default parameters
DRY_RUN=""
CORES=1
JOBS=10

# Parse command line arguments
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

echo "Submitting Step 7: Enrich Meta-Analysis with Study-Specific EAF"
echo "Cores: $CORES"
echo "Jobs: $JOBS"

if [[ -n "$DRY_RUN" ]]; then
    echo "Dry run mode - no jobs will be submitted"
fi

# Submit to Snakemake
./submit.sh \
    --snakefile rules/07_enrich_eaf.smk \
    --cores "$CORES" \
    --jobs "$JOBS" \
    $DRY_RUN

echo "Step 7 submission complete"
