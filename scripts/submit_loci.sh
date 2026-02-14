#!/bin/bash
# Submit wrapper for Step 10: Loci Identification
# This script submits jobs for all configured study combinations

# Check if dry-run flag is passed
DRY_RUN=""
if [[ "$1" == "--dry-run" ]]; then
    DRY_RUN="--dry-run"
    echo "Running in dry-run mode..."
fi

# Submit all loci identification jobs
./submit.sh --snakefile rules/10_loci_identification.smk --cores 8 --jobs 10 $DRY_RUN
