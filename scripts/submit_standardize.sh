#!/bin/bash
# Submit script for Step 0: Standardization
# Processes all datasets in parallel (one job per dataset via SLURM cluster)

# Use main submit.sh to handle cluster submission
# --jobs 7 allows up to 7 parallel jobs (one per dataset)
./submit.sh --snakefile rules/00_standardize.smk --cores 2 --jobs 7 "$@"

echo ""
echo "Standardization submitted to cluster with 7 parallel jobs."
echo "Check logs/ for job outputs."
echo "Each dataset produces:"
echo "  - results/00_standardized/{dataset}.tsv (standardized data)"
echo "  - results/00_standardized/{dataset}.done (completion marker)"
