#!/bin/bash
# Submit script for Step 0c: API remap of previously unmapped SNPs from Step 00b
# Processes all datasets in parallel (one job per dataset)

./submit.sh --snakefile rules/00c_api_map_unmapped.smk --cores 1 --jobs 7 "$@"

echo ""
echo "API remap submitted to cluster with up to 7 parallel jobs."
echo "Outputs per dataset:"
echo "  - results/snpdmapped/{dataset}_mapped.tsv"
echo "  - results/snpdmapped/{dataset}_mapped.done"
echo "  - results/snpdmapped/summary/{dataset}_api_mapping_summary.tsv"
