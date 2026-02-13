#!/bin/bash
# Submit Step 02: Filter Target SNPs
# Process all datasets in parallel (7 SLURM jobs)

./submit.sh --snakefile rules/02_filter_snps.smk --cores 1 --jobs 7
