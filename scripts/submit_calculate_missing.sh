#!/bin/bash
# Submit Step 01: Calculate Missing Columns
# Process all datasets in parallel (7 SLURM jobs)

./submit.sh --snakefile rules/01_calculate_missing.smk --cores 2 --jobs 7
