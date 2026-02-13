#!/bin/bash
# Submit Step 04: Prepare Files for METAL
# Process all datasets in parallel (7 SLURM jobs)

./submit.sh --snakefile rules/04_prepare_metal.smk --cores 1 --jobs 7
