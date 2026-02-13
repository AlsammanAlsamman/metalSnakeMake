#!/bin/bash
# Submit Step 03: LiftOver Genome Coordinates
# Process all datasets in parallel (7 SLURM jobs)

./submit.sh --snakefile rules/03_liftover.smk --cores 1 --jobs 7
