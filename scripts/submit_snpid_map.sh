#!/bin/bash
# Submit Step 00b: Map snpid (rsID) from SNPdb
# Process all datasets in parallel

./submit.sh --snakefile rules/00b_map_snpid.smk --cores 2 --jobs 7
