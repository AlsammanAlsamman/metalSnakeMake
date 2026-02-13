"""
Snakemake rule for standardizing column names and extracting required columns
Step 0: Standardization - Extract and rename columns to standard format
"""

import sys
sys.path.append("utils")
from bioconfigme import get_dataset_list, get_dataset_config

# Get all datasets from configuration
DATASETS = get_dataset_list()

# Target rule to process all datasets
rule all:
    input:
        expand("results/00_standardized/{dataset}.done", dataset=DATASETS)

# Standardization rule - one job per dataset
rule standardize_columns:
    input:
        data=lambda wildcards: get_dataset_config(wildcards.dataset)['file']
    output:
        tsv="results/00_standardized/{dataset}.tsv",
        done="results/00_standardized/{dataset}.done"
    log:
        "results/log/00_standardize/{dataset}.log"
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2
    shell:
        """
        python scripts/standardize_columns.py \
            --dataset {wildcards.dataset} \
            --output {output.tsv} \
            --done {output.done} \
            --log {log} 2>&1 | tee {log}
        """
