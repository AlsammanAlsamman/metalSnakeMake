#!/usr/bin/env python3
"""
Step 01: Calculate Missing Essential Columns
Calculate missing columns required for meta-analysis (SE, EAF, etc.)
Currently implements: SE calculation from beta, eaf, n_cases, n_controls
"""

import sys
sys.path.append("utils")

from bioconfigme import (
    get_results_dir,
    get_dataset_list,
    get_dataset_config,
    get_software_module
)

# Get configuration
RESULTS_DIR = get_results_dir()
DATASETS = get_dataset_list()

# Get R module
R_MODULE = get_software_module('r')

# Rule targets
rule all:
    input:
        expand("results/01_calculated/{dataset}.done", dataset=DATASETS)


rule calculate_missing_columns:
    """
    Calculate missing essential columns (e.g., SE from beta, eaf, sample sizes).
    If all columns present, copy file unchanged.
    """
    input:
        tsv = "results/00_standardized/{dataset}.tsv",
        done = "results/00_standardized/{dataset}.done"
    output:
        tsv = "results/01_calculated/{dataset}.tsv",
        done = "results/01_calculated/{dataset}.done"
    log:
        "results/log/01_calculate/{dataset}.log"
    params:
        dataset = "{dataset}",
        n_cases = lambda wildcards: get_dataset_config(wildcards.dataset).get('n_cases', 0),
        n_controls = lambda wildcards: get_dataset_config(wildcards.dataset).get('n_controls', 0),
        r_module = R_MODULE
    resources:
        mem_mb = 32000,
        time = "00:30:00",
        cores = 2
    shell:
        """
        set +u
        module load {params.r_module}
        set -u
        
        Rscript scripts/calculate_se.R \
            --dataset {params.dataset} \
            --input {input.tsv} \
            --output {output.tsv} \
            --done {output.done} \
            --log {log} \
            --n_cases {params.n_cases} \
            --n_controls {params.n_controls} \
            2>&1 | tee {log}
        """
