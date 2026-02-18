#!/usr/bin/env python3
"""
Step 02: Filter Target SNPs
Filter variants to target SNP list if provided, otherwise copy unchanged.
"""

import sys
sys.path.append("utils")

from bioconfigme import (
    get_results_dir,
    get_dataset_list,
    get_analysis_value,
    get_software_module
)

# Get configuration
RESULTS_DIR = get_results_dir()
DATASETS = get_dataset_list()

# Get target SNP list (optional)
TARGET_SNP_LIST = get_analysis_value(['target_snp_list'], required=False, default="")

# Get R module
R_MODULE = get_software_module('r')

# Rule targets
rule all:
    input:
        expand("results/02_filtered/{dataset}.done", dataset=DATASETS)


rule filter_target_snps:
    """
    Filter variants to target SNP list.
    If target SNP list not provided or doesn't exist, copy all variants unchanged.
    """
    input:
        tsv = "results/01_calculated/{dataset}.tsv",
        done = "results/01_calculated/{dataset}.done"
    output:
        tsv = "results/02_filtered/{dataset}.tsv",
        done = "results/02_filtered/{dataset}.done"
    log:
        "results/log/02_filter/{dataset}.log"
    params:
        dataset = "{dataset}",
        snp_list = TARGET_SNP_LIST,
        r_module = R_MODULE
    resources:
        mem_mb = 16000,
        time = "00:15:00",
        cores = 1
    shell:
        """
        set +u
        module load {params.r_module}
        set -u

        SNP_LIST_ARG=""
        if [[ -n "{params.snp_list}" ]]; then
            SNP_LIST_ARG="--snp_list {params.snp_list}"
        fi
        
        Rscript scripts/filter_target_snps.R \
            --dataset {params.dataset} \
            --input {input.tsv} \
            --output {output.tsv} \
            --done {output.done} \
            --log {log} \
            $SNP_LIST_ARG \
            2>&1 | tee {log}
        """
