#!/usr/bin/env python3
"""
Step 04: Prepare Files for METAL
Add MarkerName column with configurable format for METAL meta-analysis.
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

# Get markername format
MARKERNAME_FORMAT = get_analysis_value(['markername_format'], default='chrom:pos:snpid')

# Get R module
R_MODULE = get_software_module('r')

# Rule targets
rule all:
    input:
        expand("results/04_metal_ready/{dataset}.done", dataset=DATASETS)


rule prepare_metal:
    """
    Prepare files for METAL meta-analysis.
    Adds MarkerName column based on configurable format.
    Default format: chrom:pos:snpid
    """
    input:
        tsv = "results/03_lifted/{dataset}.tsv",
        done = "results/03_lifted/{dataset}.done"
    output:
        tsv = "results/04_metal_ready/{dataset}.tsv",
        done = "results/04_metal_ready/{dataset}.done"
    log:
        "results/log/04_prepare_metal/{dataset}.log"
    params:
        dataset = "{dataset}",
        markername_format = MARKERNAME_FORMAT,
        r_module = R_MODULE
    resources:
        mem_mb = 16000,
        time = "00:30:00",
        cores = 1
    shell:
        """
        set +u
        module load {params.r_module}
        set -u
        
        Rscript scripts/prepare_metal.R \
            --dataset {params.dataset} \
            --input {input.tsv} \
            --output {output.tsv} \
            --done {output.done} \
            --log {log} \
            --markername_format {params.markername_format} \
            2>&1 | tee {log}
        """
