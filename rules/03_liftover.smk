#!/usr/bin/env python3
"""
Step 03: LiftOver Genome Coordinates
Convert genome builds to target build using UCSC liftOver.
Datasets already in target build are copied unchanged.
"""

import sys
sys.path.append("utils")

from bioconfigme import (
    get_results_dir,
    get_dataset_list,
    get_dataset_config,
    get_target_build,
    get_liftover_chain,
    get_liftover_binary,
    get_software_module
)

# Get configuration
RESULTS_DIR = get_results_dir()
DATASETS = get_dataset_list()
TARGET_BUILD = get_target_build()

# Get liftOver module and binary
LIFTOVER_MODULE = get_software_module('liftover')
LIFTOVER_BIN = get_liftover_binary()

# Get R module (for Rscript)
R_MODULE = get_software_module('r')

# Rule targets
rule all:
    input:
        expand("results/03_lifted/{dataset}.done", dataset=DATASETS)


rule liftover_build:
    """
    LiftOver genome coordinates to target build.
    If source build matches target build, copy file unchanged.
    Otherwise, use UCSC liftOver to convert coordinates.
    Unmapped variants are removed from output.
    """
    input:
        tsv = "results/02_filtered/{dataset}.tsv",
        done = "results/02_filtered/{dataset}.done"
    output:
        tsv = "results/03_lifted/{dataset}.tsv",
        done = "results/03_lifted/{dataset}.done"
    log:
        "results/log/03_liftover/{dataset}.log"
    params:
        dataset = "{dataset}",
        source_build = lambda wildcards: get_dataset_config(wildcards.dataset).get('build', 'hg19'),
        target_build = TARGET_BUILD,
        chain_file = lambda wildcards: get_liftover_chain(
            get_dataset_config(wildcards.dataset).get('build', 'hg19'),
            TARGET_BUILD
        ) if get_dataset_config(wildcards.dataset).get('build', 'hg19') != TARGET_BUILD else "NA",
        liftover_bin = LIFTOVER_BIN,
        liftover_module = LIFTOVER_MODULE,
        r_module = R_MODULE
    resources:
        mem_mb = 16000,
        time = "01:00:00",
        cores = 1
    shell:
        """
        set +u
        module load {params.r_module}
        module load {params.liftover_module}
        set -u
        
        Rscript scripts/liftover_build.R \
            --dataset {params.dataset} \
            --input {input.tsv} \
            --output {output.tsv} \
            --done {output.done} \
            --log {log} \
            --source_build {params.source_build} \
            --target_build {params.target_build} \
            --chain_file {params.chain_file} \
            --liftover_bin {params.liftover_bin} \
            2>&1 | tee {log}
        """
