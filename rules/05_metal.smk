#!/usr/bin/env python3
"""
Step 05: METAL Meta-Analysis
Run inverse-variance weighted meta-analysis using METAL.
"""

import sys
sys.path.append("utils")

from bioconfigme import (
    get_results_dir,
    get_metal_config,
    get_metal_combinations,
    get_metal_binary,
    get_software_module
)

# Get configuration
RESULTS_DIR = get_results_dir()
METAL_CONFIG = get_metal_config()
COMBINATIONS = get_metal_combinations()

# Get METAL binary path and module
METAL_BIN = get_metal_binary()
METAL_MODULE = get_software_module('metal')

# Get METAL scheme
SCHEME = METAL_CONFIG.get('scheme', 'STDERR')
USE_EAF = METAL_CONFIG.get('use_eaf', True)

# Rule targets
rule all:
    input:
        expand("results/05_metal/{combination}/{combination}.done", 
               combination=COMBINATIONS.keys())


rule metal:
    """
    Run METAL meta-analysis for a study combination.
    Generates METAL configuration file and executes analysis.
    """
    input:
        tsv = lambda wildcards: expand("results/04_metal_ready/{dataset}.tsv", 
                                       dataset=COMBINATIONS[wildcards.combination]),
        done = lambda wildcards: expand("results/04_metal_ready/{dataset}.done",
                                        dataset=COMBINATIONS[wildcards.combination])
    output:
        done = "results/05_metal/{combination}/{combination}.done",
        tbl = "results/05_metal/{combination}/{combination}1.tbl",
        config = "results/05_metal/{combination}/{combination}.metal"
    log:
        "results/log/05_metal/{combination}.log"
    params:
        combination = "{combination}",
        datasets = lambda wildcards: COMBINATIONS[wildcards.combination],
        input_dir = "results/04_metal_ready",
        output_dir = lambda wildcards: f"results/05_metal/{wildcards.combination}",
        scheme = SCHEME,
        use_eaf = USE_EAF,
        metal_bin = METAL_BIN,
        metal_module = METAL_MODULE
    resources:
        mem_mb = 64000,
        time = "01:00:00",
        cores = 1
    shell:
        """
        set +u
        module load {params.metal_module}
        set -u
        
        python3 scripts/run_metal.py \
            --combination {params.combination} \
            --datasets {params.datasets} \
            --input_dir {params.input_dir} \
            --output_dir {params.output_dir} \
            --done {output.done} \
            --log {log} \
            --scheme {params.scheme} \
            --use_eaf {params.use_eaf} \
            --metal_bin {params.metal_bin} \
        2>&1 | tee -a {log}
        """
