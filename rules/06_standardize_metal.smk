#!/usr/bin/env python3
"""
Step 06: Standardize METAL Output
Standardize METAL meta-analysis results, split MarkerName, create Manhattan plots.
"""

import sys
sys.path.append("utils")

from bioconfigme import (
    get_results_dir,
    get_metal_combinations,
    get_analysis_value,
    get_software_module
)

# Get configuration
RESULTS_DIR = get_results_dir()
COMBINATIONS = get_metal_combinations()

# Get markername format
MARKERNAME_FORMAT = get_analysis_value(['markername_format'], default='chrom:pos:snpid')

# Get R module
R_MODULE = get_software_module('r')

# Rule targets
rule all:
    input:
        expand("results/06_metal_standardized/{combination}/{combination}.done",
               combination=COMBINATIONS.keys())


rule standardize_metal_output:
    """
    Standardize METAL meta-analysis output.
    - Rename columns to standard format
    - Split MarkerName back to constituent columns
    - Filter SNPs with p < 5e-3 and gzip
    - Generate Manhattan plots (PDF and PNG)
    - Sort by chromosome and position
    """
    input:
        tbl = "results/05_metal/{combination}/{combination}1.tbl",
        done = "results/05_metal/{combination}/{combination}.done"
    output:
        tsv = "results/06_metal_standardized/{combination}/{combination}_standardized.tsv",
        filtered = "results/06_metal_standardized/{combination}/{combination}_sig.tsv.gz",
        excel = "results/06_metal_standardized/{combination}/{combination}_sig5e5.xlsx",
        plot_pdf = "results/06_metal_standardized/{combination}/{combination}_manhattan.pdf",
        plot_png = "results/06_metal_standardized/{combination}/{combination}_manhattan.png",
        done = "results/06_metal_standardized/{combination}/{combination}.done"
    log:
        "results/log/06_standardize_metal/{combination}.log"
    params:
        combination = "{combination}",
        studies = lambda wildcards: ",".join(COMBINATIONS[wildcards.combination]),
        markername_format = MARKERNAME_FORMAT,
        r_module = R_MODULE
    resources:
        mem_mb = 32000,
        time = "01:00:00",
        cores = 1
    shell:
        """
        set +u
        module load {params.r_module}
        set -u
        
        Rscript scripts/standardize_metal_output.R \
            --input {input.tbl} \
            --output {output.tsv} \
            --output_filtered {output.filtered} \
            --output_excel {output.excel} \
            --plot_pdf {output.plot_pdf} \
            --plot_png {output.plot_png} \
            --combination {params.combination} \
            --studies "{params.studies}" \
            --markername_format {params.markername_format} \
            --log {log} \
        2>&1 | tee -a {log}
        
        echo "Standardization completed for {params.combination}" > {output.done}
        echo "Input: {input.tbl}" >> {output.done}
        echo "Output: {output.tsv}" >> {output.done}
        echo "Filtered (p<5e-3): {output.filtered}" >> {output.done}
        echo "Excel (p<5e-5): {output.excel}" >> {output.done}
        echo "Plots: {output.plot_pdf}, {output.plot_png}" >> {output.done}
        """
