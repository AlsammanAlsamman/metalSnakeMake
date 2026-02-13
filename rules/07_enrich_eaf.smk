#!/usr/bin/env python3
"""
Step 07: Enrich Meta-Analysis with Study-Specific EAF
Add study-specific effect allele frequencies to meta-analysis results.
Align EAF values to meta effect allele (flip to 1-EAF when needed).
Calculate meta EAF as mean of aligned study EAFs.
"""

import sys
sys.path.append("utils")

from bioconfigme import (
    get_results_dir,
    get_metal_combinations,
    get_software_module
)

# Get configuration
RESULTS_DIR = get_results_dir()
COMBINATIONS = get_metal_combinations()

# Get R module
R_MODULE = get_software_module('r')

# Rule targets
rule all:
    input:
        expand("results/07_enriched/{combination}/{combination}.done",
               combination=COMBINATIONS.keys())


rule enrich_meta_eaf:
    """
    Enrich meta-analysis results with study-specific EAF values.
    - Merge Step 6 standardized output with Step 4 metal_ready files
    - Align study EAF to meta effect allele (flip to 1-EAF when needed)
    - Add eaf_{STUDY} columns for each study
    - Calculate meta eaf as mean of aligned study EAFs (excluding NAs)
    - Generate Excel output with formatting
    """
    input:
        meta = "results/06_metal_standardized/{combination}/{combination}_standardized.tsv",
        done = "results/06_metal_standardized/{combination}/{combination}.done"
    output:
        tsv = "results/07_enriched/{combination}/{combination}_enriched.tsv",
        excel = "results/07_enriched/{combination}/{combination}_enriched.xlsx",
        done = "results/07_enriched/{combination}/{combination}.done"
    log:
        "results/log/07_enrich_eaf/{combination}.log"
    params:
        combination = "{combination}",
        studies = lambda wildcards: ",".join(COMBINATIONS[wildcards.combination]),
        r_module = R_MODULE
    resources:
        mem_mb = 64000,
        time = "01:00:00",
        cores = 1
    shell:
        """
        set +u
        module load {params.r_module}
        set -u
        
        Rscript scripts/enrich_meta_eaf.R \
            --meta {input.meta} \
            --output {output.tsv} \
            --output_excel {output.excel} \
            --combination {params.combination} \
            --studies "{params.studies}" \
            --log {log} \
        2>&1 | tee -a {log}
        
        echo "EAF enrichment completed for {params.combination}" > {output.done}
        echo "Input: {input.meta}" >> {output.done}
        echo "Output: {output.tsv}" >> {output.done}
        echo "Excel: {output.excel}" >> {output.done}
        echo "Studies: {params.studies}" >> {output.done}
        """
