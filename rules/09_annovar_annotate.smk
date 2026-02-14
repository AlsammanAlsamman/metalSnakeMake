#!/usr/bin/env python3
"""
Step 09: Annotate Meta-Analysis with ANNOVAR
Add functional gene-based annotations using ANNOVAR.
Uses target_build from analysis.yml (hg38 or hg19).
"""

import sys
sys.path.append("utils")

from bioconfigme import (
    get_results_dir,
    get_metal_combinations,
    get_software_module,
    get_analysis_value
)

# Get configuration
RESULTS_DIR = get_results_dir()
COMBINATIONS = get_metal_combinations()

# Get target build (hg38 or hg19)
TARGET_BUILD = get_analysis_value(['target_build'])

# Get ANNOVAR module (perl is usually available by default)
ANNOVAR_MODULE = get_software_module('annovar')

# ANNOVAR enrichment parameters from config
ANNOVAR_HUMANDB = get_analysis_value(['enrichment', 'annovar', 'humandb'])
ANNOVAR_PROTOCOLS = ",".join(get_analysis_value(['enrichment', 'annovar', 'protocols']))
ANNOVAR_OPERATIONS = ",".join(get_analysis_value(['enrichment', 'annovar', 'operations']))
ANNOVAR_EXCEL_PVAL = get_analysis_value(['enrichment', 'annovar', 'excel_pvalue_threshold'])

# Get R module
R_MODULE = get_software_module('r')

# Rule targets
rule all:
    input:
        expand("results/09_annovar_annotated/{combination}/{combination}.done",
               combination=COMBINATIONS.keys())


rule annovar_annotate:
    """
    Annotate meta-analysis results with ANNOVAR gene-based annotations.
    - Convert to ANNOVAR input format (chr, start, end, ref, alt)
    - Run ANNOVAR table_annovar.pl with refGene protocol
    - Merge annotations back to enriched file
    - Generate TSV with all SNPs + annotations
    - Generate Excel with SNPs p < 5e-3
    """
    input:
        meta = "results/07_enriched/{combination}/{combination}_enriched.tsv",
        done = "results/07_enriched/{combination}/{combination}.done"
    output:
        tsv = "results/09_annovar_annotated/{combination}/{combination}_annotated.tsv",
        excel = "results/09_annovar_annotated/{combination}/{combination}_annotated.xlsx",
        done = "results/09_annovar_annotated/{combination}/{combination}.done"
    log:
        "results/log/09_annovar_annotate/{combination}.log"
    params:
        combination = "{combination}",
        annovar_module = ANNOVAR_MODULE,
        r_module = R_MODULE,
        buildver = TARGET_BUILD,
        humandb = ANNOVAR_HUMANDB,
        protocols = ANNOVAR_PROTOCOLS,
        operations = ANNOVAR_OPERATIONS,
        excel_pval = ANNOVAR_EXCEL_PVAL,
        outdir = "results/09_annovar_annotated/{combination}"
    resources:
        mem_mb = 16000,
        time = "01:00:00",
        cores = 1
    shell:
        """
        set +u
        module load {params.annovar_module}
        module load {params.r_module}
        set -u
        
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        
        # Run annotation and merging script
        bash scripts/annotate_annovar.sh \
            --input {input.meta} \
            --output {output.tsv} \
            --output_excel {output.excel} \
            --combination {params.combination} \
            --buildver {params.buildver} \
            --humandb {params.humandb} \
            --protocols {params.protocols} \
            --operations {params.operations} \
            --excel_pval {params.excel_pval} \
            --outdir {params.outdir} \
            --log {log} \
        2>&1 | tee -a {log}
        
        echo "ANNOVAR annotation completed for {params.combination}" > {output.done}
        echo "Input: {input.meta}" >> {output.done}
        echo "Output: {output.tsv}" >> {output.done}
        echo "Excel: {output.excel}" >> {output.done}
        """
