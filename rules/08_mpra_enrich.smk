#!/usr/bin/env python3
"""
Step 08: Enrich ANNOVAR-Annotated Meta-Analysis with MPRA Functional Data
Match ANNOVAR-annotated SNPs with MPRA database by chr+pos.
Aggregate multiple MPRA records per SNP.
Filter MPRA for fdr < 0.05 before matching.
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

# Get R module
R_MODULE = get_software_module('r')

# MPRA enrichment parameters from config
MPRA_FILE = get_analysis_value(['enrichment', 'mpra', 'file'])
MPRA_FDR_THRESHOLD = get_analysis_value(['enrichment', 'mpra', 'fdr_threshold'])
MPRA_EXCEL_PVAL = get_analysis_value(['enrichment', 'mpra', 'excel_pvalue_threshold'])

# Rule targets
rule all:
    input:
        expand("results/08_mpra_enriched/{combination}/{combination}.done",
               combination=COMBINATIONS.keys())


rule enrich_mpra:
    """
    Enrich ANNOVAR-annotated results with MPRA functional data.
    - Filter MPRA for fdr < 0.05
    - Match by chr+pos (ignore alleles)
    - Aggregate multiple MPRA records per SNP:
      * Unique comma-separated diseases/celllines (sorted)
      * Max/min absolute log2FC
      * Count of MPRA records
    - Generate TSV with all SNPs (NA if no MPRA match)
    - Generate Excel with SNPs p < 5e-3
    """
    input:
        meta = "results/09_annovar_annotated/{combination}/{combination}_annotated.tsv",
        done = "results/09_annovar_annotated/{combination}/{combination}.done",
        mpra = MPRA_FILE
    output:
        tsv = "results/08_mpra_enriched/{combination}/{combination}_mpra_enriched.tsv",
        excel = "results/08_mpra_enriched/{combination}/{combination}_mpra_enriched.xlsx",
        done = "results/08_mpra_enriched/{combination}/{combination}.done"
    log:
        "results/log/08_mpra_enrich/{combination}.log"
    params:
        combination = "{combination}",
        r_module = R_MODULE,
        fdr_threshold = MPRA_FDR_THRESHOLD,
        excel_pval = MPRA_EXCEL_PVAL
    resources:
        mem_mb = 32000,
        time = "01:00:00",
        cores = 1
    shell:
        """
        set +u
        module load {params.r_module}
        set -u
        
        Rscript scripts/enrich_mpra.R \
            --input {input.meta} \
            --mpra {input.mpra} \
            --output {output.tsv} \
            --output_excel {output.excel} \
            --combination {params.combination} \
            --fdr_threshold {params.fdr_threshold} \
            --excel_pval {params.excel_pval} \
            --log {log} \
        2>&1 | tee -a {log}
        
        echo "MPRA enrichment completed for {params.combination}" > {output.done}
        echo "Input: {input.meta}" >> {output.done}
        echo "Output: {output.tsv}" >> {output.done}
        echo "Excel: {output.excel}" >> {output.done}
        """
