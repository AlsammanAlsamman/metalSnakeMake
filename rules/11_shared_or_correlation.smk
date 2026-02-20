#!/usr/bin/env python3
"""
Step 11: Shared SNP OR Correlation Across Datasets
Select shared significant SNPs (p < 5e-5), harmonize alleles/beta,
convert beta to OR, and evaluate cross-dataset OR concordance.
"""

import os
import sys
sys.path.append("utils")

from bioconfigme import (
    get_results_dir,
    get_software_module,
    get_metal_combinations
)


RESULTS_DIR = get_results_dir()
COMBINATIONS = get_metal_combinations()
R_MODULE = get_software_module('r')


rule all:
    input:
        expand(
            os.path.join(RESULTS_DIR, "11_shared_or", "{combination}", "{combination}.done"),
            combination=COMBINATIONS.keys()
        )


rule shared_or_correlation:
    input:
        tsv=lambda wildcards: expand(
            os.path.join(RESULTS_DIR, "04_metal_ready", "{dataset}.tsv"),
            dataset=COMBINATIONS[wildcards.combination]
        ),
        done=lambda wildcards: expand(
            os.path.join(RESULTS_DIR, "04_metal_ready", "{dataset}.done"),
            dataset=COMBINATIONS[wildcards.combination]
        )
    output:
        pdf=os.path.join(RESULTS_DIR, "11_shared_or", "{combination}", "{combination}_or_scatter.pdf"),
        summary=os.path.join(RESULTS_DIR, "11_shared_or", "{combination}", "{combination}_summary.tsv"),
        table=os.path.join(RESULTS_DIR, "11_shared_or", "{combination}", "{combination}_harmonized_or.tsv"),
        done=os.path.join(RESULTS_DIR, "11_shared_or", "{combination}", "{combination}.done")
    log:
        os.path.join(RESULTS_DIR, "log", "11_shared_or", "{combination}.log")
    params:
        combination="{combination}",
        datasets=lambda wildcards: ",".join(COMBINATIONS[wildcards.combination]),
        input_dir=os.path.join(RESULTS_DIR, "04_metal_ready"),
        r_module=R_MODULE
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2
    shell:
        """
        mkdir -p $(dirname {output.pdf}) $(dirname {log})

        set +u
        module load {params.r_module}
        set -u

        Rscript scripts/shared_or_correlation.R \
            --combination {params.combination} \
            --datasets {params.datasets} \
            --input_dir {params.input_dir} \
            --output_pdf {output.pdf} \
            --output_summary {output.summary} \
            --output_table {output.table} \
            --done {output.done} \
            --log {log} \
            2>&1 | tee -a {log}
        """
