#!/usr/bin/env python3
"""
Step 00c: Recover rsIDs for previously unmapped variants using Ensembl REST API.
Reads Step 00b unmapped variants, queries API per SNP, merges recovered rsIDs back into mapped table.
"""

import sys
sys.path.append("utils")

from bioconfigme import (
    get_results_dir,
    get_dataset_list,
    get_dataset_config,
    get_software_module
)

RESULTS_DIR = get_results_dir()
DATASETS = get_dataset_list()
R_MODULE = get_software_module('r')

rule all:
    input:
        expand(f"{RESULTS_DIR}/snpdmapped/{{dataset}}_mapped.done", dataset=DATASETS)


rule api_map_unmapped_snps:
    """
    Read Step 00b output, detect fallback snpids (chr:pos:*:*), query Ensembl API, and update rsIDs in place.
    """
    input:
        mapped_tsv = f"{RESULTS_DIR}/snpdmapped/{{dataset}}.tsv",
        mapped_done = f"{RESULTS_DIR}/snpdmapped/{{dataset}}.done"
    output:
        tsv = f"{RESULTS_DIR}/snpdmapped/{{dataset}}_mapped.tsv",
        done = f"{RESULTS_DIR}/snpdmapped/{{dataset}}_mapped.done",
        summary = f"{RESULTS_DIR}/snpdmapped/summary/{{dataset}}_api_mapping_summary.tsv"
    log:
        f"{RESULTS_DIR}/log/00c_api_map_unmapped/{{dataset}}.log"
    params:
        dataset = "{dataset}",
        build = lambda wildcards: get_dataset_config(wildcards.dataset).get('build', 'hg38'),
        r_module = R_MODULE,
        rate_per_second = 10,
        max_queries = 0
    threads: 1
    resources:
        mem_mb = 16000,
        time = "2:00:00"
    shell:
        """
        set +u
        module load {params.r_module}
        set -u

        if ! command -v curl >/dev/null 2>&1; then
            echo "Error: curl is not available in PATH on this node." >&2
            exit 1
        fi

        Rscript scripts/remap_fallback_snpid_with_ensembl_api.R \
            --dataset {params.dataset} \
            --build {params.build} \
            --input {input.mapped_tsv} \
            --output {output.tsv} \
            --done {output.done} \
            --summary {output.summary} \
            --log {log} \
            --rate-per-second {params.rate_per_second} \
            --max-queries {params.max_queries} \
            2>&1 | tee {log}
        """
