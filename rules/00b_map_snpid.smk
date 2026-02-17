#!/usr/bin/env python3
"""
Step 00b: Map SNP IDs (rsID) from SNPdb
Maps variants from standardized datasets to rsIDs using chr+pos and allele harmonization.
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
BCFTOOLS_MODULE = get_software_module('bcftools')

rule all:
    input:
        expand(f"{RESULTS_DIR}/snpdmapped/{{dataset}}.done", dataset=DATASETS)


rule map_snpid_from_snpdb:
    """
    Map rsIDs using SNPdb with chromosome-split, position-restricted lookup and allele harmonization.
    """
    input:
        tsv = "results/00_standardized/{dataset}.tsv",
        done = "results/00_standardized/{dataset}.done"
    output:
        tsv = f"{RESULTS_DIR}/snpdmapped/{{dataset}}.tsv",
        done = f"{RESULTS_DIR}/snpdmapped/{{dataset}}.done",
        summary = f"{RESULTS_DIR}/snpdmapped/summary/{{dataset}}_mapping_summary.tsv"
    log:
        f"{RESULTS_DIR}/log/00b_map_snpid/{{dataset}}.log"
    params:
        dataset = "{dataset}",
        build = lambda wildcards: get_dataset_config(wildcards.dataset).get('build', 'hg38'),
        snpdb_root = "resources/SNPdb",
        temp_dir = f"{RESULTS_DIR}/snpdmapped/tmp/{{dataset}}",
        r_module = R_MODULE,
        bcftools_module = BCFTOOLS_MODULE
    resources:
        mem_mb = 32000,
        time = "00:30:00",
        cores = 2
    shell:
        """
        set +u
        module load {params.r_module}
        module load {params.bcftools_module}
        set -u

        Rscript scripts/map_snpid_from_snpdb.R \
            --dataset {params.dataset} \
            --build {params.build} \
            --input {input.tsv} \
            --output {output.tsv} \
            --done {output.done} \
            --summary {output.summary} \
            --log {log} \
            --snpdb-root {params.snpdb_root} \
            --temp-dir {params.temp_dir} \
            --threads {resources.cores} \
            2>&1 | tee {log}
        """