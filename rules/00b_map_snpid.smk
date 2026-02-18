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
BCFTOOLS_MODULE = get_software_module('bcftools', required=False, default='')
CHR_ORDER = [str(i) for i in range(1, 24)] + ["X", "Y"]

rule all:
    input:
        expand(f"{RESULTS_DIR}/snpdmapped/{{dataset}}.done", dataset=DATASETS)


rule map_snpid_from_snpdb_chr:
    """
    Map rsIDs for one chromosome; each chromosome is its own job with persistent outputs.
    """
    input:
        tsv = "results/00_standardized/{dataset}.tsv",
        done = "results/00_standardized/{dataset}.done"
    output:
        tsv = f"{RESULTS_DIR}/snpdmapped/{{dataset}}/mapping/chr_{{chrom}}.mapped.tsv",
        done = f"{RESULTS_DIR}/snpdmapped/{{dataset}}/mapping/chr_{{chrom}}.done"
    log:
        f"{RESULTS_DIR}/log/00b_map_snpid/{{dataset}}_chr{{chrom}}.log"
    params:
        dataset = "{dataset}",
        chrom = "{chrom}",
        build = lambda wildcards: get_dataset_config(wildcards.dataset).get('build', 'hg38'),
        snpdb_root = "resources/SNPdb",
        mapping_dir = f"{RESULTS_DIR}/snpdmapped/{{dataset}}/mapping",
        r_module = R_MODULE,
        bcftools_module = BCFTOOLS_MODULE
    threads: 1
    resources:
        mem_mb = 32000,
        time = "1:00:00"
    shell:
        """
        set +u
        module load {params.r_module}
        if [[ -n "{params.bcftools_module}" ]]; then
            module load {params.bcftools_module} || module --ignore_cache load {params.bcftools_module} || true
        fi
        set -u

        if ! command -v bcftools >/dev/null 2>&1; then
            echo "Error: bcftools is not available in PATH on this node. Update configs/software.yml for bcftools module/path." >&2
            exit 1
        fi

        Rscript scripts/map_snpid_from_snpdb_chr.R \
            --dataset {params.dataset} \
            --chrom {params.chrom} \
            --build {params.build} \
            --input {input.tsv} \
            --output {output.tsv} \
            --done {output.done} \
            --log {log} \
            --snpdb-root {params.snpdb_root} \
            --mapping-dir {params.mapping_dir} \
            2>&1 | tee {log}
        """


rule merge_snpid_from_snpdb:
    """
    Merge per-chromosome mapped outputs in fixed chromosome order.
    """
    input:
        per_chr_done = lambda wildcards: expand(
            f"{RESULTS_DIR}/snpdmapped/{{dataset}}/mapping/chr_{{chrom}}.done",
            dataset=[wildcards.dataset],
            chrom=CHR_ORDER
        ),
        per_chr_tsv = lambda wildcards: expand(
            f"{RESULTS_DIR}/snpdmapped/{{dataset}}/mapping/chr_{{chrom}}.mapped.tsv",
            dataset=[wildcards.dataset],
            chrom=CHR_ORDER
        )
    output:
        tsv = f"{RESULTS_DIR}/snpdmapped/{{dataset}}.tsv",
        done = f"{RESULTS_DIR}/snpdmapped/{{dataset}}.done",
        summary = f"{RESULTS_DIR}/snpdmapped/summary/{{dataset}}_mapping_summary.tsv"
    log:
        f"{RESULTS_DIR}/log/00b_map_snpid/{{dataset}}_merge.log"
    params:
        dataset = "{dataset}",
        chrom_order = ",".join(CHR_ORDER),
        inputs_csv = lambda wildcards, input: ",".join(list(input.per_chr_tsv)),
        r_module = R_MODULE
    threads: 1
    resources:
        mem_mb = 32000,
        time = "0:30:00"
    shell:
        """
        set +u
        module load {params.r_module}
        set -u

        Rscript scripts/merge_snpid_chr_outputs.R \
            --dataset {params.dataset} \
            --chrom-order "{params.chrom_order}" \
            --inputs-csv "{params.inputs_csv}" \
            --output {output.tsv} \
            --done {output.done} \
            --summary {output.summary} \
            --log {log} \
            2>&1 | tee {log}
        """