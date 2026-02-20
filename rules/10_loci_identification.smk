"""
Step 10: Genomic Loci Identification

This rule identifies independent genomic risk loci from meta-analysis results
using LD-based clumping with a reference panel. It:
1. Subsets the reference panel to SNPs present in the meta-analysis (cached for efficiency)
2. Performs LD-based clumping to identify independent significant variants
3. Merges correlated signals into genomic loci based on LD and distance

Input: Standardized Metal meta-analysis results from Step 6
Output: GenomicRiskLoci.txt with independent loci definitions
"""

import sys
import os
sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value

# Load configuration
results_dir = get_results_dir()

# Get loci identification parameters for each combination
loci_config = get_analysis_value(['loci_identification'], default={})

# Get reference panels configuration
reference_panels = get_analysis_value(['reference_panels'], default={}, required=False)
if not reference_panels:
    # Try from software.yml if not in analysis.yml
    import yaml
    with open('configs/software.yml', 'r') as f:
        software_config = yaml.safe_load(f)
        reference_panels = software_config.get('reference_panels', {})

# Default target rule
rule all:
    input:
        expand(
            os.path.join(results_dir, "10_loci", "{combination}", "{combination}_loci.done"),
            combination=loci_config.keys()
        )

rule loci_identification:
    input:
        meta = os.path.join(results_dir, "06_metal_standardized", "{combination}", "{combination}_standardized.tsv")
    output:
        done = os.path.join(results_dir, "10_loci", "{combination}", "{combination}_loci.done"),
        loci = os.path.join(results_dir, "10_loci", "{combination}", "GenomicRiskLoci.txt")
    params:
        outdir = lambda wildcards: os.path.join(results_dir, "10_loci", wildcards.combination),
        temp_dir = lambda wildcards: os.path.join(results_dir, "10_loci_temp"),
        ref_panel = lambda wildcards: loci_config.get(wildcards.combination, {}).get("reference_panel", ""),
        ref_bfile = lambda wildcards: reference_panels.get(
            loci_config.get(wildcards.combination, {}).get("reference_panel", ""), {}
        ).get("bfile", ""),
        target_build = get_analysis_value(['target_build'], default='hg38'),
        leadP = lambda wildcards: loci_config.get(wildcards.combination, {}).get("leadP", 5e-8),
        gwasP = lambda wildcards: loci_config.get(wildcards.combination, {}).get("gwasP", 1e-5),
        r2 = lambda wildcards: loci_config.get(wildcards.combination, {}).get("r2", 0.6),
        r2_2 = lambda wildcards: loci_config.get(wildcards.combination, {}).get("r2_2", 0.1),
        mergeDist = lambda wildcards: loci_config.get(wildcards.combination, {}).get("mergeDist", 250),
        maf = lambda wildcards: loci_config.get(wildcards.combination, {}).get("maf", 0.01),
        panel_match_p_threshold = lambda wildcards: loci_config.get(wildcards.combination, {}).get("panel_match_p_threshold", 5e-3),
        windowKb = lambda wildcards: loci_config.get(wildcards.combination, {}).get("windowKb", 500),
        threads = lambda wildcards: loci_config.get(wildcards.combination, {}).get("threads", 8),
        refSNPs = lambda wildcards: loci_config.get(wildcards.combination, {}).get("refSNPs", 1),
        plink_module = get_software_module("plink"),
        r_module = get_software_module("r")
    log:
        out = os.path.join(results_dir, "log", "identify_loci_{combination}.out"),
        err = os.path.join(results_dir, "log", "identify_loci_{combination}.err")
    resources:
        mem_mb = 64000,
        time = "02:00:00",
        cpus_per_task = 8
    shell:
        """
        # Load required modules
        module load {params.plink_module} 2>/dev/null || true
        module load {params.r_module}
        
        # Create output and temp directories
        mkdir -p {params.outdir}
        mkdir -p {params.temp_dir}
        
        # Run loci identification R script
        Rscript scripts/identify_loci.R \
            --meta {input.meta} \
            --ref-bfile {params.ref_bfile} \
            --ref-name {params.ref_panel} \
            --target-build {params.target_build} \
            --combination {wildcards.combination} \
            --temp-dir {params.temp_dir} \
            --outdir {params.outdir} \
            --leadP {params.leadP} \
            --gwasP {params.gwasP} \
            --r2 {params.r2} \
            --r2-2 {params.r2_2} \
            --mergeDist {params.mergeDist} \
            --maf {params.maf} \
            --panel-match-p-threshold {params.panel_match_p_threshold} \
            --windowKb {params.windowKb} \
            --threads {params.threads} \
            --refSNPs {params.refSNPs} \
            --plink plink \
            > {log.out} 2> {log.err}
        
        # Create done marker
        if [ -f {output.loci} ]; then
            touch {output.done}
            echo "Loci identification completed for {wildcards.combination}" >> {output.done}
        else
            echo "ERROR: GenomicRiskLoci.txt not created" >&2
            exit 1
        fi
        """
