#!/bin/bash

# Standarizing
./submit.sh --snakefile rules/00_standardize.smk --cores 2 --jobs 7
./submit.sh --snakefile rules/01_calculate_missing.smk --cores 2 --jobs 7
./submit.sh --snakefile rules/02_filter_snps.smk --cores 1 --jobs 7
./submit.sh --snakefile rules/03_liftover.smk --cores 1 --jobs 7

./submit.sh --snakefile rules/04_prepare_metal.smk --cores 1 --jobs 7
./submit.sh --snakefile rules/11_shared_or_correlation.smk --cores 1 --jobs 10
./submit.sh --snakefile rules/05_metal.smk --cores 1 --jobs 10
./submit.sh --snakefile rules/06_standardize_metal.smk --cores 1 --jobs 10
./submit.sh --snakefile rules/09_annovar_annotate.smk --cores 1 --jobs 10



./submit.sh --snakefile rules/00b_map_snpid.smk --cores 1 --jobs 22 /s/nath-lab/alsamman/____MyCodes____/metalSnakeMake/results/snpdmapped/TWN.done


###############################
# METAL Meta-Analysis Pipeline Execution Steps
# Each step can be run independently - previous steps must complete successfully (check .done files)

# Step -1: Split SNPdb by chromosome - Pre-split build VCFs for faster parallel rsID annotation
# ./utility_scripts/split_snpdb_by_chr.sh --input-dir resources/SNPdb --output-root resources/SNPdb/by_chr --workers 16


# Step -1 (alternative with overwrite of existing split files):
# ./utility_scripts/split_snpdb_by_chr.sh --input-dir resources/SNPdb --output-root resources/SNPdb/by_chr --workers 16 --force

# Step 0: Standardization - Extract required columns and standardize column names
# Process all 7 datasets in parallel (7 SLURM jobs)
./scripts/submit_standardize.sh





# Step 0 (alternative - all datasets via command line):
# ./submit.sh --snakefile rules/00_standardize.smk --cores 2 --jobs 7


# Step 0b: Map SNP IDs (rsID) from SNPdb - chromosome split + allele harmonization
# Process all datasets with per-chromosome jobs.
# Outputs are persistent under: results/snpdmapped/{dataset}/mapping/chr_{chrom}.mapped.tsv + .done
# Re-run is resumable because completed chromosome .done files are skipped.
bash ./scripts/submit_snpid_map.sh


# Step 0b (alternative - all datasets via command line):
# ./submit.sh --snakefile rules/00b_map_snpid.smk --cores 1 --jobs 22

# Step 0b (Asian only):
# ./submit.sh --snakefile rules/00b_map_snpid.smk --cores 1 --jobs 22 /s/nath-lab/alsamman/____MyCodes____/metalSnakeMake/results/snpdmapped/Asian.done


# Step 0c: API remap for SNPs still unmapped after Step 0b
# Reads results/snpdmapped/{dataset}.tsv, finds fallback snpids (chr:pos:ea:nea), and remaps via API
# Output: results/snpdmapped/{dataset}_mapped.tsv
./scripts/submit_snpid_api_map.sh


# Step 0c (alternative - all datasets via command line):
# ./submit.sh --snakefile rules/00c_api_map_unmapped.smk --cores 1 --jobs 7

# Step 0c (Asian only):
# ./submit.sh --snakefile rules/00c_api_map_unmapped.smk --cores 1 /s/nath-lab/alsamman/____MyCodes____/metalSnakeMake/results/snpdmapped/Asian_mapped.done


# Step 1: Calculate Missing Columns - Calculate SE from beta, eaf, n_cases, n_controls
# Process all datasets in parallel (7 SLURM jobs)
./scripts/submit_calculate_missing.sh

# Step 1 (alternative - all datasets via command line):
# ./submit.sh --snakefile rules/01_calculate_missing.smk --cores 2 --jobs 7


# Step 2: Filter Target SNPs - Filter to target SNP list or copy all variants
# Process all datasets in parallel (7 SLURM jobs)
./scripts/submit_filter_snps.sh

# Step 2 (alternative - all datasets via command line):
# ./submit.sh --snakefile rules/02_filter_snps.smk --cores 1 --jobs 7


# Step 3: LiftOver Genome Coordinates - Convert builds to target build (hg38)
# Process all datasets in parallel (7 SLURM jobs)
./scripts/submit_liftover.sh

# Step 3 (alternative - all datasets via command line):
# ./submit.sh --snakefile rules/03_liftover.smk --cores 1 --jobs 7


# Step 4: Prepare for METAL - Add MarkerName column for meta-analysis
# Process all datasets in parallel (7 SLURM jobs)
./scripts/submit_prepare_metal.sh

# Step 4 (alternative - all datasets via command line):
# ./submit.sh --snakefile rules/04_prepare_metal.smk --cores 1 --jobs 7


# Step 5: METAL Meta-Analysis - Run inverse-variance weighted meta-analysis
# Process all study combinations in parallel (1 job per combination)
./submit_metal.sh

# Step 5 (alternative - all combinations via command line):
# ./submit.sh --snakefile rules/05_metal.smk --cores 1 --jobs 10

# Step 5 (dry run to check configuration):
# ./submit_metal.sh --dry-run


# Step 6: Standardize METAL Output - Standardize results and generate Manhattan plots
# Process all combinations in parallel (1 job per combination)
./scripts/submit_standardize_metal.sh

# Step 6 (alternative - all combinations via command line):
# ./submit.sh --snakefile rules/06_standardize_metal.smk --cores 1 --jobs 10

# Step 6 (dry run to check configuration):
# ./scripts/submit_standardize_metal.sh --dry-run


# Step 7: Enrich Meta-Analysis with Study-Specific EAF - Add study EAFs and calculate meta EAF
# Process all combinations in parallel (1 job per combination)
./scripts/submit_enrich_eaf.sh

# Step 7 (alternative - all combinations via command line):
# ./submit.sh --snakefile rules/07_enrich_eaf.smk --cores 1 --jobs 10

# Step 7 (dry run to check configuration):
# ./scripts/submit_enrich_eaf.sh --dry-run


# Step 8: Enrich Meta-Analysis with MPRA Functional Data - Match SNPs with MPRA by chr+pos
# Process all combinations in parallel (1 job per combination)
./submit.sh --snakefile rules/08_mpra_enrich.smk --cores 1 --jobs 10

# Step 8 (dry run to check configuration):
# ./submit.sh --snakefile rules/08_mpra_enrich.smk --cores 1 --jobs 10 --dry-run


# Step 9: Annotate Meta-Analysis with ANNOVAR - Add functional gene-based annotations
# Process all combinations in parallel (1 job per combination)
./submit.sh --snakefile rules/09_annovar_annotate.smk --cores 1 --jobs 10

# Step 9 (dry run to check configuration):
# ./submit.sh --snakefile rules/09_annovar_annotate.smk --cores 1 --jobs 10 --dry-run


# Step 10: Identify Genomic Loci - LD-based clumping to identify independent risk loci
# Process all combinations in parallel (1 job per combination, 8 cores each)
./submit.sh --snakefile rules/10_loci_identification.smk --cores 8 --jobs 10

# Step 10 (single combination test):
# ./submit.sh --snakefile rules/10_loci_identification.smk results/10_loci/hisp/hisp_loci.done

# Step 10 (dry run to check configuration):
# ./submit.sh --snakefile rules/10_loci_identification.smk --cores 8 --jobs 10 --dry-run


# Step 11: Shared SNP OR Correlation - harmonize alleles/beta and assess OR concordance
./submit.sh --snakefile rules/11_shared_or_correlation.smk --cores 1 --jobs 10

# Step 11 (single combination test):
# ./submit.sh --snakefile rules/11_shared_or_correlation.smk /s/nath-lab/alsamman/____MyCodes____/metalSnakeMake/results/11_shared_or/Asian_EUR_Hisp/Asian_EUR_Hisp.done