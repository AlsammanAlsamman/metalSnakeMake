#!/bin/bash
#
# Submit Step 05: METAL Meta-Analysis
# Runs meta-analysis for all study combinations defined in analysis.yml
#

# Create logs directory
mkdir -p logs results/log/05_metal

# Print header
echo "================================================================"
echo "Step 05: METAL Meta-Analysis"
echo "================================================================"
echo ""

# Parse command line arguments
DRY_RUN=""
EXTRA_FLAGS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run|-n)
            DRY_RUN="--dry-run"
            shift
            ;;
        *)
            EXTRA_FLAGS="$EXTRA_FLAGS $1"
            shift
            ;;
    esac
done

# Set up environment
module load slurm
module load python/3.7.0

# Configuration
SNAKEFILE="rules/05_metal.smk"
CORES=1
JOBS=10

# Print configuration
echo "Snakefile: $SNAKEFILE"
echo "Cores: $CORES"
echo "Max parallel jobs: $JOBS"
echo ""

if [[ -n "$DRY_RUN" ]]; then
    echo "DRY RUN MODE - No jobs will be submitted"
    echo ""
fi

# Submit the workflow to SLURM
snakemake \
    --cluster "sbatch \
        --mem={resources.mem_mb} \
        --time={resources.time} \
        --job-name={rule}_{wildcards.combination} \
        --cpus-per-task={resources.cores} \
        --output=logs/{rule}_{wildcards.combination}_%j.out \
        --error=logs/{rule}_{wildcards.combination}_%j.err" \
    --jobs $JOBS \
    --latency-wait 120 \
    --keep-going \
    --rerun-incomplete \
    --snakefile "$SNAKEFILE" \
    --cores $CORES \
    --verbose \
    --printshellcmds \
    --stats logs/snakemake_metal_stats.json \
    $DRY_RUN \
    $EXTRA_FLAGS

# Print completion message
echo ""
echo "================================================================"
echo "METAL meta-analysis jobs submitted"
echo "================================================================"
echo ""
echo "Monitor jobs: squeue -u $USER"
echo "Check logs: ls -lh logs/*metal*"
echo "Results: ls -lh results/05_metal/*/"
echo ""
