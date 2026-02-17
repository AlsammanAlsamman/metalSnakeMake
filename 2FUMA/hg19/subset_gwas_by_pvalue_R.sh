#!/bin/bash
#SBATCH --job-name=subset_gwas
#SBATCH --output=log/subset_gwas_%j.out
#SBATCH --error=log/subset_gwas_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=128G
# subset_gwas_by_pvalue_R.sh
# R-based version with more detailed statistics and sorted output

ml R
set -e

PVAL_THRESHOLD="${1:-0.005}"
OUTPUT_SUFFIX=".p5e3"

echo "=========================================="
echo "Creating GWAS subsets with p-value < $PVAL_THRESHOLD"
echo "=========================================="

INPUT_FILES=(
#    "MVP.hg19.tsv"
#    "MVP_shared.hg19.tsv"
    "hisp_mvp_standardized.hg19.tsv"
#    "hisp_standardized.hg19.tsv"
)

for INPUT_FILE in "${INPUT_FILES[@]}"; do
    if [ ! -f "$INPUT_FILE" ]; then
        echo "Warning: $INPUT_FILE not found, skipping..."
        continue
    fi
    
    echo "----------------------------------------"
    echo "Processing: $INPUT_FILE"
    echo "----------------------------------------"
    
    BASENAME=$(basename "$INPUT_FILE" .tsv)
    OUTPUT_FILE="${BASENAME}${OUTPUT_SUFFIX}.tsv"
    STATS_FILE="${BASENAME}${OUTPUT_SUFFIX}.stats.txt"
    
    # Use R for more robust p-value handling and sorting
    Rscript -e "
        library(data.table)
        
        # Read the file
        cat('Reading $INPUT_FILE...\n')
        df <- fread('$INPUT_FILE', stringsAsFactors=FALSE, data.table=FALSE)
        
        cat('Total rows:', nrow(df), '\n')
        
        # Find chromosome and position columns
        chr_cols <- grep('^chr|^chrom', names(df), ignore.case=TRUE, value=TRUE)
        pos_cols <- grep('^pos|^position', names(df), ignore.case=TRUE, value=TRUE)
        
        if (length(chr_cols) == 0) {
            cat('Warning: Could not find chromosome column, using first column\n')
            chr_col <- names(df)[1]
        } else {
            chr_col <- chr_cols[1]
        }
        
        if (length(pos_cols) == 0) {
            cat('Warning: Could not find position column, using second column\n')
            pos_col <- names(df)[2]
        } else {
            pos_col <- pos_cols[1]
        }
        
        cat('Using chromosome column:', chr_col, '\n')
        cat('Using position column:', pos_col, '\n')
        
        # Find p-value column
        pval_cols <- grep('^p$|^pval$|^p_value$|^p-val', names(df), ignore.case=TRUE, value=TRUE)
        if (length(pval_cols) == 0) {
            pval_cols <- grep('p', names(df), ignore.case=TRUE, value=TRUE)
        }
        
        if (length(pval_cols) == 0) {
            stop('Could not find p-value column')
        }
        
        pval_col <- pval_cols[1]
        cat('Using p-value column:', pval_col, '\n')
        
        # Convert p-value to numeric (handles scientific notation)
        df[[pval_col]] <- as.numeric(df[[pval_col]])
        
        # Remove rows with NA p-values
        df <- df[!is.na(df[[pval_col]]), ]
        cat('Rows after removing NA p-values:', nrow(df), '\n')
        
        # Subset based on p-value threshold
        df_subset <- df[df[[pval_col]] < $PVAL_THRESHOLD, ]
        cat('Rows in subset:', nrow(df_subset), '\n')
        
        # Sort by chromosome and position
        if (nrow(df_subset) > 0) {
            cat('Sorting by chromosome and position...\n')
            
            # Function to extract numeric chromosome
            df_subset\$chr_num <- as.numeric(gsub('chr|chrom|chromosome|chr', '', 
                                                   df_subset[[chr_col]], ignore.case=TRUE))
            
            # Handle non-numeric chromosomes (X, Y, MT) by assigning large numbers
            df_subset\$chr_num[grepl('X|x', df_subset[[chr_col]])] <- 23
            df_subset\$chr_num[grepl('Y|y', df_subset[[chr_col]])] <- 24
            df_subset\$chr_num[grepl('M|MT|chrM', df_subset[[chr_col]])] <- 25
            
            # Sort by chromosome number then position
            df_subset <- df_subset[order(df_subset\$chr_num, df_subset[[pos_col]]), ]
            
            # Remove temporary column
            df_subset\$chr_num <- NULL
        }
        
        # Write subset
        cat('Writing sorted subset to $OUTPUT_FILE...\n')
        fwrite(df_subset, '$OUTPUT_FILE', sep='\t', na='NA', quote=FALSE)
        
        # Generate statistics
        sink('$STATS_FILE')
        cat('=== Subset Statistics for $INPUT_FILE ===\n')
        cat('P-value threshold: $PVAL_THRESHOLD\n')
        cat('Total variants in original:', nrow(df), '\n')
        cat('Variants in subset:', nrow(df_subset), '\n')
        cat('Percentage:', sprintf('%.2f%%', nrow(df_subset)/nrow(df)*100), '\n\n')
        
        if (nrow(df_subset) > 0) {
            cat('P-value distribution in subset:\n')
            print(summary(df_subset[[pval_col]]))
            
            cat('\nTop 10 most significant SNPs (by p-value):\n')
            top10 <- head(df_subset[order(df_subset[[pval_col]]), ], 10)
            print(top10[, c(chr_col, pos_col, pval_col, names(top10)[1:min(3, ncol(top10))])])
            
            cat('\nFirst 10 SNPs in genomic order:\n')
            first10 <- head(df_subset, 10)
            print(first10[, c(chr_col, pos_col, pval_col, names(first10)[1:min(3, ncol(first10))])])
        } else {
            cat('No variants meet the p-value threshold\n')
        }
        sink()
        
        cat('\nStatistics saved to: $STATS_FILE\n')
    "
    
    # Count results
    if [ -f "$OUTPUT_FILE" ]; then
        SUBSET_COUNT=$(($(wc -l < "$OUTPUT_FILE") - 1))
        echo "✓ Created subset with $SUBSET_COUNT SNPs"
        echo "  Output: $OUTPUT_FILE (sorted by chr:pos)"
        echo "  Stats: $STATS_FILE"
    fi
    
    echo ""
done

echo "=========================================="
echo "All subsets created successfully!"
echo "=========================================="