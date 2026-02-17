#!/bin/bash
#SBATCH --job-name=liftover_gwas
#SBATCH --output=liftover_%j.out
#SBATCH --error=liftover_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=2

# liftover_gwas.sh
# This script converts GWAS summary statistics from hg38 to hg19 coordinates
# Usage: sbatch liftover_gwas.sh

ml R
ml ucsc

set -e  # Exit on error
set -o pipefail  # Exit if any command in a pipe fails

echo "=========================================="
echo "Starting liftOver conversion at: $(date)"
echo "=========================================="

# Configuration
CHAIN_FILE="hg38ToHg19.over.chain.gz"
INPUT_FILES=(
    # "MVP.tsv"
    # "MVP_shared.tsv" 
    # "hisp_mvp_metas_standardized.tsv"
    "hisp_mvp_standardized.tsv"
)
OUTPUT_SUFFIX=".hg19"
TEMP_DIR="tmp_liftover_$$"

# Create temporary directory
mkdir -p "$TEMP_DIR"
echo "Temporary directory: $TEMP_DIR"

# Check if chain file exists
if [ ! -f "$CHAIN_FILE" ]; then
    echo "Error: Chain file $CHAIN_FILE not found!"
    exit 1
fi

# Function to detect chromosome and position columns
detect_coord_columns() {
    local file=$1
    Rscript -e "
        library(data.table)
        df <- fread('$file', nrows=1)
        cols <- names(df)
        
        # Find chromosome column
        chr_col <- grep('^chr|chrom|chromosome', cols, ignore.case=TRUE, value=TRUE)
        if (length(chr_col) == 0) chr_col <- 'chrom'
        
        # Find position column
        pos_col <- grep('^pos|position|bp|coord', cols, ignore.case=TRUE, value=TRUE)
        if (length(pos_col) == 0) pos_col <- 'pos'
        
        cat(paste(chr_col[1], pos_col[1], sep='\n'))
    "
}

# Process each input file
for INPUT_FILE in "${INPUT_FILES[@]}"; do
    if [ ! -f "$INPUT_FILE" ]; then
        echo "Warning: $INPUT_FILE not found, skipping..."
        continue
    fi
    
    echo "=========================================="
    echo "Processing: $INPUT_FILE at $(date)"
    echo "=========================================="
    
    BASENAME=$(basename "$INPUT_FILE" .tsv)
    OUTPUT_FILE="${BASENAME}${OUTPUT_SUFFIX}.tsv"
    BED_FILE="$TEMP_DIR/${BASENAME}.bed"
    LIFTED_BED="$TEMP_DIR/${BASENAME}.lifted.bed"
    UNMAPPED_FILE="$TEMP_DIR/${BASENAME}.unmapped.txt"
    
    # Detect chromosome and position columns
    COORDS=$(detect_coord_columns "$INPUT_FILE")
    CHR_COL=$(echo "$COORDS" | head -1)
    POS_COL=$(echo "$COORDS" | tail -1)
    
    echo "Detected columns - Chromosome: $CHR_COL, Position: $POS_COL"
    
    # Step 1: Create BED file from input (adding 0-based start)
    echo "Creating BED file..."
    Rscript -e "
        library(data.table)
        
        # Read input file
        cat('Reading $INPUT_FILE...\n')
        df <- fread('$INPUT_FILE', stringsAsFactors=FALSE, data.table=FALSE)
        
        cat('Total rows:', nrow(df), '\n')
        
        # Convert position to integer (handle scientific notation)
        cat('Converting positions to integers...\n')
        df[['$POS_COL']] <- as.integer(round(df[['$POS_COL']]))
        
        # Remove any rows with NA positions
        df <- df[!is.na(df[['$POS_COL']]), ]
        cat('Rows after removing NA positions:', nrow(df), '\n')
        
        # Create BED format (0-based start, 1-based end)
        bed <- data.frame(
            chr = df[['$CHR_COL']],
            start = as.integer(df[['$POS_COL']] - 1),  # 0-based start
            end = as.integer(df[['$POS_COL']]),         # 1-based end
            id = 1:nrow(df),                             # Row ID for mapping back
            stringsAsFactors = FALSE
        )
        
        # Ensure chromosome has 'chr' prefix for liftOver
        bed\$chr <- ifelse(grepl('^chr', bed\$chr), bed\$chr, paste0('chr', bed\$chr))
        
        # Remove any rows with invalid positions (negative or zero)
        bed <- bed[bed\$start >= 0 & bed\$end > bed\$start, ]
        cat('Rows after removing invalid positions:', nrow(bed), '\n')
        
        # Write BED file with integer formatting
        fwrite(bed, '$BED_FILE', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE, 
               scipen=999)  # Prevent scientific notation
        cat('BED file created with', nrow(bed), 'rows\n')
    "
    
    # Check if BED file was created and has content
    if [ ! -s "$BED_FILE" ]; then
        echo "Error: BED file is empty for $INPUT_FILE"
        continue
    fi
    
    # Step 2: Run liftOver
    echo "Running liftOver..."
    if [ -f "$LIFTED_BED" ]; then
        rm "$LIFTED_BED" "$UNMAPPED_FILE" 2>/dev/null || true
    fi
    
    liftOver "$BED_FILE" "$CHAIN_FILE" "$LIFTED_BED" "$UNMAPPED_FILE"
    
    # Check liftOver results
    if [ ! -f "$LIFTED_BED" ]; then
        echo "Error: liftOver failed for $INPUT_FILE"
        continue
    fi
    
    LIFTED_COUNT=$(wc -l < "$LIFTED_BED")
    UNMAPPED_COUNT=$(wc -l < "$UNMAPPED_FILE" 2>/dev/null || echo 0)
    
    echo "liftOver complete:"
    echo "  Successfully lifted: $LIFTED_COUNT variants"
    echo "  Failed to lift: $UNMAPPED_COUNT variants"
    
    # Step 3: Merge lifted coordinates back to original data
    echo "Creating output file with hg19 coordinates..."
    Rscript -e "
        library(data.table)
        
        # Read original data (keep all original values including scientific notation)
        cat('Reading original file...\n')
        original <- fread('$INPUT_FILE', stringsAsFactors=FALSE, data.table=FALSE)
        
        # Add ID column to original
        original\$id <- 1:nrow(original)
        
        # Read lifted BED file
        if (file.exists('$LIFTED_BED') && file.size('$LIFTED_BED') > 0) {
            cat('Reading lifted coordinates...\n')
            lifted <- fread('$LIFTED_BED', stringsAsFactors=FALSE, data.table=FALSE,
                           col.names=c('chr_hg19', 'start_hg19', 'end_hg19', 'id'))
            
            # Calculate new position (end is 1-based position)
            lifted\$pos_hg19 <- lifted\$end_hg19
            
            # Remove 'chr' prefix if original didn't have it
            if (!grepl('^chr', as.character(original[1,'$CHR_COL']))) {
                lifted\$chr_hg19 <- gsub('^chr', '', lifted\$chr_hg19)
            }
            
            # Merge lifted coordinates back
            cat('Merging coordinates...\n')
            result <- merge(original, lifted[, c('id', 'chr_hg19', 'pos_hg19')], 
                           by='id', all.x=TRUE, sort=FALSE)
        } else {
            result <- original
            result\$chr_hg19 <- NA
            result\$pos_hg19 <- NA
        }
        
        # Update chromosome and position columns (only for successfully lifted variants)
        result[['$CHR_COL']] <- ifelse(!is.na(result\$chr_hg19), 
                                        result\$chr_hg19, 
                                        result[['$CHR_COL']])
        result[['$POS_COL']] <- ifelse(!is.na(result\$pos_hg19), 
                                        result\$pos_hg19, 
                                        result[['$POS_COL']])
        
        # Remove temporary columns
        result\$id <- NULL
        result\$chr_hg19 <- NULL
        result\$pos_hg19 <- NULL
        
        # Write output file
        cat('Writing output to $OUTPUT_FILE...\n')
        fwrite(result, '$OUTPUT_FILE', sep='\t', na='NA', quote=FALSE)
        
        # Print summary
        cat('\n=== Summary for $INPUT_FILE ===\n')
        cat('Total variants in original:', nrow(original), '\n')
        cat('Successfully converted:', sum(!is.na(result[['$POS_COL']]) & 
                                           result[['$POS_COL']] != original[['$POS_COL']]), '\n')
        
        # Show sample of converted variants
        converted <- which(!is.na(result[['$POS_COL']]) & 
                          result[['$POS_COL']] != original[['$POS_COL']])
        if (length(converted) > 0) {
            cat('\nSample of converted variants (first 5):\n')
            for (i in head(converted, 5)) {
                cat(sprintf('  %s: %s %s -> %s %s\n', 
                           original[i,'$CHR_COL'],
                           format(original[i,'$POS_COL'], scientific=FALSE),
                           '->',
                           result[i,'$CHR_COL'],
                           format(result[i,'$POS_COL'], scientific=FALSE)))
            }
        }
    "
    
    echo "✓ Output file created: $OUTPUT_FILE"
    echo ""
    
done

# Clean up
echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo "=========================================="
echo "All conversions completed at: $(date)"
echo "=========================================="

# Final summary
echo -e "\n=== Final Output Files ==="
for INPUT_FILE in "${INPUT_FILES[@]}"; do
    if [ -f "$INPUT_FILE" ]; then
        BASENAME=$(basename "$INPUT_FILE" .tsv)
        OUTPUT_FILE="${BASENAME}${OUTPUT_SUFFIX}.tsv"
        if [ -f "$OUTPUT_FILE" ]; then
            ORIG_COUNT=$(wc -l < "$INPUT_FILE")
            NEW_COUNT=$(wc -l < "$OUTPUT_FILE")
            echo "$OUTPUT_FILE: $NEW_COUNT lines (original: $ORIG_COUNT)"
        fi
    fi
done

echo -e "\nNote: Variants that failed liftOver kept their original coordinates"