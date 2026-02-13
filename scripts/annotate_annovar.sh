#!/bin/bash
# ANNOVAR Annotation Wrapper Script
# Converts meta-analysis to ANNOVAR format, runs annotation, and merges results

set -e
set -o pipefail

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input) INPUT="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --output_excel) OUTPUT_EXCEL="$2"; shift 2 ;;
        --combination) COMBINATION="$2"; shift 2 ;;
        --buildver) BUILDVER="$2"; shift 2 ;;
        --humandb) HUMANDB="$2"; shift 2 ;;
        --protocols) PROTOCOLS="$2"; shift 2 ;;
        --operations) OPERATIONS="$2"; shift 2 ;;
        --excel_pval) EXCEL_PVAL="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --log) LOG="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Log function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"
}

log "=== ANNOVAR Annotation Pipeline ==="
log "Combination: $COMBINATION"
log "Input: $INPUT"
log "Build version: $BUILDVER"
log "Protocols: $PROTOCOLS"
log "Operations: $OPERATIONS"

# Create temporary directory
TMPDIR="${OUTDIR}/tmp"
mkdir -p "$TMPDIR"

# Step 1: Convert to ANNOVAR input format (chr, start, end, ref, alt)
ANNOVAR_INPUT="${TMPDIR}/${COMBINATION}.avinput"
log "Converting to ANNOVAR input format..."

# Extract chr, pos, pos, ref(nea), alt(ea) and create unique key
# ANNOVAR format: chr start end ref alt otherinfo
awk 'BEGIN {FS=OFS="\t"} 
NR==1 {
    # Find column indices
    for(i=1; i<=NF; i++) {
        col[$i] = i
    }
    next
}
{
    # Extract columns: chrom, pos, pos, nea(ref), ea(alt), markername
    chrom = $col["chrom"]
    pos = $col["pos"]
    nea = $col["nea"]  # ref allele
    ea = $col["ea"]    # alt allele
    markername = $col["markername"]
    
    # Print ANNOVAR format: chr start end ref alt key
    print chrom, pos, pos, nea, ea, markername
}' "$INPUT" > "$ANNOVAR_INPUT"

VARIANT_COUNT=$(wc -l < "$ANNOVAR_INPUT")
log "Variants to annotate: $VARIANT_COUNT"

# Step 2: Run ANNOVAR table_annovar.pl
ANNOVAR_PREFIX="${TMPDIR}/${COMBINATION}"
log "Running ANNOVAR table_annovar.pl..."

table_annovar.pl \
    "$ANNOVAR_INPUT" \
    "$HUMANDB" \
    -buildver "$BUILDVER" \
    -out "$ANNOVAR_PREFIX" \
    -remove \
    -protocol "$PROTOCOLS" \
    -operation "$OPERATIONS" \
    -nastring . \
    -polish \
    2>&1 | tee -a "$LOG"

# Check if annotation completed
ANNOVAR_OUTPUT="${ANNOVAR_PREFIX}.${BUILDVER}_multianno.txt"
if [[ ! -f "$ANNOVAR_OUTPUT" ]]; then
    log "ERROR: ANNOVAR output not found: $ANNOVAR_OUTPUT"
    exit 1
fi

log "ANNOVAR annotation completed: $ANNOVAR_OUTPUT"

# Step 3: Merge ANNOVAR results with original enriched file using R
log "Merging ANNOVAR annotations with enriched file..."

Rscript scripts/merge_annovar.R \
    --input "$INPUT" \
    --annovar "$ANNOVAR_OUTPUT" \
    --output "$OUTPUT" \
    --output_excel "$OUTPUT_EXCEL" \
    --excel_pval "$EXCEL_PVAL" \
    --combination "$COMBINATION" \
    --log "$LOG"

# Clean up temporary files
log "Cleaning up temporary files..."
rm -rf "$TMPDIR"

log "=== ANNOVAR Annotation Complete ==="
