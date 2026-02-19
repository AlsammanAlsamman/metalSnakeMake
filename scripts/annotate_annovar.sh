#!/bin/bash
# ANNOVAR Annotation Wrapper Script
# Converts meta-analysis to ANNOVAR format, runs annotation, and merges results

# Debug output to stderr
echo "=== ANNOVAR script starting ===" >&2
echo "Arguments: $@" >&2

set -e
set -o pipefail

# Parse arguments FIRST
while [[ $# -gt 0 ]]; do
    case $1 in
        --input) INPUT="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --output_excel) OUTPUT_EXCEL="$2"; shift 2 ;;
        --combination) COMBINATION="$2"; shift 2 ;;
        --studies) STUDIES="$2"; shift 2 ;;
        --buildver) BUILDVER="$2"; shift 2 ;;
        --humandb) HUMANDB="$2"; shift 2 ;;
        --protocols) PROTOCOLS="$2"; shift 2 ;;
        --operations) OPERATIONS="$2"; shift 2 ;;
        --excel_pval) EXCEL_PVAL="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --log) LOG="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

echo "Arguments parsed successfully" >&2
echo "LOG=$LOG" >&2

# Ensure log directory exists
mkdir -p "$(dirname "$LOG")"
echo "Log directory created" >&2

# Log function (defined AFTER LOG variable is set)
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"
}

log "=== ANNOVAR Annotation Pipeline ==="
log "Combination: $COMBINATION"
log "Studies: $STUDIES"
log "Input: $INPUT"
log "Build version: $BUILDVER"
log "Protocols: $PROTOCOLS"
log "Operations: $OPERATIONS"

# Check if ANNOVAR is available
if ! command -v table_annovar.pl &> /dev/null; then
    log "ERROR: table_annovar.pl not found in PATH"
    log "Please ensure ANNOVAR module is loaded correctly"
    exit 1
fi

log "ANNOVAR found: $(which table_annovar.pl)"

# Check if humandb directory exists
if [[ ! -d "$HUMANDB" ]]; then
    log "ERROR: ANNOVAR humandb directory not found: $HUMANDB"
    log "Please download ANNOVAR databases first"
    exit 1
fi

log "ANNOVAR humandb: $HUMANDB"

# Check if input file exists
if [[ ! -f "$INPUT" ]]; then
    log "ERROR: Input file not found: $INPUT"
    exit 1
fi

# Create temporary directory
TMPDIR="${OUTDIR}/tmp"
mkdir -p "$TMPDIR"

# Step 1: Convert to ANNOVAR input format and create markername mapping
ANNOVAR_INPUT="${TMPDIR}/${COMBINATION}.avinput"
MARKERNAME_MAP="${OUTDIR}/${COMBINATION}_markername_map.txt"
log "Converting to ANNOVAR input format..."

# Create ANNOVAR input file (chr, start, end, ref, alt)
awk 'BEGIN {FS=OFS="\t"} 
NR==1 {
    for(i=1; i<=NF; i++) col[$i] = i
    next
}
{
    print $col["chrom"], $col["pos"], $col["pos"], $col["nea"], $col["ea"]
}' "$INPUT" > "$ANNOVAR_INPUT"

# Create markername mapping file (chr:pos:ref:alt -> markername)
awk 'BEGIN {FS=OFS="\t"} 
NR==1 {
    for(i=1; i<=NF; i++) col[$i] = i
    next
}
{
    merge_key = $col["chrom"] ":" $col["pos"] ":" $col["nea"] ":" $col["ea"]
    print merge_key, $col["markername"]
}' "$INPUT" > "$MARKERNAME_MAP"

VARIANT_COUNT=$(wc -l < "$ANNOVAR_INPUT")
log "Variants to annotate: $VARIANT_COUNT"
log "Markername mapping created: $MARKERNAME_MAP"

# Step 2: Run ANNOVAR table_annovar.pl
ANNOVAR_PREFIX="${TMPDIR}/${COMBINATION}"
log "Running ANNOVAR table_annovar.pl..."

table_annovar.pl \
    "$ANNOVAR_INPUT" \
    "$HUMANDB" \
    -buildver "$BUILDVER" \
    -out "$ANNOVAR_PREFIX" \
    -protocol "$PROTOCOLS" \
    -operation "$OPERATIONS" \
    -nastring . \
    -polish \
    --thread 1 \
    --maxgenethread 1 \
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
    --markername_map "$MARKERNAME_MAP" \
    --output "$OUTPUT" \
    --output_excel "$OUTPUT_EXCEL" \
    --excel_pval "$EXCEL_PVAL" \
    --combination "$COMBINATION" \
    --studies "$STUDIES" \
    --log "$LOG"

# Step 4: Create Excel file from annotated TSV (separate process to avoid sink issues)
log "Creating Excel output..."

Rscript scripts/create_excel_v2.R \
    --input "$OUTPUT" \
    --output "$OUTPUT_EXCEL" \
    --pval_threshold "$EXCEL_PVAL" \
    --combination "$COMBINATION" \
    2>&1 | tee -a "$LOG"

# Clean up temporary files
log "Cleaning up temporary files..."
rm -rf "$TMPDIR"

log "=== ANNOVAR Annotation Complete ==="
