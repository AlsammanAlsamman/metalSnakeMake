#!/usr/bin/env bash
#SBATCH --job-name=filter_twn_p
#SBATCH --output=slurm-filter_twn_p_%j.out
#SBATCH --error=slurm-filter_twn_p_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

set -euo pipefail

timestamp() {
  date "+%Y-%m-%d %H:%M:%S"
}

log() {
  echo "[$(timestamp)] $*"
}

on_error() {
  local exit_code=$?
  local line_no=$1
  local cmd=${2:-"<unknown>"}
  echo "[$(timestamp)] ERROR: command failed at line $line_no: $cmd (exit=$exit_code)" >&2
  exit "$exit_code"
}

trap 'on_error ${LINENO} "${BASH_COMMAND}"' ERR

# Usage:
#   sbatch utility_scripts/filter_twn_pvalues.sbatch.sh
#   sbatch utility_scripts/filter_twn_pvalues.sbatch.sh input/TWN.txt

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

INPUT_FILE="${1:-input/TWN.txt}"
if [[ "$INPUT_FILE" != /* ]]; then
  INPUT_FILE="$PROJECT_ROOT/$INPUT_FILE"
fi

OUT_DIR="$PROJECT_ROOT/utility_out"
RUN_LOG="$OUT_DIR/filter_twn_pvalues_runtime.log"

OUT_P01_TSV="$OUT_DIR/TWN_p_lt_0.1.tsv"
OUT_P5E3_TSV="$OUT_DIR/TWN_p_lt_5e-3.tsv"
OUT_P01_GZ="$OUT_P01_TSV.gz"
OUT_P5E3_GZ="$OUT_P5E3_TSV.gz"

mkdir -p "$OUT_DIR"

# Stream all output to terminal and runtime log
touch "$RUN_LOG"
exec > >(tee -a "$RUN_LOG") 2>&1

log "Starting TWN p-value filtering job"
log "SLURM_JOB_ID: ${SLURM_JOB_ID:-not_set}"
log "Host: $(hostname)"
log "Working directory: $(pwd)"

if [[ ! -f "$INPUT_FILE" ]]; then
  log "Error: input file not found: $INPUT_FILE"
  exit 1
fi

if ! command -v gzip >/dev/null 2>&1; then
  log "Error: gzip not found in PATH"
  exit 1
fi

log "Project root: $PROJECT_ROOT"
log "Input file: $INPUT_FILE"
log "Output dir: $OUT_DIR"
log "Counting input rows..."
input_rows=$(wc -l < "$INPUT_FILE")
log "Input rows (including header): $input_rows"
log "Filtering with thresholds p<0.1 and p<5e-3"

awk -F'\t' -v OFS='\t' \
  -v out01="$OUT_P01_TSV" \
  -v out5e3="$OUT_P5E3_TSV" '
NR==1 {
  pcol = 0
  for (i = 1; i <= NF; i++) {
    if ($i == "p.value") {
      pcol = i
      break
    }
  }
  if (pcol == 0) {
    print "Error: required column p.value not found in header" > "/dev/stderr"
    exit 2
  }
  print $0 > out01
  print $0 > out5e3
  next
}
{
  if (NR % 1000000 == 0) {
    print "[progress] processed rows: " NR > "/dev/stderr"
  }
  p = $pcol + 0
  if (p < 0.1) {
    print $0 >> out01
  }
  if (p < 5e-3) {
    print $0 >> out5e3
  }
}
' "$INPUT_FILE" 2>> "$RUN_LOG"

log "Compression started"

gzip -c "$OUT_P01_TSV" > "$OUT_P01_GZ"
gzip -c "$OUT_P5E3_TSV" > "$OUT_P5E3_GZ"

rows_p01=$(( $(wc -l < "$OUT_P01_TSV") - 1 ))
rows_p5e3=$(( $(wc -l < "$OUT_P5E3_TSV") - 1 ))

log "Rows p<0.1: $rows_p01"
log "Rows p<5e-3: $rows_p5e3"
log "Finished successfully"

echo "Done. Created files:"
echo "  - $OUT_P01_TSV"
echo "  - $OUT_P01_GZ"
echo "  - $OUT_P5E3_TSV"
echo "  - $OUT_P5E3_GZ"
echo "  - $RUN_LOG"
