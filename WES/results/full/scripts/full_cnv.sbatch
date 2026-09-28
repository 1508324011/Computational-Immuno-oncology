#!/bin/bash
#SBATCH --job-name=full_cnv
#SBATCH --partition=cn-long
#SBATCH --account=bjx131_g1
#SBATCH --qos=bjx131cnl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --output=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.out
#SBATCH --error=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.err
#
# FULL-data round: CNVkit, roadmap section 10 kept verbatim (amplicon mode).
# Depends on full_prep_OC + full_prep_PBMC having produced BQSR BAMs.

# ------------------------------------------------------------------ 0. Environment (roadmap 0)
export SHARE=/lustre1/share
export REF=$SHARE/references
export TARGETS=$REF/S07604514_AllTracks_V6_60_hg38.bed
export REFERENCE=$REF/hg38.fa
export BASE=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920
export OUT=$BASE

source "$SHARE/miniconda3/etc/profile.d/conda.sh"
conda activate wes
set -uo pipefail

mkdir -p "$OUT"/{4_cnv,log}
DEVIATIONS="$OUT/log/deviations.log"
touch "$DEVIATIONS"
note()  { echo "[DEVIATION $(date '+%F %T')][cnv] $*" | tee -a "$DEVIATIONS"; }
stage() { echo "==== [$(date '+%F %T')][cnv] STAGE $* ===="; }
die()   { echo "[FATAL][cnv] $*" | tee -a "$DEVIATIONS" >&2; exit 1; }

TUMOR=$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam
NORMAL=$OUT/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam

# ------------------------------------------------------------------ 10.1 Prepare target BED (roadmap 10.1, verbatim)
stage "10.1 prepare target BED"
CNVOUT=$OUT/4_cnv
mkdir -p "$CNVOUT"
ln -sf "$TARGETS" "$CNVOUT/targets.bed"

# ------------------------------------------------------------------ 10.2 CNVkit batch (roadmap 10.2, verbatim: -m amplicon)
stage "10.2 cnvkit.py batch"
if cnvkit.py batch \
  "$TUMOR" \
  -n "$NORMAL" \
  -f "$REFERENCE" \
  -t "$CNVOUT/targets.bed" \
  -m amplicon \
  --output-reference "$CNVOUT/reference.cnn" \
  -d "$CNVOUT/cnvkit_out" \
  --scatter \
  --diagram \
  -p 8; then
  echo "CNVkit batch done"
else
  note "10.2 cnvkit.py batch failed (amplicon mode on hybrid-capture WES panel; see job log)"
fi

# ------------------------------------------------------------------ summary
stage "cnv summary"
SUM="$OUT/log/summary_cnv.log"
{
  echo "=== 10 CNVkit (full data) ==="
  ls -lh "$CNVOUT/cnvkit_out" 2>/dev/null
  for f in "$CNVOUT"/cnvkit_out/*.cnr "$CNVOUT"/cnvkit_out/*.cns "$CNVOUT"/cnvkit_out/*.call.cns; do
    [[ -e "$f" ]] || continue
    printf "%s rows: %s\n" "$(basename "$f")" "$(wc -l < "$f")"
  done
  CNS=$(ls "$CNVOUT"/cnvkit_out/*.call.cns 2>/dev/null | head -1)
  if [[ -n "${CNS:-}" && -s "$CNS" ]]; then
    echo "--- segments by called cn ---"
    awk -F'\t' 'NR>1 {c[$6]++} END {for (k in c) print "cn="k": "c[k]" segments"}' "$CNS"
    echo "--- segments with |log2| >= 0.3 ---"
    awk -F'\t' 'NR>1 && ($5<=-0.3 || $5>=0.3)' "$CNS" | head -30
  fi
} > "$SUM" 2>&1
echo "summary written to $SUM"
stage "cnv done"
