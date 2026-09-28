#!/bin/bash
#SBATCH --job-name=full_cnv_hybrid
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
# CORRECTED-EDITION extra (issue #4): CNVkit in HYBRID mode (correct for the
# SureSelect V6 hybrid-capture panel) with the cluster's access file.
# Companion to full_cnv.sbatch which runs the roadmap's literal amplicon mode.

export SHARE=/lustre1/share
export REF=$SHARE/references
export TARGETS=$REF/S07604514_AllTracks_V6_60_hg38.bed
export REFERENCE=$REF/hg38.fa
export ACCESS=$REF/access-5kb.hg38.bed
export BASE=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920
export OUT=$BASE

source "$SHARE/miniconda3/etc/profile.d/conda.sh"
conda activate wes
set -uo pipefail

mkdir -p "$OUT"/{4_cnv,log}
note()  { echo "[CNV-HYBRID $(date '+%F %T')] $*"; }
die()   { echo "[FATAL][cnv-hybrid] $*" >&2; exit 1; }

TUMOR=$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam
NORMAL=$OUT/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam
CNVOUT=$OUT/4_cnv
HYB=$CNVOUT/cnvkit_out_hybrid
mkdir -p "$HYB"

note "target bins + antitarget bins"
cnvkit.py target "$TARGETS" --annotate "$REF/refGene.hg38.txt" -o "$CNVOUT/targets.annotated.bed" \
  || note "target --annotate failed (no refGene txt at $REF); continuing with plain targets"
TARGET_BED="$CNVOUT/targets.annotated.bed"
[[ -s "$TARGET_BED" ]] || TARGET_BED="$TARGETS"

if [[ -s "$ACCESS" ]]; then
  cnvkit.py access "$REFERENCE" -o "$CNVOUT/access.bed" -g "$ACCESS" \
    || cp "$ACCESS" "$CNVOUT/access.bed"
else
  cnvkit.py access "$REFERENCE" -o "$CNVOUT/access.bed" || die "access generation failed"
fi

note "batch -m hybrid"
cnvkit.py batch \
  "$TUMOR" -n "$NORMAL" \
  -f "$REFERENCE" \
  -t "$TARGET_BED" \
  -g "$CNVOUT/access.bed" \
  -m hybrid \
  --output-reference "$CNVOUT/reference_hybrid.cnn" \
  -d "$HYB" \
  --scatter \
  --diagram \
  -p 8 || die "cnvkit batch hybrid failed"

note "summary"
SUM="$OUT/log/summary_cnv_hybrid.log"
{
  echo "=== 10 CNVkit CORRECTED (hybrid mode, access file) ==="
  ls -lh "$HYB"
  CNS=$(ls "$HYB"/*.call.cns 2>/dev/null | head -1)
  if [[ -n "$CNS" && -s "$CNS" ]]; then
    echo "--- segments by called cn ---"
    awk -F'\t' 'NR>1 {c[$6]++} END {for (k in c) print "cn="k": "c[k]" segments"}' "$CNS"
    echo "--- segments with |log2| >= 0.3 (top 30 by weight) ---"
    awk -F'\t' 'NR>1 && ($5<=-0.3 || $5>=0.3)' "$CNS" | sort -t$'\t' -k9 -gr | head -30
    echo "--- HGSOC classic genes grep (from diagram cns genes col if present) ---"
  fi
} > "$SUM" 2>&1
echo "summary written to $SUM"
note "done"
