#!/bin/bash
#SBATCH --job-name=full_hsmetrics
#SBATCH --partition=cn-long
#SBATCH --account=bjx131_g1
#SBATCH --qos=bjx131cnl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.out
#SBATCH --error=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.err
#
# CORRECTED-EDITION extra (issue #9): CollectHsMetrics for both rounds —
# the on-target depth numbers for the §11 checklist.
# Full-round BAMs (this run dir) + demo-round BAMs (wes_run_cluster_20260914).

export SHARE=/lustre1/share
export REF=$SHARE/references
export TARGETS=$REF/S07604514_AllTracks_V6_60_hg38.bed
export REFERENCE=$REF/hg38.fa
export BASE=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920
export OUT=$BASE
DEMO=/lustre1/user/bjx131_pkuhpc/wes_run_cluster_20260914

source "$SHARE/miniconda3/etc/profile.d/conda.sh"
conda activate wes
set -uo pipefail

mkdir -p "$OUT"/{log,qc}
note() { echo "[HS $(date '+%F %T')] $*"; }

# CollectHsMetrics (Picard) 要求带 SAM 头的 interval_list；先把 BED 转换过去
if [[ ! -s "$OUT/qc/targets.interval_list" ]]; then
  gatk --java-options "-Xmx4g" BedToIntervalList \
    -I "$TARGETS" \
    -O "$OUT/qc/targets.interval_list" \
    -SD "$REF/hg38.dict" || note "BedToIntervalList failed"
fi
IL="$OUT/qc/targets.interval_list"

run_hs () {  # $1=bam $2=label
  gatk --java-options "-Xmx8g" CollectHsMetrics \
    -I "$1" \
    -O "$OUT/qc/$2.hs_metrics.txt" \
    --TARGET_INTERVALS "$IL" \
    --BAIT_INTERVALS "$IL" \
    --REFERENCE_SEQUENCE "$REFERENCE" || note "CollectHsMetrics failed for $2"
  note "done $2"
}

run_hs "$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam"       "full_OC"
run_hs "$OUT/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam"  "full_PBMC"
run_hs "$DEMO/1_alignment/OC/OC_sorted.markdup.BQSR.bam"      "demo_OC"
run_hs "$DEMO/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam" "demo_PBMC"

# demo joint-VCF stats (ck9 het/hom) — tiny, safe here
bcftools stats "$DEMO/2_variant_call/OC_blood_variants.vcf" > "$OUT/qc/demo_joint_stats.txt" 2>&1 \
  || note "bcftools stats on demo joint vcf failed"

{
  echo "=== MEAN_TARGET_COVERAGE / PCT 30X / on-target % ==="
  for f in "$OUT"/qc/*.hs_metrics.txt; do
    printf "%s: " "$(basename "$f")"
    awk -F'\t' 'NR==5 {printf "meanCov=%.1f PCT_30X=%.4f PCT_100X=%.4f ON_TARGET=%.4f FOLD_ENRICH=%.1f\n", $13, $29, $33, $42, $44}' "$f"
  done
  echo "=== demo joint het/hom ==="
  grep -E "^PSC|^SN" "$OUT/qc/demo_joint_stats.txt" | head -20
} > "$OUT/log/summary_hsmetrics.log" 2>&1
echo "summary written"
