#!/bin/bash
#SBATCH --job-name=full_somatic
#SBATCH --partition=cn-long
#SBATCH --account=bjx131_g1
#SBATCH --qos=bjx131cnl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=16:00:00
#SBATCH --output=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.out
#SBATCH --error=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.err
#
# FULL-data round: somatic TMB, roadmap section 9 kept verbatim.
# Depends on full_prep_OC + full_prep_PBMC having produced BQSR BAMs.

# ------------------------------------------------------------------ 0. Environment (roadmap 0)
export SHARE=/lustre1/share
export REF=$SHARE/references
export TARGETS=$REF/S07604514_AllTracks_V6_60_hg38.bed
export REFERENCE=$REF/hg38.fa
export DBSNP=$REF/dbsnp_138.hg38.vcf
export INDELS=$REF/Mills_and_1000G_gold_standard.indels.hg38.vcf
export HIGHCONF=$REF/1000G_phase1.snps.high_confidence.hg38.vcf
export HAPMAP=$REF/hapmap_3.3.hg38.vcf
export OMNI=$REF/1000G_omni2.5.hg38.vcf
export ANNOVAR=$SHARE/annovar_new
export HUMANDB=$ANNOVAR/humandb
export PATH=$ANNOVAR:$PATH
export BASE=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920
export OUT=$BASE

source "$SHARE/miniconda3/etc/profile.d/conda.sh"
conda activate wes
set -uo pipefail

mkdir -p "$OUT"/{5_somatic_tmb,log}
DEVIATIONS="$OUT/log/deviations.log"
touch "$DEVIATIONS"
note()  { echo "[DEVIATION $(date '+%F %T')][somatic] $*" | tee -a "$DEVIATIONS"; }
stage() { echo "==== [$(date '+%F %T')][somatic] STAGE $* ===="; }
die()   { echo "[FATAL][somatic] $*" | tee -a "$DEVIATIONS" >&2; exit 1; }

# ------------------------------------------------------------------ 9.1 Mutect2 (roadmap 9.1, verbatim)
stage "9.1 Mutect2"
TUMOR=$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam
NORMAL=$OUT/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam
SOMOUT=$OUT/5_somatic_tmb/

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" Mutect2 \
  -R "$REFERENCE" \
  -I "$TUMOR" \
  -I "$NORMAL" \
  --normal-sample blood \
  -L "$TARGETS" \
  --germline-resource "$HIGHCONF" \
  -O "$SOMOUT/somatic.raw.vcf.gz" \
  --native-pair-hmm-threads 8 \
  --f1r2-tar-gz "$SOMOUT/f1r2.tar.gz" || die "9.1 Mutect2 failed"

# ------------------------------------------------------------------ 9.2 orientation bias (roadmap 9.2, verbatim)
stage "9.2 LearnReadOrientationModel + FilterMutectCalls"
gatk --java-options "-Xmx8g -XX:ParallelGCThreads=4" LearnReadOrientationModel \
  -I "$SOMOUT/f1r2.tar.gz" \
  -O "$SOMOUT/orientation_model.tar.gz" || die "9.2 LearnReadOrientationModel failed"

gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" FilterMutectCalls \
  -R "$REFERENCE" \
  -V "$SOMOUT/somatic.raw.vcf.gz" \
  --stats "$SOMOUT/somatic.raw.vcf.gz.stats" \
  --ob-priors "$SOMOUT/orientation_model.tar.gz" \
  -O "$SOMOUT/somatic.filtered.vcf.gz" || die "9.2 FilterMutectCalls failed"

bcftools view -f PASS \
  -Oz \
  -o "$SOMOUT/somatic.pass.vcf.gz" \
  "$SOMOUT/somatic.filtered.vcf.gz"

bcftools index -t "$SOMOUT/somatic.pass.vcf.gz"

# ------------------------------------------------------------------ 9.3 TMB (roadmap 9.3, verbatim)
stage "9.3 TMB"
TARGET_SIZE_MB=$(awk 'BEGIN{FS="\t"} {sum += ($3 - $2)} END {printf "%.6f", sum / 1000000}' "$TARGETS")
PASS_VARIANTS=$(bcftools view -H -f PASS "$SOMOUT/somatic.filtered.vcf.gz" | wc -l)

awk -v n="$PASS_VARIANTS" -v mb="$TARGET_SIZE_MB" \
  'BEGIN {printf "somatic TMB = %.6f variants/Mb\n", n / mb}'

# ------------------------------------------------------------------ summary
stage "somatic summary"
SUM="$OUT/log/summary_somatic.log"
{
  echo "=== 9 somatic (full data) ==="
  printf "raw records: ";      bcftools view -H "$SOMOUT/somatic.raw.vcf.gz" | wc -l
  printf "filtered records: ";  bcftools view -H "$SOMOUT/somatic.filtered.vcf.gz" | wc -l
  printf "PASS records: ";     bcftools view -H -f PASS "$SOMOUT/somatic.filtered.vcf.gz" | wc -l
  echo "filter breakdown:"
  bcftools query -f '%FILTER\n' "$SOMOUT/somatic.filtered.vcf.gz" 2>/dev/null | sort | uniq -c | sort -rn
  printf "target Mb: %s\n" "$TARGET_SIZE_MB"
  awk -v n="$PASS_VARIANTS" -v mb="$TARGET_SIZE_MB" 'BEGIN {printf "somatic TMB = %.6f variants/Mb\n", n / mb}'
  echo "PASS variants list (chr:pos ref>alt, first 50):"
  bcftools query -f '%CHROM:%POS %REF>%ALT\n' -f PASS "$SOMOUT/somatic.filtered.vcf.gz" 2>/dev/null | head -50
} > "$SUM" 2>&1
echo "summary written to $SUM"
stage "somatic done"
