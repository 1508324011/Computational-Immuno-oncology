#!/bin/bash
#SBATCH --job-name=full_fix_vqsr
#SBATCH --partition=cn-long
#SBATCH --account=bjx131_g1
#SBATCH --qos=bjx131cnl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --output=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.out
#SBATCH --error=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.err
#
# CORRECTED-EDITION extra (issue #1 + #2): sequential ApplyVQSR (each variant
# filtered exactly once by its own model) + ANNOVAR with the INSTALLED rich
# protocol names. Uses the recal models produced by full_joint.sbatch.

export SHARE=/lustre1/share
export REF=$SHARE/references
export TARGETS=$REF/S07604514_AllTracks_V6_60_hg38.bed
export REFERENCE=$REF/hg38.fa
export DBSNP=$REF/dbsnp_138.hg38.vcf
export ANNOVAR=$SHARE/annovar_new
export HUMANDB=$ANNOVAR/humandb
export PATH=$ANNOVAR:$PATH
export BASE=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920
export OUT=$BASE

source "$SHARE/miniconda3/etc/profile.d/conda.sh"
conda activate wes
set -uo pipefail

mkdir -p "$OUT"/{2_variant_call,3_annotation,log}
note()  { echo "[FIX-VQSR $(date '+%F %T')] $*"; }
die()   { echo "[FATAL][fix-vqsr] $*" >&2; exit 1; }

V=$OUT/2_variant_call
RAW="$V/OC_blood_variants.vcf"
SNP_RECAL="$V/OC_blood.snps.recal";      SNP_TR="$V/OC_blood.snps.tranches"
IND_RECAL="$V/OC_blood.indel.recal";     IND_TR="$V/OC_blood.indel.tranches"

note "sequential ApplyVQSR: step 1 (SNP model on full callset)"
gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" ApplyVQSR \
  -R "$REFERENCE" \
  -V "$RAW" \
  -O "$V/OC_blood.step1.snpVQSR.vcf" \
  --truth-sensitivity-filter-level 99.5 \
  --tranches-file "$SNP_TR" \
  --recal-file "$SNP_RECAL" \
  --mode SNP || die "ApplyVQSR SNP step failed"

note "sequential ApplyVQSR: step 2 (INDEL model on step-1 output)"
gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" ApplyVQSR \
  -R "$REFERENCE" \
  -V "$V/OC_blood.step1.snpVQSR.vcf" \
  -O "$V/OC_blood.all.sequentialVQSR.vcf" \
  --truth-sensitivity-filter-level 99.0 \
  --tranches-file "$IND_TR" \
  --recal-file "$IND_RECAL" \
  --mode INDEL || die "ApplyVQSR INDEL step failed"

FINAL="$V/OC_blood.all.sequentialVQSR.vcf"
printf "sequential-VQSR total records: "; bcftools view -H "$FINAL" | wc -l
printf "sequential-VQSR PASS:         "; bcftools view -H -f PASS "$FINAL" | wc -l
printf "unique sites:                  "; bcftools query -f '%CHROM:%POS\n' "$FINAL" | sort -u | wc -l

note "ANNOVAR on the corrected (sequential) VCF"
perl "$ANNOVAR/convert2annovar.pl" \
  -format vcf4old \
  "$FINAL" \
  -includeinfo \
  -comment \
  -outfile "$OUT/3_annotation/OC_blood.seqVQSR.avinput" || die "convert2annovar failed"

perl "$ANNOVAR/table_annovar.pl" \
  "$OUT/3_annotation/OC_blood.seqVQSR.avinput" \
  "$HUMANDB" \
  -buildver hg38 \
  -out "$OUT/3_annotation/OC_blood_anno_rich_installed" \
  -remove \
  -protocol refGeneWithVer,cytoBand,gnomad41_exome,dbnsfp54a,avsnp151 \
  -operation g,r,f,f,f || note "rich protocol with installed names failed"

note "summary"
SUM="$OUT/log/summary_fix_vqsr.log"
{
  echo "=== 7 CORRECTED: sequential ApplyVQSR (full data) ==="
  printf "total records: ";  bcftools view -H "$FINAL" | wc -l
  printf "PASS: ";            bcftools view -H -f PASS "$FINAL" | wc -l
  printf "unique sites: ";    bcftools query -f '%CHROM:%POS\n' "$FINAL" | sort -u | wc -l
  echo "--- Ti/Tv (bcftools stats) ---"
  bcftools stats -f PASS "$FINAL" 2>/dev/null | grep -E "^TSTV|^SN" | head -8
  echo "--- filter breakdown ---"
  bcftools query -f '%FILTER\n' "$FINAL" 2>/dev/null | sort | uniq -c | sort -rn | head -10
  echo "=== 8 CORRECTED: ANNOVAR rich (installed names) ==="
  RICH="$OUT/3_annotation/OC_blood_anno_rich_installed.hg38_multianno.txt"
  if [[ -s "$RICH" ]]; then
    printf "multianno rows: "; awk 'END{print NR}' "$RICH"
    echo "multianno columns:"; head -1 "$RICH" | tr '\t' '|'
    printf "AF<=0.01 rare variants (gnomad41 AF column approx): "
    awk -F'\t' 'NR>1 && $NF!="-" && $NF!="" && $NF+0<=0.01 {n++} END{print n+0}' "$RICH"
  else
    echo "rich multianno missing"
  fi
  printf "avinput lines: "; grep -vc '^#' "$OUT/3_annotation/OC_blood.seqVQSR.avinput" || true
} > "$SUM" 2>&1
echo "summary written to $SUM"
note "done"
