#!/bin/bash
#SBATCH --job-name=full_joint
#SBATCH --partition=cn-long
#SBATCH --account=bjx131_g1
#SBATCH --qos=bjx131cnl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.out
#SBATCH --error=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.err
#
# FULL-data round: germline joint steps, roadmap sections 6.3-8 kept verbatim.
# Depends on full_prep_OC + full_prep_PBMC having produced BQSR BAMs and GVCFs.

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

mkdir -p "$OUT"/{2_variant_call,3_annotation,log}
DEVIATIONS="$OUT/log/deviations.log"
touch "$DEVIATIONS"
note()  { echo "[DEVIATION $(date '+%F %T')][joint] $*" | tee -a "$DEVIATIONS"; }
stage() { echo "==== [$(date '+%F %T')][joint] STAGE $* ===="; }
die()   { echo "[FATAL][joint] $*" | tee -a "$DEVIATIONS" >&2; exit 1; }

# ------------------------------------------------------------------ 6.3 Combine and genotype GVCFs (roadmap 6.3, verbatim)
stage "6.3 CombineGVCFs + GenotypeGVCFs"
gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" CombineGVCFs \
  -R "$REFERENCE" \
  -V "$OUT/1_alignment/OC/OC.g.vcf" \
  -V "$OUT/1_alignment/PBMC/blood.g.vcf" \
  -O "$OUT/2_variant_call/combined.g.vcf" || die "6.3 CombineGVCFs failed"

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" GenotypeGVCFs \
  -R "$REFERENCE" \
  -V "$OUT/2_variant_call/combined.g.vcf" \
  -D "$DBSNP" \
  -O "$OUT/2_variant_call/OC_blood_variants.vcf" || die "6.3 GenotypeGVCFs failed"

# ------------------------------------------------------------------ 6.4 VariantAnnotator (roadmap 6.4, verbatim)
stage "6.4 VariantAnnotator"
gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" VariantAnnotator \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam" \
  -V "$OUT/2_variant_call/OC_blood_variants.vcf" \
  -D "$DBSNP" \
  -O "$OUT/2_variant_call/OC_blood_variants.anno.vcf" \
  || note "6.4 VariantAnnotator failed (roadmap notes it may be redundant in GATK4)"

# ------------------------------------------------------------------ 7. VQSR (roadmap 7.1/7.2/7.3, verbatim)
stage "7.1 SNP VQSR"
VQSR_SNP_OK=no
if gatk --java-options "-Xmx20g -XX:ParallelGCThreads=4" VariantRecalibrator \
  -R "$REFERENCE" \
  -V "$OUT/2_variant_call/OC_blood_variants.vcf" \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "$HAPMAP" \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 "$OMNI" \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 "$HIGHCONF" \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -mode SNP \
  -O "$OUT/2_variant_call/OC_blood.snps.recal" \
  --tranches-file "$OUT/2_variant_call/OC_blood.snps.tranches" \
  --rscript-file "$OUT/2_variant_call/OC_blood.snps.plots.R" \
  --dont-run-rscript; then

  if gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" ApplyVQSR \
    -R "$REFERENCE" \
    -V "$OUT/2_variant_call/OC_blood_variants.vcf" \
    -O "$OUT/2_variant_call/OC_blood.snps.VQSR.vcf" \
    --truth-sensitivity-filter-level 99.5 \
    --tranches-file "$OUT/2_variant_call/OC_blood.snps.tranches" \
    --recal-file "$OUT/2_variant_call/OC_blood.snps.recal" \
    --mode SNP; then VQSR_SNP_OK=yes; fi
else
  note "7.1 SNP VariantRecalibrator failed"
fi

stage "7.2 INDEL VQSR"
VQSR_INDEL_OK=no
if gatk --java-options "-Xmx20g -XX:ParallelGCThreads=4" VariantRecalibrator \
  -R "$REFERENCE" \
  -V "$OUT/2_variant_call/OC_blood_variants.vcf" \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 "$INDELS" \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
  -an QD -an ReadPosRankSum -an FS -an SOR \
  -mode INDEL \
  --max-gaussians 4 \
  -O "$OUT/2_variant_call/OC_blood.indel.recal" \
  --tranches-file "$OUT/2_variant_call/OC_blood.indel.tranches" \
  --rscript-file "$OUT/2_variant_call/OC_blood.indel.plots.R" \
  --dont-run-rscript; then

  if gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" ApplyVQSR \
    -R "$REFERENCE" \
    -V "$OUT/2_variant_call/OC_blood_variants.vcf" \
    -O "$OUT/2_variant_call/OC_blood.indel.VQSR.vcf" \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file "$OUT/2_variant_call/OC_blood.indel.tranches" \
    --recal-file "$OUT/2_variant_call/OC_blood.indel.recal" \
    --mode INDEL; then VQSR_INDEL_OK=yes; fi
else
  note "7.2 INDEL VariantRecalibrator failed"
fi

stage "7.3 MergeVcfs"
VQSR_MERGED=no
if [[ "$VQSR_SNP_OK" == yes && "$VQSR_INDEL_OK" == yes ]]; then
  gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" MergeVcfs \
    -O "$OUT/2_variant_call/OC_blood.all.VQSR.vcf" \
    -I "$OUT/2_variant_call/OC_blood.indel.VQSR.vcf" \
    -I "$OUT/2_variant_call/OC_blood.snps.VQSR.vcf" && VQSR_MERGED=yes
  printf "PASS variants in merged VCF: "
  bcftools view -H -f PASS "$OUT/2_variant_call/OC_blood.all.VQSR.vcf" | wc -l
else
  note "7.3 MergeVcfs skipped: SNP ok=$VQSR_SNP_OK INDEL ok=$VQSR_INDEL_OK"
fi

# ------------------------------------------------------------------ 8. ANNOVAR (roadmap 8.1-8.3, verbatim)
stage "8.1 convert2annovar"
mkdir -p "$OUT/3_annotation"
ANNO_IN="$OUT/2_variant_call/OC_blood.all.VQSR.vcf"
if [[ ! -s "$ANNO_IN" ]]; then
  note "8.x fallback: VQSR merged VCF missing -> annotating OC_blood_variants.anno.vcf (pre-VQSR joint callset) instead"
  ANNO_IN="$OUT/2_variant_call/OC_blood_variants.anno.vcf"
  [[ -s "$ANNO_IN" ]] || ANNO_IN="$OUT/2_variant_call/OC_blood_variants.vcf"
fi
echo "ANNOVAR input: $ANNO_IN"
perl "$ANNOVAR/convert2annovar.pl" \
  -format vcf4old \
  "$ANNO_IN" \
  -includeinfo \
  -comment \
  -outfile "$OUT/3_annotation/OC_blood.all.avinput" \
  || die "8.1 convert2annovar failed even with fallback input"

stage "8.2 gene annotation"
perl "$ANNOVAR/annotate_variation.pl" \
  -geneanno \
  -buildver hg38 \
  -dbtype refGeneWithVer \
  "$OUT/3_annotation/OC_blood.all.avinput" \
  "$HUMANDB" || die "8.2 annotate_variation failed"

stage "8.3 table_annovar (basic protocol)"
perl "$ANNOVAR/table_annovar.pl" \
  "$OUT/3_annotation/OC_blood.all.avinput" \
  "$HUMANDB" \
  -buildver hg38 \
  -out "$OUT/3_annotation/OC_blood_anno" \
  -remove \
  -protocol refGeneWithVer \
  -operation g || die "8.3 table_annovar basic failed"

stage "8.3b table_annovar (richer protocol, verbatim - names not in humandb)"
if perl "$ANNOVAR/table_annovar.pl" \
  "$OUT/3_annotation/OC_blood.all.avinput" \
  "$HUMANDB" \
  -buildver hg38 \
  -out "$OUT/3_annotation/OC_blood_anno_rich" \
  -remove \
  -protocol refGeneWithVer,cytoBand,dbnsfp30a,exac03,gnomad_genome,avsnp147 \
  -operation g,r,f,f,f,f; then
  note "8.3b rich protocol succeeded (unexpected - humandb names matched)"
else
  note "8.3b rich protocol failed as anticipated: humandb lacks dbnsfp30a/exac03/gnomad_genome/avsnp147 (installed: dbnsfp54a, gnomad211_exome, gnomad41_exome, avsnp151); basic protocol output stands"
fi

# ------------------------------------------------------------------ summary (addition for step-11 checklist)
stage "joint summary"
SUM="$OUT/log/summary_joint.log"
{
  echo "=== 6.3 joint callset counts ==="
  printf "OC_blood_variants.vcf records: "; bcftools view -H "$OUT/2_variant_call/OC_blood_variants.vcf" | wc -l
  if [[ -s "$OUT/2_variant_call/OC_blood_variants.anno.vcf" ]]; then
    printf "anno.vcf records: "; bcftools view -H "$OUT/2_variant_call/OC_blood_variants.anno.vcf" | wc -l
  fi
  echo "=== 7 VQSR status ==="; echo "SNP ok=$VQSR_SNP_OK INDEL ok=$VQSR_INDEL_OK merged=$VQSR_MERGED"
  if [[ "$VQSR_MERGED" == yes ]]; then
    printf "PASS in merged VQSR vcf: "; bcftools view -H -f PASS "$OUT/2_variant_call/OC_blood.all.VQSR.vcf" | wc -l
    echo "--- SNP tranche tail ---"; tail -3 "$OUT/2_variant_call/OC_blood.snps.tranches" 2>/dev/null
    echo "--- INDEL tranche tail ---"; tail -3 "$OUT/2_variant_call/OC_blood.indel.tranches" 2>/dev/null
    echo "--- Ti/Tv (bcftools stats) ---"
    bcftools stats -f PASS "$OUT/2_variant_call/OC_blood.all.VQSR.vcf" 2>/dev/null | grep -E "^TSTV|^SN" | head -8
  fi
  echo "=== 8 ANNOVAR ==="
  printf "avinput lines: "; grep -vc '^#' "$OUT/3_annotation/OC_blood.all.avinput" || true
  printf "multianno rows: "; awk 'END{print NR}' "$OUT/3_annotation/OC_blood_anno.hg38_multianno.txt"
  echo "multianno columns:"; head -1 "$OUT/3_annotation/OC_blood_anno.hg38_multianno.txt" | tr '\t' '|'
} > "$SUM" 2>&1
echo "summary written to $SUM"
stage "joint done"
