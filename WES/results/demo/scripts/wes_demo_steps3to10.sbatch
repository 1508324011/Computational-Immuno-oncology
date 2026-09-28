#!/bin/bash
#SBATCH --job-name=wes_demo
#SBATCH --partition=cn-long
#SBATCH --account=bjx131_g1
#SBATCH --qos=bjx131cnl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/lustre1/user/bjx131_pkuhpc/wes_run_cluster_20260914/log/%x_%j.out
#SBATCH --error=/lustre1/user/bjx131_pkuhpc/wes_run_cluster_20260914/log/%x_%j.err
#
# Demo round on downsampled data_new — roadmap commands kept verbatim (sections 3-10).
# Only additions: SLURM header, env activation, expected-failure recording (deviations.log),
# and a final summary block for the step-11 checklist. fastp (section 2) was already done.

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
export DATA=$SHARE/data_new
export BASE=/lustre1/user/bjx131_pkuhpc/wes_run_cluster_20260914   # continue previous run (roadmap 0 note)
export OUT=$BASE

source "$SHARE/miniconda3/etc/profile.d/conda.sh"
conda activate wes

set -uo pipefail   # (after conda activation: conda's activate scripts are not nounset-safe)

mkdir -p "$OUT"/{0_fastq,1_alignment/OC,1_alignment/PBMC,2_variant_call,3_annotation,4_cnv,5_somatic_tmb}
mkdir -p "$OUT"/{log,resources}

DEVIATIONS="$OUT/log/deviations.log"
: > "$DEVIATIONS"
note()  { echo "[DEVIATION $(date '+%F %T')] $*" | tee -a "$DEVIATIONS"; }
stage() { echo "==== [$(date '+%F %T')] STAGE $* ===="; }
die()   { echo "[FATAL] $*" | tee -a "$DEVIATIONS" >&2; exit 1; }

# ------------------------------------------------------------------ preflight
stage "preflight"
for f in "$REFERENCE"{,.fai,.amb,.ann,.bwt,.pac,.sa} "$REF/hg38.dict" "$TARGETS" \
        "$DBSNP"{,.idx} "$INDELS"{,.idx} "$HIGHCONF"{,.idx} "$HAPMAP"{,.idx} "$OMNI"{,.idx} \
        "$OUT/0_fastq/OC_R1.clean.fq.gz" "$OUT/0_fastq/OC_R2.clean.fq.gz" \
        "$OUT/0_fastq/PBMC_R1.clean.fq.gz" "$OUT/0_fastq/PBMC_R2.clean.fq.gz" \
        "$HUMANDB/hg38_refGeneWithVer.txt" "$HUMANDB/hg38_refGeneWithVerMrna.fa"; do
  [[ -e "$f" ]] || die "preflight: missing $f"
done
echo "preflight OK: node $(hostname), $(nproc) CPUs"

# ------------------------------------------------------------------ 3. Alignment with BWA-MEM
stage "3.1 bwa mem"
bwa mem \
  -t 8 \
  -M \
  -R "@RG\tID:OC\tPL:illumina\tSM:OC\tLB:WES_OC" \
  "$REFERENCE" \
  "$OUT/0_fastq/OC_R1.clean.fq.gz" \
  "$OUT/0_fastq/OC_R2.clean.fq.gz" \
  > "$OUT/1_alignment/OC/OC.sam" || die "3.1 OC alignment failed"

bwa mem \
  -t 8 \
  -M \
  -R "@RG\tID:PBMC\tPL:illumina\tSM:blood\tLB:WES_PBMC" \
  "$REFERENCE" \
  "$OUT/0_fastq/PBMC_R1.clean.fq.gz" \
  "$OUT/0_fastq/PBMC_R2.clean.fq.gz" \
  > "$OUT/1_alignment/PBMC/blood.sam" || die "3.1 PBMC alignment failed"

stage "3.2 SAM to BAM"
samtools view -bS -@ 8 \
  "$OUT/1_alignment/OC/OC.sam" \
  > "$OUT/1_alignment/OC/OC.bam"

samtools view -bS -@ 8 \
  "$OUT/1_alignment/PBMC/blood.sam" \
  > "$OUT/1_alignment/PBMC/blood.bam"

stage "3.3 sort and index"
samtools sort -@ 8 -m 4G \
  "$OUT/1_alignment/OC/OC.bam" \
  -o "$OUT/1_alignment/OC/OC_sorted.bam"

samtools sort -@ 8 -m 4G \
  "$OUT/1_alignment/PBMC/blood.bam" \
  -o "$OUT/1_alignment/PBMC/blood_sorted.bam"

samtools index -@ 8 "$OUT/1_alignment/OC/OC_sorted.bam"
samtools index -@ 8 "$OUT/1_alignment/PBMC/blood_sorted.bam"

stage "3.4 flagstat"
samtools flagstat -@ 8 "$OUT/1_alignment/OC/OC_sorted.bam" \
  > "$OUT/1_alignment/OC/OC.flagstat.txt"

samtools flagstat -@ 8 "$OUT/1_alignment/PBMC/blood_sorted.bam" \
  > "$OUT/1_alignment/PBMC/blood.flagstat.txt"

[[ -s "$OUT/1_alignment/OC/OC_sorted.bam" ]] || die "3.x OC_sorted.bam missing"
[[ -s "$OUT/1_alignment/PBMC/blood_sorted.bam" ]] || die "3.x blood_sorted.bam missing"

# ------------------------------------------------------------------ 4. Mark PCR duplicates
stage "4.1 MarkDuplicates"
gatk --java-options "-Xmx8g" MarkDuplicates \
  -I "$OUT/1_alignment/OC/OC_sorted.bam" \
  -O "$OUT/1_alignment/OC/OC_sorted.markdup.bam" \
  -M "$OUT/1_alignment/OC/OC_markdup_metrics.txt" \
  --CREATE_INDEX true || die "4.1 OC MarkDuplicates failed"

gatk --java-options "-Xmx8g" MarkDuplicates \
  -I "$OUT/1_alignment/PBMC/blood_sorted.bam" \
  -O "$OUT/1_alignment/PBMC/blood_sorted.markdup.bam" \
  -M "$OUT/1_alignment/PBMC/blood_markdup_metrics.txt" \
  --CREATE_INDEX true || die "4.1 PBMC MarkDuplicates failed"

# ------------------------------------------------------------------ 5. BQSR
stage "5.1 BaseRecalibrator"
gatk --java-options "-Xmx20g -XX:ParallelGCThreads=4" BaseRecalibrator \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/OC/OC_sorted.markdup.bam" \
  --known-sites "$DBSNP" \
  --known-sites "$INDELS" \
  --known-sites "$HIGHCONF" \
  -O "$OUT/1_alignment/OC/OC.recal_data.table" || die "5.1 OC BaseRecalibrator failed"

gatk --java-options "-Xmx20g -XX:ParallelGCThreads=4" BaseRecalibrator \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/PBMC/blood_sorted.markdup.bam" \
  --known-sites "$DBSNP" \
  --known-sites "$INDELS" \
  --known-sites "$HIGHCONF" \
  -O "$OUT/1_alignment/PBMC/blood.recal_data.table" || die "5.1 PBMC BaseRecalibrator failed"

stage "5.2 ApplyBQSR"
gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" ApplyBQSR \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/OC/OC_sorted.markdup.bam" \
  --bqsr-recal-file "$OUT/1_alignment/OC/OC.recal_data.table" \
  -O "$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam" || die "5.2 OC ApplyBQSR failed"

gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" ApplyBQSR \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/PBMC/blood_sorted.markdup.bam" \
  --bqsr-recal-file "$OUT/1_alignment/PBMC/blood.recal_data.table" \
  -O "$OUT/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam" || die "5.2 PBMC ApplyBQSR failed"

# ------------------------------------------------------------------ 6. Germline variant calling
stage "6.2 HaplotypeCaller GVCF"
gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" HaplotypeCaller \
  -ERC GVCF \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam" \
  -D "$DBSNP" \
  -L "$TARGETS" \
  --native-pair-hmm-threads 8 \
  -O "$OUT/1_alignment/OC/OC.g.vcf" || die "6.2 OC HaplotypeCaller failed"

gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" HaplotypeCaller \
  -ERC GVCF \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam" \
  -D "$DBSNP" \
  -L "$TARGETS" \
  --native-pair-hmm-threads 8 \
  -O "$OUT/1_alignment/PBMC/blood.g.vcf" || die "6.2 PBMC HaplotypeCaller failed"

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

stage "6.4 VariantAnnotator"
gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" VariantAnnotator \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam" \
  -V "$OUT/2_variant_call/OC_blood_variants.vcf" \
  -D "$DBSNP" \
  -O "$OUT/2_variant_call/OC_blood_variants.anno.vcf" \
  || note "6.4 VariantAnnotator failed (roadmap notes it may be redundant in GATK4)"

# ------------------------------------------------------------------ 7. VQSR  (EXPECTED TO FAIL on downsampled data)
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
  note "7.1 SNP VariantRecalibrator failed (expected on ~0.5x downsampled data: too few variants to fit Gaussian mixture model)"
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
  note "7.2 INDEL VariantRecalibrator failed (expected on ~0.5x downsampled data)"
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

# ------------------------------------------------------------------ 8. ANNOVAR
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

stage "8.3b table_annovar (richer protocol, verbatim - EXPECTED TO FAIL)"
# roadmap prints this as an option "if frequency and prediction databases are available";
# the names it uses are NOT present in the installed humandb (see deviations). Kept verbatim.
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

# ------------------------------------------------------------------ 9. Somatic TMB (Mutect2)
stage "9.1 Mutect2"
TUMOR=$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam
NORMAL=$OUT/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam
SOMOUT=$OUT/5_somatic_tmb/
mkdir -p "$SOMOUT"

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

stage "9.2 orientation bias model + FilterMutectCalls"
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

stage "9.3 TMB"
TARGET_SIZE_MB=$(awk 'BEGIN{FS="\t"} {sum += ($3 - $2)} END {printf "%.6f", sum / 1000000}' "$TARGETS")
PASS_VARIANTS=$(bcftools view -H -f PASS "$SOMOUT/somatic.filtered.vcf.gz" | wc -l)

awk -v n="$PASS_VARIANTS" -v mb="$TARGET_SIZE_MB" \
  'BEGIN {printf "somatic TMB = %.6f variants/Mb\n", n / mb}'

# ------------------------------------------------------------------ 10. CNVkit
stage "10.1 prepare target BED"
CNVOUT=$OUT/4_cnv
mkdir -p "$CNVOUT"
ln -sf "$TARGETS" "$CNVOUT/targets.bed"

stage "10.2 CNVkit batch (amplicon mode, verbatim)"
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

# ------------------------------------------------------------------ 11. Summary (step-11 checklist data; addition, not in roadmap)
stage "11 summary"
SUM="$OUT/log/summary_demo.log"
{
  echo "=== node/resources ==="; hostname; nproc
  echo "=== 3.4 flagstat (OC) ===";   cat "$OUT/1_alignment/OC/OC.flagstat.txt"
  echo "=== 3.4 flagstat (PBMC) ==="; cat "$OUT/1_alignment/PBMC/blood.flagstat.txt"
  echo "=== 4.2 duplication (OC) ===";   grep -A 6 '^LIBRARY' "$OUT/1_alignment/OC/OC_markdup_metrics.txt" | head -8
  echo "=== 4.2 duplication (PBMC) ==="; grep -A 6 '^LIBRARY' "$OUT/1_alignment/PBMC/blood_markdup_metrics.txt" | head -8
  echo "=== 6.3 joint callset counts ==="
  printf "OC_blood_variants.vcf records: "; bcftools view -H "$OUT/2_variant_call/OC_blood_variants.vcf" | wc -l
  if [[ -s "$OUT/2_variant_call/OC_blood_variants.anno.vcf" ]]; then
    printf "anno.vcf records: "; bcftools view -H "$OUT/2_variant_call/OC_blood_variants.anno.vcf" | wc -l
  fi
  echo "=== 7 VQSR status ==="; echo "SNP ok=$VQSR_SNP_OK INDEL ok=$VQSR_INDEL_OK merged=$VQSR_MERGED"
  if [[ "$VQSR_MERGED" == yes ]]; then
    printf "PASS in merged VQSR vcf: "; bcftools view -H -f PASS "$OUT/2_variant_call/OC_blood.all.VQSR.vcf" | wc -l
  fi
  echo "=== 8 ANNOVAR ==="
  printf "avinput lines: "; grep -vc '^#' "$OUT/3_annotation/OC_blood.all.avinput" || true
  printf "multianno rows: "; awk 'END{print NR}' "$OUT/3_annotation/OC_blood_anno.hg38_multianno.txt"
  echo "=== 9 somatic ==="
  printf "raw records: ";      bcftools view -H "$SOMOUT/somatic.raw.vcf.gz" | wc -l
  printf "filtered records: "; bcftools view -H "$SOMOUT/somatic.filtered.vcf.gz" | wc -l
  printf "PASS records: ";     bcftools view -H -f PASS "$SOMOUT/somatic.filtered.vcf.gz" | wc -l
  printf "target Mb: %s\n" "$TARGET_SIZE_MB"
  awk -v n="$PASS_VARIANTS" -v mb="$TARGET_SIZE_MB" 'BEGIN {printf "somatic TMB = %.6f variants/Mb\n", n / mb}'
  echo "=== 10 CNVkit ==="
  for f in "$CNVOUT"/cnvkit_out/*.cnr "$CNVOUT"/cnvkit_out/*.cns "$CNVOUT"/cnvkit_out/*.call.cns "$CNVOUT"/cnvkit_out/*.png "$CNVOUT"/cnvkit_out/*.pdf; do
    [[ -e "$f" ]] && { printf "%s rows: " "$(basename "$f")"; wc -l < "$f"; }
  done
  echo "=== deviations ==="; cat "$DEVIATIONS"
} > "$SUM" 2>&1
echo "summary written to $SUM"

stage "done"
echo "ALL DONE. Deviations logged in $DEVIATIONS"
