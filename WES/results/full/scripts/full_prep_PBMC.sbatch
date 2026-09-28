#!/bin/bash
#SBATCH --job-name=full_prep_PBMC
#SBATCH --partition=cn-long
#SBATCH --account=bjx131_g1
#SBATCH --qos=bjx131cnl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=0
#SBATCH --time=16:00:00
#SBATCH --output=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.out
#SBATCH --error=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920/log/%x_%j.err
#
# FULL-data round, PBMC sample prep: roadmap sections 2-6.2 kept verbatim
# EXCEPT the agreed deviation: bwa mem | samtools sort pipe replaces
# roadmap 3.1+3.2+3.3 (SAM -> BAM -> sort). Final blood_sorted.bam is identical.

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
export DATA=$SHARE/data/OC_WES_all
export BASE=/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920
export OUT=$BASE

source "$SHARE/miniconda3/etc/profile.d/conda.sh"
conda activate wes
set -uo pipefail

DEVIATIONS="$OUT/log/deviations.log"
mkdir -p "$OUT"/{0_fastq/PBMC,1_alignment/PBMC,2_variant_call,3_annotation,log,resources}
touch "$DEVIATIONS"
note()  { echo "[DEVIATION $(date '+%F %T')][PBMC] $*" | tee -a "$DEVIATIONS"; }
stage() { echo "==== [$(date '+%F %T')][PBMC] STAGE $* ===="; }
die()   { echo "[FATAL][PBMC] $*" | tee -a "$DEVIATIONS" >&2; exit 1; }

# ------------------------------------------------------------------ 0.2 symlinks (roadmap 0.2, verbatim for PBMC)
ln -sf "$DATA/PBMC/blood_R1.fq.gz" "$OUT/0_fastq/PBMC_R1.fq.gz"
ln -sf "$DATA/PBMC/blood_R2.fq.gz" "$OUT/0_fastq/PBMC_R2.fq.gz"

# ------------------------------------------------------------------ 2.2 fastp (roadmap 2.2, verbatim)
stage "2.2 fastp"
fastp \
  -i "$OUT/0_fastq/PBMC_R1.fq.gz" \
  -I "$OUT/0_fastq/PBMC_R2.fq.gz" \
  -o "$OUT/0_fastq/PBMC_R1.clean.fq.gz" \
  -O "$OUT/0_fastq/PBMC_R2.clean.fq.gz" \
  --detect_adapter_for_pe \
  --cut_right \
  --cut_right_window_size 4 \
  --cut_right_mean_quality 20 \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --n_base_limit 5 \
  --length_required 50 \
  --thread 8 \
  --html "$OUT/0_fastq/PBMC/blood.fastp.html" \
  --json "$OUT/0_fastq/PBMC/blood.fastp.json" || die "2.2 fastp failed"

# ------------------------------------------------------------------ 3. Alignment
stage "3.1-3.3 bwa mem | samtools sort  (DEVIATION: piped, replaces SAM/BAM/sort steps)"
note "3.1-3.3 alignment piped: 'bwa mem | samtools sort' replaces roadmap literal SAM->BAM->sort (agreed efficiency deviation; final sorted BAM identical)"
bwa mem \
  -t 8 \
  -M \
  -R "@RG\tID:PBMC\tPL:illumina\tSM:blood\tLB:WES_PBMC" \
  "$REFERENCE" \
  "$OUT/0_fastq/PBMC_R1.clean.fq.gz" \
  "$OUT/0_fastq/PBMC_R2.clean.fq.gz" \
  | samtools sort -@ 8 -m 4G \
    -o "$OUT/1_alignment/PBMC/blood_sorted.bam" \
    - || die "3.x PBMC alignment/sort pipe failed"

samtools index -@ 8 "$OUT/1_alignment/PBMC/blood_sorted.bam"

stage "3.4 flagstat"
samtools flagstat -@ 8 "$OUT/1_alignment/PBMC/blood_sorted.bam" \
  > "$OUT/1_alignment/PBMC/blood.flagstat.txt"

# ------------------------------------------------------------------ 4. Mark PCR duplicates (roadmap 4.1, verbatim)
stage "4.1 MarkDuplicates"
gatk --java-options "-Xmx8g" MarkDuplicates \
  -I "$OUT/1_alignment/PBMC/blood_sorted.bam" \
  -O "$OUT/1_alignment/PBMC/blood_sorted.markdup.bam" \
  -M "$OUT/1_alignment/PBMC/blood_markdup_metrics.txt" \
  --CREATE_INDEX true || die "4.1 PBMC MarkDuplicates failed"

# ------------------------------------------------------------------ 5. BQSR (roadmap 5.1/5.2, verbatim)
stage "5.1 BaseRecalibrator"
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
  -I "$OUT/1_alignment/PBMC/blood_sorted.markdup.bam" \
  --bqsr-recal-file "$OUT/1_alignment/PBMC/blood.recal_data.table" \
  -O "$OUT/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam" || die "5.2 PBMC ApplyBQSR failed"

# ------------------------------------------------------------------ 6.2 HaplotypeCaller (roadmap 6.2, verbatim)
stage "6.2 HaplotypeCaller GVCF"
gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" HaplotypeCaller \
  -ERC GVCF \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam" \
  -D "$DBSNP" \
  -L "$TARGETS" \
  --native-pair-hmm-threads 8 \
  -O "$OUT/1_alignment/PBMC/blood.g.vcf" || die "6.2 PBMC HaplotypeCaller failed"

stage "PBMC prep done"
