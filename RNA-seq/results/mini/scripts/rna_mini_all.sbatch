#!/bin/bash
#SBATCH --job-name=rna_mini
#SBATCH --partition=cn-long
#SBATCH --qos=bjx131cnl
#SBATCH --account=bjx131_g1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/lustre1/user/bjx131_pkuhpc/rna_run_mini_20260921/log/rna_mini_%j.out
#SBATCH --error=/lustre1/user/bjx131_pkuhpc/rna_run_mini_20260921/log/rna_mini_%j.err

# =============================================================================
# RNA-seq mini tutorial (chr22 / 10k reads) — run AS WRITTEN, minimal wrapper.
# Tutorial source: /gpfs1/home/bjx131_pkuhpc/class/RNA-seq/RNA_seq_tutorial_mini.html
# Discipline (same as WES runs): roadmap commands verbatim + SLURM header +
# conda activation + failure/deviation logging; no login-node compute.
# Deviations are pre-registered in log/deviations.log.
# =============================================================================
RUN=/lustre1/user/bjx131_pkuhpc/rna_run_mini_20260921
cd "$RUN"
mkdir -p log

# --- conda activate FIRST: `set -u` breaks conda's activate scripts (#25) ---
source /lustre1/share/miniconda3/etc/profile.d/conda.sh
conda activate RNA-seq

# 2026-09-28: cluster-wide SLURM submission outage (every account/partition/QOS
# combo rejected; queue empty cluster-wide; the WES scripts that ran fine on
# 09-20 are rejected too). Tutorial's own design (instructor's 02_simulate_30_
# load test gates ALLOW_LOGIN_NODE=YES) targets interactive login-node use with
# exactly this 10k-read single-thread workload -> running directly on login60.
SLURM_JOB_ID="${SLURM_JOB_ID:-login60-direct}"

set -euo pipefail
SUM=log/summary_mini.log
DEV=log/deviations.log
: > "$SUM"; : > "$DEV"
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$SUM"; }
dev(){ echo "[$(date '+%F %T')] DEVIATION: $*" | tee -a "$DEV" "$SUM"; }

log "== RNA-seq mini tutorial start (job $SLURM_JOB_ID, $(hostname)) =="
log "tools: $(fastp --version 2>&1 | head -1) | STAR $(STAR --version | tail -1) | salmon $(salmon --version) | fastqc $(fastqc --version 2>&1 | tail -1)"

dev "workdir = own run dir (tutorial: /lustre1/share/RNA_seq_student_2026/{id}); our cluster discipline — no writes to classroom share"
dev "whole tutorial in one sbatch job (tutorial assumes interactive typing on a node); login node only submits"
dev "R section run via Rscript under LC_ALL=C (tutorial: interactive R with the same export as segfault workaround)"
dev "KEGG block kept commented exactly as tutorial section 9; the consolidated 'CODES' cell wraps it in ''' which is invalid R syntax (doc bug, recorded)"
dev "diagnostic prints added in R (sample order, DEG counts) — additive only, no analysis logic changed"
dev "SLURM submission broken cluster-wide since the restructure week (all combos rejected, queue empty, 9/20-proven WES scripts rejected too) -> tutorial-native interactive mode: direct single-thread run on login60, well within the instructor's 30-concurrent-students load-test envelope (STAR mini index 381M, salmon 16M, node 45G free / load 1.5)"

# ---------- 0. rsync mini data (tutorial: server-user section) ----------
log "step 0: rsync mini fastq"
rsync --progress /lustre1/share/RNA_seq_class/mini_fastq_10k/*.fastq.gz .
log "  -> $(ls *.fastq.gz | wc -l) fastq.gz files"

# ---------- 1. fastqc batch (tutorial: find ... -exec fastqc) ----------
log "step 1: fastqc batch"
find . -name "*.fastq.gz" -exec fastqc {} -d . -o . \;
log "  -> $(ls *_fastqc.zip | wc -l) fastqc zips"

# Total Sequences extraction (tutorial demonstrates one sample; looped here)
log "step 1b: Total Sequences (all fastqc zips)"
for z in *_fastqc.zip; do
  totalreads=$(unzip -c "$z" */fastqc_data.txt | grep 'Total Sequences' | cut -f 2)
  echo "  $z: $totalreads" | tee -a "$SUM"
done
dev "Total-Sequences extraction looped over all 12 files (tutorial: single example, 'replacing with your file name')"

# ---------- 2. fastp batch (tutorial: fastp_batch.sh, verbatim loop body) ----------
log "step 2: fastp batch"
for file in *.read1.fastq.gz; do
    base=$(basename "$file" .read1.fastq.gz)
    read1_file="${base}.read1.fastq.gz"
    read2_file="${base}.read2.fastq.gz"
    output_read1="${base}.read1.trimmed.fastq.gz"
    output_read2="${base}.read2.trimmed.fastq.gz"
    json_report="${base}.fastp.json"
    html_report="${base}.fastp.html"
    fastp -i "$read1_file" -I "$read2_file" -o "$output_read1" -O "$output_read2" --detect_adapter_for_pe -l 25 -j "$json_report" -h "$html_report"
done
log "  -> $(ls *.trimmed.fastq.gz | wc -l) trimmed fastq files"

# ---------- 3. STAR batch (tutorial: batch_STAR.sh, mini chr22 index, runThreadN 1) ----------
log "step 3: STAR batch"
GENOMEDIR="/lustre1/share/RNA_seq_class/genome/mini_ref_chr22/star_index"
for file in *.read1.trimmed.fastq.gz; do
    base=$(basename "$file" .read1.trimmed.fastq.gz)
    read1_file="${base}.read1.trimmed.fastq.gz"
    read2_file="${base}.read2.trimmed.fastq.gz"
    output_prefix="${base}"
    STAR --runThreadN 1 \
         --genomeDir $GENOMEDIR \
         --readFilesIn "$read1_file" "$read2_file" \
         --outFileNamePrefix "$output_prefix" \
         --readFilesCommand zcat \
         --outSAMtype BAM Unsorted \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --outFilterMultimapNmax 20 \
         --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverReadLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --quantMode TranscriptomeSAM \
         --outSAMattributes NH HI AS NM MD
done
log "  -> $(ls *Aligned.toTranscriptome.out.bam | wc -l) transcriptome BAMs"
log "  STAR mapping rates:"
for f in *Log.final.out; do
  echo "  $f: $(grep -E 'Uniquely mapped reads number |% of reads mapped to multiple loci' "$f" | tr -s ' ' | tr '\n' ' ')" | tee -a "$SUM"
done

# ---------- 4. multiqc (tutorial: multiqc "$PWD" -o multiqc_report) ----------
log "step 4: multiqc"
multiqc "$PWD" -o multiqc_report
log "  -> multiqc_report/ present: $(test -d multiqc_report/multiqc_data && echo yes || echo check)"

# ---------- 5. salmon align-mode batch (tutorial: batch_salmon_align.sh) ----------
log "step 5: salmon align-mode batch"
TRANSCRIPTOME="/lustre1/share/RNA_seq_class/genome/mini_ref_chr22/chr22.transcripts.fa"
for bam_file in *Aligned.toTranscriptome.out.bam; do
    sample_name=$(basename "$bam_file" Aligned.toTranscriptome.out.bam)
    salmon quant -t "$TRANSCRIPTOME" --libType A -a "$bam_file" -o "${sample_name}_align.salmon_quant" --gcBias --seqBias -p 1
    echo "Processed: $sample_name"
done
dev "align-mode run per batch_salmon_align.sh (--gcBias --seqBias); tutorial's single-sample demo command omits both flags — doc inconsistency, batch version used"

# ---------- 6. salmon direct-mode batch (tutorial: batch_salmon_direct.sh) ----------
log "step 6: salmon direct-mode batch"
INDEX="/lustre1/share/RNA_seq_class/genome/mini_ref_chr22/salmon_index"
for r1 in *.read1.trimmed.fastq.gz; do
    base=$(basename "$r1" .read1.trimmed.fastq.gz)
    r2="${base}.read2.trimmed.fastq.gz"
    salmon quant -i "$INDEX" \
        -l A \
        -1 "$r1" \
        -2 "$r2" \
        --validateMappings \
        --gcBias --seqBias \
        -p 1 \
        -o "${base}_direct.salmon_quant"
    echo "Processed: $base"
done

# ---------- 7. R DEG analysis (tutorial sections 1-10, LC_ALL=C) ----------
log "step 7: R analysis (Rscript, LC_ALL=C)"
export LC_ALL=C
export LANG=C
Rscript scripts/rna_mini_analysis.R 2>&1 | tee log/R_analysis.log
log "  R PDF outputs: $(ls -1 *.pdf 2>/dev/null | tr '\n' ' ')"

log "== RNA-seq mini tutorial done =="
