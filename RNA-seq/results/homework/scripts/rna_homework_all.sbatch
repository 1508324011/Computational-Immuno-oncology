#!/bin/bash
#SBATCH --job-name=rna_hw
#SBATCH --partition=cn-long
#SBATCH --qos=bjx131cnl
#SBATCH --account=bjx131_g1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --output=/lustre1/user/bjx131_pkuhpc/rna_run_homework_20261009/log/rna_hw_%j.out
#SBATCH --error=/lustre1/user/bjx131_pkuhpc/rna_run_homework_20261009/log/rna_hw_%j.err

# =============================================================================
# RNA-seq HOMEWORK (full-size reads, full GRCh38 STAR index) — tutorial-native
# serial run on the LOGIN NODE, per user directive 2026-10-09.
# Tutorial source: /gpfs1/home/bjx131_pkuhpc/class/RNA-seq/RNA_seq_homework.html
# Commands verbatim + minimal wrapping (SLURM header kept for provenance only;
# this round is NOT submitted via sbatch). Same discipline as WES/mini rounds:
# failures/deviations recorded in log/deviations.log.
# =============================================================================
RUN=/lustre1/user/bjx131_pkuhpc/rna_run_homework_20261009
cd "$RUN"
mkdir -p log

# --- conda activate FIRST: `set -u` breaks conda's activate scripts (#25) ---
source /lustre1/share/miniconda3/etc/profile.d/conda.sh
conda activate RNA-seq

SLURM_JOB_ID="${SLURM_JOB_ID:-login60-direct}"
set -euo pipefail
SUM=log/summary_homework.log
DEV=log/deviations.log
: > "$SUM"; : > "$DEV"
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$SUM"; }
dev(){ echo "[$(date '+%F %T')] DEVIATION: $*" | tee -a "$DEV" "$SUM"; }
wait_mem(){ # restart-safety guard: full GRCh38 STAR index needs ~31G RAM
  local need_mb=$1 avail
  while :; do
    avail=$(awk '/MemAvailable/ {print int($2/1024)}' /proc/meminfo)
    if [ "$avail" -ge "$need_mb" ]; then break; fi
    log "  waiting for memory: ${avail}G avail < ${need_mb/1024}G needed; sleep 120"
    sleep 120
  done
}

log "== RNA-seq HOMEWORK start ($SLURM_JOB_ID, $(hostname)) =="
log "tools: $(fastp --version 2>&1 | head -1) | STAR $(STAR --version | tail -1) | salmon $(salmon --version) | fastqc $(fastqc --version 2>&1 | tail -1)"

dev "workdir = own run dir (tutorial: /lustre1/share/RNA_seq_student_2026_homework/{id}); no writes to classroom share"
dev "tutorial-native interactive mode: serial single-process run on login node, per user directive 2026-10-09 — SLURM had RECOVERED by probe time (test-only alloc OK on cn-long) but the tutorial's own design is direct interactive execution; recorded, not used this round"
dev "single consolidated runner script (tutorial assumes interactive typing)"
dev "restart guards added: steps whose outputs already exist are skipped (additive wrapper for a multi-hour serial run)"
dev "STAR memory guard added before each invocation: wait until MemAvailable >= 37G (index SA 25G + Genome 3.2G + SAindex 1.6G ~ 31G RAM on a shared login node)"
dev "STAR --runThreadN 4 per the homework batch_STAR.sh (single-sample demo cell says 1; batch cell is the homework path)"
dev "batch_STAR.sh in the doc lost the '#' of its shebang (!/bin/bash) — cosmetic doc bug, recorded; loop body used verbatim"
dev "multiqc: doc cell has placeholder path /lustre1/share/RNA_seq_student/yourfold; used \$PWD like the mini version"
dev "Total-Sequences extraction looped over all 12 files (tutorial: single example)"
dev "single-sample salmon align cell has a mangled filename (HBR/UHR mix in name); the batch_salmon_align.sh version is correct and used"
dev "salmon quasi mode runs on RAW reads (glob '*read1.fastq.gz' matches only untrimmed files) — inconsistent with mini version's direct mode on trimmed reads; run AS WRITTEN, recorded"
dev "R section: step-by-step path used (salmon_files filtered by grepl('quasi')); the consolidated CODES cell still filters '_direct' (mini naming) AND wraps KEGG in ''' (invalid R) — unusable for homework outputs, doc bug recorded"
dev "KEGG block kept commented exactly as written"
dev "diagnostic prints added in R (sample order, kept transcripts, DEG counts) — additive only"
dev "teacher reference class.R contains salmon_files[-2] (drops 2nd sample) — teacher-side quirk, NOT replicated here; noted for the corrected edition"

# ---------- 0. rsync full-size data ----------
log "step 0: rsync homework fastq (112M source)"
rsync --progress /lustre1/share/RNA_seq/RNA_seq_files/*.fastq.gz .
log "  -> $(ls *.fastq.gz | wc -l) fastq.gz files"

# ---------- 1. fastqc batch ----------
log "step 1: fastqc batch (serial as written)"
find . -name "*.fastq.gz" -exec fastqc {} -d . -o . \;
log "  -> $(ls *_fastqc.zip | wc -l) fastqc zips"

log "step 1b: Total Sequences (all fastqc zips)"
for z in *_fastqc.zip; do
  totalreads=$(unzip -c "$z" */fastqc_data.txt | grep 'Total Sequences' | cut -f 2)
  echo "  $z: $totalreads" | tee -a "$SUM"
done

# ---------- 2. fastp batch (verbatim) ----------
log "step 2: fastp batch"
for file in *.read1.fastq.gz; do
    base=$(basename "$file" .read1.fastq.gz)
    [ -f "${base}.fastp.json" ] && { log "  skip $base (done)"; continue; }
    read1_file="${base}.read1.fastq.gz"
    read2_file="${base}.read2.fastq.gz"
    output_read1="${base}.read1.trimmed.fastq.gz"
    output_read2="${base}.read2.trimmed.fastq.gz"
    json_report="${base}.fastp.json"
    html_report="${base}.fastp.html"
    fastp -i "$read1_file" -I "$read2_file" -o "$output_read1" -O "$output_read2" --detect_adapter_for_pe -l 25 -j "$json_report" -h "$html_report"
done
log "  -> $(ls *.trimmed.fastq.gz | wc -l) trimmed fastq files"

# ---------- 3. STAR batch (verbatim loop body, runThreadN 4, memory guard) ----------
log "step 3: STAR batch (full GRCh38 index ~31G, serial, guarded)"
GENOMEDIR="/lustre1/share/RNA_seq/genome"
for file in *.read1.trimmed.fastq.gz; do
    base=$(basename "$file" .read1.trimmed.fastq.gz)
    [ -f "${base}Log.final.out" ] && { log "  skip $base (done)"; continue; }
    read1_file="${base}.read1.trimmed.fastq.gz"
    read2_file="${base}.read2.trimmed.fastq.gz"
    output_prefix="${base}"
    wait_mem 37888
    log "  STAR $base (MemAvailable $(awk '/MemAvailable/ {print int($2/1048576)}' /proc/meminfo)G)"
    STAR --runThreadN 4 \
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

# ---------- 4. multiqc ----------
log "step 4: multiqc"
multiqc "$PWD" -o multiqc_report
log "  -> multiqc_report/ present: $(test -d multiqc_report/multiqc_data && echo yes || echo check)"

# ---------- 5. salmon align-mode batch (verbatim) ----------
log "step 5: salmon align-mode batch"
TRANSCRIPTOME="/lustre1/share/RNA_seq/genome/GRCh38_no_alt_analysis_set_gencode.v36.transcripts.fa"
for bam_file in *Aligned.toTranscriptome.out.bam; do
    sample_name=$(basename "$bam_file" Aligned.toTranscriptome.out.bam)
    [ -f "${sample_name}_align.salmon_quant/quant.sf" ] && { log "  skip $sample_name (done)"; continue; }
    salmon quant -t "$TRANSCRIPTOME" --libType A -a "$bam_file" -o "${sample_name}_align.salmon_quant" --gcBias --seqBias
    echo "Processed: $sample_name" | tee -a "$SUM"
done

# ---------- 6. salmon quasi-mode batch on RAW reads (verbatim, as written) ----------
log "step 6: salmon quasi-mode batch (raw reads, -p 8, as written)"
INDEX="/lustre1/share/RNA_seq/genome/salmon_index"
for r1 in *read1.fastq.gz; do
    base=$(basename "$r1" .read1.fastq.gz)
    [ -f "${base}_quasi.salmon_quant/quant.sf" ] && { log "  skip $base (done)"; continue; }
    r2="${base}.read2.fastq.gz"
    salmon quant -i $INDEX \
        -l A \
        -1 "$r1" \
        -2 "$r2" \
        --validateMappings \
        --gcBias --seqBias \
        -p 8 \
        -o "${base}_quasi.salmon_quant"
    echo "Processed: $base" | tee -a "$SUM"
done

# ---------- 7. R DEG analysis (step-by-step blocks 21-31, LC_ALL=C) ----------
log "step 7: R analysis (Rscript, LC_ALL=C, quasi pattern)"
export LC_ALL=C
export LANG=C
Rscript scripts/rna_homework_analysis.R 2>&1 | tee log/R_analysis.log
log "  R PDF outputs: $(ls -1 *.pdf 2>/dev/null | tr '\n' ' ')"

log "== RNA-seq HOMEWORK done =="
