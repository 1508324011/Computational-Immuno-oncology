# RNA-seq 全流程教程 · 矫正增强版（Corrected & Annotated Edition）

: Document version 1.0 · 2026-10-09 · reproduced from `RNA_seq_tutorial_mini.html` and `RNA_seq_homework.html` (course originals), with a two-round full-pipeline reproduction (mini 2026-09-21, homework/full 2026-10-09) and every deviation verified, recorded, and corrected.

## 0. 关于本版（About this edition）

This edition is the RNA-seq counterpart of the WES corrected roadmap (`WES/ref/WES_roadmap_corrected.md`). It mirrors the course originals command-for-command, then adds what the originals do not carry: per-step rationale, real measured numbers from two complete runs, a pitfall registry with root causes, and corrections that were actually executed and verified.

**Sources (course originals, kept untouched):**

| File | Role |
| --- | --- |
| `RNA_seq_tutorial_mini.html` | 9-step mini pipeline (10k reads/sample), the in-class executable path |
| `RNA_seq_homework.html` | Same pipeline at full read depth on the full GRCh38 index |
| `rna_seq.pdf` | 14-page theory slides (RNA-seq principles, STAR/Salmon positioning, R basics) |
| `RNA_seq_code_explaination.pdf` | 9-page line-by-line command commentary |

**The two-round reproduction** (both rounds: commands verbatim, minimal wrapping, every deviation logged):

| Round | Date | Data | STAR index | Where it ran | Runtime |
| --- | --- | --- | --- | --- | --- |
| mini | 2026-09-21 | 12 fastq × 10,000 reads (20 MB) | chr22 mini (381 MB) | login node, serial | ~4 min |
| homework (full) | 2026-10-09 | 12 fastq, 118k–227k read pairs (112 MB) | full GRCh38 (~31 GB RAM) | login node, serial | 13 min |

**Where the rounds ran, and why.** The tutorial's native mode is interactive execution on the login node — the two HTML documents contain no `sbatch`/`srun`/SLURM references at all; the only cluster-touching command is the initial `ssh`. The teacher's own load test (`02_simulate_30_students_pipeline.sh`, artifacts under `/lustre1/share/RNA_seq_class/loadtest/mini30_20260928_*`) simulated 30 students running this pipeline **concurrently on the login node**, which is direct evidence that the course design expects login-node execution. The mini round ran during a campus-wide SLURM submission outage (every `sbatch` rejected with `Invalid account or account/partition combination`; `sacct` empty since 9/15; the WES scripts that had run successfully on 9/20 were also rejected — cluster-side reconfiguration, recorded in the round's deviations log). The homework round ran on the login node by explicit directive, although SLURM had recovered by probe time (test-only allocation succeeded on cn-long). Both choices are recorded in each round's `log/deviations.log`.

**What "corrected" means here** (same discipline as the WES edition): run the printed commands as written, log every pitfall encountered or predicted, then fix and verify the fixes with real outputs — the fix evidence is in this repository, not just described. Numbers in this document are measured values from the two runs, not estimates.

**How to read each step:** every pipeline step carries the same block set — *Original* (the printed command, verbatim), *Rationale* (what the tool and its key parameters actually do), *What we ran* (measured numbers from both rounds), *Pitfalls* (verified, with evidence), *Interpretation*, *Tuning*, and *Extensions* (where to go next when this step is routine).

## 1. 问题清单速览（Issues at a glance）

Verified pitfalls in the course originals. Full evidence and corrections are in the step sections; the two most consequential (RNA-08/10) have dedicated sections.

| ID | Where | Symptom / root cause | Correction |
| --- | --- | --- | --- |
| RNA-01 | R sample_info | Sample identity is **positional end-to-end**: `condition` assumes `list.files()` returns HBR×3 then UHR×3 (alphabetical luck), and unnamed files give `colnames = NULL`, so edgeR names samples Sample1–6 | Name samples explicitly; assert order before building the design matrix |
| RNA-02 | mini salmon single-sample | Demo cell omits `--gcBias --seqBias` but the batch script includes them → systematic correction difference between demo and batch outputs | Use the batch parameters everywhere |
| RNA-03 | homework salmon single-sample | Align-mode demo cell has a mangled filename (UHR sample name + a pasted HBR `_ERCC-Mix2...chr22` suffix) → file not found | Batch script is correct; use it |
| RNA-04 | homework batch_salmon_quasi | `for r1 in *read1.fastq.gz` matches only **raw** reads — quasi mode quantifies **untrimmed** data (mini direct mode used trimmed) | Run as written, record; or point the glob at `*.read1.trimmed.fastq.gz` |
| RNA-05 | homework batch_STAR.sh | Shebang printed as `!/bin/bash` (missing `#`) | Cosmetic; run via `bash script.sh` |
| RNA-06 | homework multiqc cell | Placeholder path `/lustre1/share/RNA_seq_student/yourfold` | Point at the real working directory |
| RNA-07 | consolidated "CODES" cell (both docs) | KEGG block wrapped in `'''` — **not an R comment**; parses as one giant string literal (silently no-ops) and any apostrophe inside would close the string early, turning the rest into live code | The step-by-step cells correctly use `#`; treat the consolidated cell as non-executable |
| RNA-08 | homework consolidated "CODES" cell | Still filters `"_direct"` (mini naming) — homework outputs are `_quasi` → zero files found, downstream failure | Use the step-by-step path (`grepl("quasi")`) |
| RNA-09 | homework heatmap block | `gene_map` keyed by **versioned** `transcript_id`, looked up with **version-stripped** rownames → every label misses → all 50 rows render as "NA", silently | Strip versions from the map keys as well (one line); verified 0/50 NA |
| RNA-10 | teacher reference outputs | Teacher's `heatmap_DEG.pdf` row labels are **versioned ENST IDs** — neither gene names nor the all-NA the printed code produces: the reference run used a third code variant | Compare against reference with care; see Appendix B |
| RNA-11 | teacher `class.R` | `salmon_files <- salmon_files[-2]` drops the 2nd sample — cleanup of teacher's own duplicate directory, not a method | Do not replicate; be aware when diffing against reference |
| RNA-12 | interpretation pitfall | Quasi/direct mode maps only 30–37% of reads (selective alignment is conservative); align mode reports 100% (input BAM is pre-aligned) — the two are not comparable denominators | Compare rates only within a mode |
| RNA-13 | both docs | `--runThreadN 1` in single-sample demos vs `4` in the batch script — 1 is a 30-student classroom load decision, not a technical limit | Solo runs may raise threads (mind the shared login node) |
| RNA-14 | R startup | `address (nil) ... Segmentation fault` in R on this cluster unless `LC_ALL=C`/`LANG=C` are exported first | Export both before `R`/`Rscript` (tutorial notes this; keep it) |
| RNA-15 | volcano plot | `ggrepel: N unlabeled data points (too many overlaps)` — informational (9 in mini, 280 in full), labels are just decluttered | Benign; raise `max.overlaps` if more labels are wanted |

## 2. 环境与数据（Environment and data）

**Original (setup cells, both docs):**

```bash
ssh .....
mkdir /lustre1/share/RNA_seq_student_2026_homework/{studentId}_{name}
cd /lustre1/share/RNA_seq_student_2026_homework/{studentId}_{name}
conda activate RNA-seq
```

**Rationale.** The shared environment (`/lustre1/share/miniconda3/envs/RNA-seq`) carries every tool the pipeline needs, pinned by the course: fastp 0.23.4, FastQC 0.12.1, STAR 2.7.11b, salmon 1.10.3, multiqc, samtools, gffread, and R 4.4.0 with edgeR/limma/tximport/ggplot2/ggrepel/pheatmap/clusterProfiler/org.Hs.eg.db/biomaRt. Activating it (instead of installing per-student) is what makes a 30-student classroom feasible and what the teacher's load test validates.

**The dataset is SEQC — the only RNA-seq set with a truth label.** Six paired-end libraries: **HBR ×3** (Human Brain Reference) and **UHR ×3** (Universal Human Reference RNA), each spiked with ERCC controls at two mixes (HBR = Mix2, UHR = Mix1). "Differential expression" here has a known answer: brain-enriched genes are high in HBR, and UHR — a pooled-tissue mix that includes immune-cell RNA — carries high immunoglobulin transcripts. Every downstream sanity check in this edition leans on that truth (section 10). The `chr22`/`ErccTranscripts` in filenames marks that the shipped FASTQs were pre-filtered to chr22 (+ ERCC) content.

**Teacher-built assets (do not rebuild — the originals mark index-building cells `Do NOT RUN`):**

| Asset | Path | Size | Used by |
| --- | --- | --- | --- |
| chr22 mini STAR index | `/lustre1/share/RNA_seq_class/genome/mini_ref_chr22/star_index` | 381 MB | mini round STAR |
| chr22 mini salmon index | `.../mini_ref_chr22/salmon_index` | small | mini round direct mode |
| chr22 transcriptome fasta | `.../mini_ref_chr22/chr22.transcripts.fa` | small | mini round align mode |
| full GRCh38 STAR index | `/lustre1/share/RNA_seq/genome` | SA 24.9 GB + Genome 3.2 GB + SAindex 1.6 GB (≈31 GB RAM to load) | homework round STAR |
| full transcriptome fasta | `/lustre1/share/RNA_seq/genome/GRCh38_no_alt_analysis_set_gencode.v36.transcripts.fa` | 368 MB | homework align mode |
| full salmon index (k31) | `/lustre1/share/RNA_seq/genome/salmon_index` | 801 MB | homework quasi mode |
| gene annotation RDS | `/lustre1/share/RNA_seq_class/genome/gene_anno_gencode.v36.rds` | — | R annotation (transcript_id → gene_name) |

**Data scale (measured, both rounds):**

| Round | Files | Reads per sample (pairs) | fastp kept |
| --- | --- | --- | --- |
| mini | 12 × 10,000 reads | 10,000 reads (5,000 pairs) | — (trimmed lightly) |
| homework | 12 fastq.gz (112 MB) | HBR 118,571–144,826; UHR 162,373–227,392 | 97.4–98.7% |

## 3. 数据获取与原始质控（rsync + fastqc）

**Original (homework; mini is identical in shape):**

```bash
rsync --progress /lustre1/share/RNA_seq/RNA_seq_files/*.fastq.gz .
find . -name "*.fastq.gz" -exec fastqc {} -d . -o . \;
totalreads=$(unzip -c HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2_fastqc.zip */fastqc_data.txt | grep 'Total Sequences' | cut -f 2)
echo "Total Sequences: $totalreads"
```

**Rationale.** `rsync` (not `cp`) gives progress and resumability over Lustre. FastQC per file yields the standard module battery (per-base quality, adapter traces, duplication, GC); the homework note `DO NOT UNZIP` matters twice over — the reads must stay gzipped for `--readFilesCommand zcat` later, and `unzip -c` reads `fastqc_data.txt` **out of the zip without extracting anything**.

**What we ran.** Both rounds: 12 fastqc zips. `Total Sequences` per file confirmed the scales above (mini: 10,000 per file; homework: read1 = read2 counts, 118,571…227,392 per file).

**Pitfalls.** None at this step beyond the obvious (running fastqc from the wrong directory silently finds nothing; the `find -exec` runs serially — fine at this size).

**Tuning.** FastQC accepts `-t` threads; at 12 files it is unnecessary. **Extensions.** multiqc (section 6) aggregates these reports; `fastq_screen` would add contamination profiling.

## 4. 接头与质量修剪（fastp）

**Original (batch cell; single-sample demo cell is the same command one-off):**

```bash
for file in *.read1.fastq.gz
do
    base=$(basename "$file" .read1.fastq.gz)
    read1_file="${base}.read1.fastq.gz"
    read2_file="${base}.read2.fastq.gz"
    output_read1="${base}.read1.trimmed.fastq.gz"
    output_read2="${base}.read2.trimmed.fastq.gz"
    json_report="${base}.fastp.json"
    html_report="${base}.fastp.html"
    fastp -i "$read1_file" -I "$read2_file" -o "$output_read1" -O "$output_read2" --detect_adapter_for_pe -l 25 -j "$json_report" -h "$html_report"
done
```

**Rationale.** `--detect_adapter_for_pe` infers the adapter pair from PE read overlap — no adapter sequence needs to be known ahead. `-l 25` discards reads shorter than 25 bp after trimming (STAR degrades on very short inputs). The per-sample JSON/HTML pair is what multiqc will consume.

**What we ran (homework round, all six samples):**

| Sample | Reads in | Reads out | Kept |
| --- | --- | --- | --- |
| HBR_Rep1 | 237,142 | 234,120 | 98.7% |
| HBR_Rep2 | 289,652 | 285,044 | 98.4% |
| HBR_Rep3 | 259,572 | 256,166 | 98.7% |
| UHR_Rep1 | 454,784 | 444,834 | 97.8% |
| UHR_Rep2 | 324,746 | 317,040 | 97.6% |
| UHR_Rep3 | 370,884 | 361,336 | 97.4% |

**Interpretation.** ~97–99% survival means the libraries are clean: adapter content low, quality high — consistent with a commercial reference RNA lot. UHR trims slightly more than HBR (pool composition, slightly more short fragments).

**Pitfalls.** None new; note that the trimmed filenames (`*.read1.trimmed.fastq.gz`) are what STAR's loop matches, while the homework quasi-mode loop matches the **raw** names (RNA-04) — the two globs encode two different preprocessing intents.

**Tuning.** `--qualified_quality_phasing`/`--cut_front` etc. exist but SEQC libraries need none of them. **Extensions.** cutadapt/trim-galore are equivalent alternatives; `--overlap_len` tightening helps when PE overlap detection over-trims.

## 5. 剪接比对（STAR）

**Original (homework batch script; single-sample demo identical except `--runThreadN 1`):**

```bash
GENOMEDIR="/lustre1/share/RNA_seq/genome"
for file in *.read1.trimmed.fastq.gz
do
    base=$(basename "$file" .read1.trimmed.fastq.gz)
    read1_file="${base}.read1.trimmed.fastq.gz"
    read2_file="${base}.read2.trimmed.fastq.gz"
    output_prefix="${base}"
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
```

**Rationale (parameters that matter).** STAR is splice-aware: it seeds with exact k-mers (default 50 for 100 bp reads), extends, and scores spliced alignments against a splice-junction database built into the index. `--outFilterType BySJout` keeps an intron-containing alignment only if its junction is supported by at least one *other* read (or the annotated SJ database) — the main junk-intron filter. `--alignSJoverhangMin 8` requires ≥8 bases on each side of a *novel* junction. `--outFilterMultimapNmax 20` allows up to 20 loci (immunoglobulin and repetitive regions genuinely multi-map). `--outFilterMismatchNmax 999` + `--outFilterMismatchNoverReadLmax 0.04` together cap mismatches at 4% of read length. Intron bounds 20…1,000,000 and mates-gap 1 Mb are the GENCODE-compatible envelope. `--quantMode TranscriptomeSAM` is **the step's load-bearing choice**: it additionally emits `<prefix>Aligned.toTranscriptome.out.bam`, coordinates converted to transcript space — the exact input salmon's align mode consumes (section 7). `--outSAMattributes NH HI AS NM MD` emits multimapping hierarchy and mismatch tags that downstream tools expect.

**What we ran (unique / multi / unmapped-too-short, per sample):**

| Round | Sample | Unique % | Multi-loci % | Too short % |
| --- | --- | --- | --- | --- |
| mini | HBR_Rep1 | 53.25 | 0.68 | 12.48 |
| mini | HBR_Rep2 | 53.03 | 0.68 | 11.98 |
| mini | HBR_Rep3 | 53.33 | 0.48 | 12.90 |
| mini | UHR_Rep1 | 46.76 | 2.02 | 11.86 |
| mini | UHR_Rep2 | 49.94 | 2.08 | 11.49 |
| mini | UHR_Rep3 | 45.06 | 1.90 | 12.87 |
| homework | HBR_Rep1 | 51.82 | 1.09 | 47.08 |
| homework | HBR_Rep2 | 51.70 | 1.11 | 47.19 |
| homework | HBR_Rep3 | 51.77 | 1.06 | 47.17 |
| homework | UHR_Rep1 | 45.30 | 2.77 | 51.92 |
| homework | UHR_Rep2 | 48.47 | 3.02 | 48.50 |
| homework | UHR_Rep3 | 44.56 | 2.77 | 52.66 |

**Interpretation.** Three signals, all stable across rounds and depths: (1) **HBR maps ~5 points better than UHR** — the brain reference is a single-tissue RNA with fewer immunoglobulin/repetitive transcripts; UHR's pool includes IG genes that genuinely multi-map (UHR multi-loci 1.9–3.0% vs HBR 0.5–1.1%). (2) **Unique rates are identical between the 10k-read mini and the 100k+ read full round** (HBR 53.0–53.3 vs 51.7–51.8; UHR 45.1–49.9 vs 44.6–48.5) — the profile is a property of the reads, not of depth. (3) The homework round's ~47–53% "unmapped: too short" (vs 11–13% in mini, same reads, different index) is a property of the **full-genome index + BySJout combination** on this pre-filtered library content: reads whose best full-genome alignment falls below the filter thresholds are dropped as "too short". The mini chr22 index reproduces the course's own numbers; the full index trades coverage for stringency. Record, don't panic: the uniquely-mapped fraction that feeds quantification is consistent.

**Pitfalls.** RNA-05 (shebang `!/bin/bash` in the printed batch script — cosmetic). RNA-13 (thread count 1 vs 4 — classroom load design; our homework round used 4 with a memory guard: before every STAR invocation, wait until `MemAvailable` ≥ 37 GB, because loading the full index costs ~31 GB and the login node is shared). The mini index needs only ~0.4 GB — that round needed no guard.

**Tuning.** `--runThreadN` up on a quiet node; `--twopassMode Basic` improves novel-junction sensitivity at ~2× cost; `--outSAMtype BAM SortedByCoordinate` saves a samtools sort if genome-coordinate tools come next.

**Extensions.** samtools `sort`/`index` + `flagstat` for BAM QC; RSeQC (`infer_experiment`, `inner_distance`) for strandedness and insert profiles; `--quantMode GeneCounts` for a quick raw count matrix alongside the transcriptome BAM.

## 6. 汇总质控（multiqc）

**Original (homework cell; mini uses the real working directory):**

```bash
multiqc /lustre1/share/RNA_seq_student/yourfold
```

**Rationale.** multiqc walks a directory tree, recognizes tool outputs it knows (fastqc zips, fastp JSON/HTML, STAR `Log.final.out`, salmon `quant.sf`, …), and merges them into one report with per-sample columns — the fastest way to see all twelve fastq files and six samples side by side.

**What we ran.** `multiqc "$PWD" -o multiqc_report` (RNA-06: the printed path is a placeholder — `yourfold` is literally "your folder"). Both rounds produced the aggregated report; the STAR and fastp numbers quoted in this document are cross-checkable there.

**Pitfalls.** RNA-06 (placeholder path). Also note: multiqc only reports what it finds under the given path — pointing it at the wrong directory "succeeds" with an empty report rather than failing.

**Tuning.** `-f` to rebuild, `--cl` to embed the command line in the report, a `multiqc_config.yaml` to pin module order for teaching. **Extensions.** `--pdf` export; running multiqc at the end and at the mid-pipeline point (post-STAR) gives before/after views.

## 7. 定量 I：比对模式（salmon align mode）

**Original (homework batch cell — the correct one; see RNA-03 for the single-sample cell's mangled filename):**

```bash
TRANSCRIPTOME="/lustre1/share/RNA_seq/genome/GRCh38_no_alt_analysis_set_gencode.v36.transcripts.fa"
for bam_file in *Aligned.toTranscriptome.out.bam; do
    sample_name=$(basename "$bam_file" Aligned.toTranscriptome.out.bam)
    salmon quant -t "$TRANSCRIPTOME" --libType A -a "$bam_file" -o "${sample_name}_align.salmon_quant" --gcBias --seqBias
    echo "Processed: $sample_name"
done
```

**Rationale.** This is *alignment-based* quantification: STAR already decided where each read goes; salmon's align mode (`-a`) takes the transcriptome-coordinate BAM, assigns multi-mapping reads fractionally by EM, and estimates abundance per transcript. `--gcBias --seqBias` correct observed sequence-content and GC biases from the alignments themselves. `-t` is the transcriptome fasta the BAM's transcript names refer to; `--libType A` lets salmon infer the library type.

**What we ran.** All six samples in both rounds; **100.0% of fragments "mapped" in every sample** — which is exactly what it sounds like and not a quality number: the input BAM contains only fragments STAR already aligned, so the denominator is "aligned fragments", not "reads in the library" (RNA-12).

**Pitfalls.** RNA-03 (the single-sample demo cell's filename is a UHR/HBR paste-mix — copy it and it fails at file-not-found). RNA-02 (the mini single-sample cell omits `--gcBias --seqBias` while its batch cell includes them — quantifying one sample by the demo command and the rest by the batch script produces systematically different corrections across samples).

**Interpretation.** Align mode inherits every STAR decision (junction filters, multimap caps, the too-short drops). Its counts are "STAR's view of the transcriptome"; quasi mode (section 8) is salmon's own, more conservative view — that contrast is the point of running both.

**Tuning.** `-p` threads; `--posBias` if positional attenuation is suspected. **Extensions.** RSEM is the classical alternative for the same BAM; `--numBootstraps` adds uncertainty replicates tximport can consume.

## 8. 定量 II：选择性比对模式（salmon quasi / `--validateMappings`）

**Original (homework batch cell):**

```bash
INDEX="/lustre1/share/RNA_seq/genome/salmon_index"
for r1 in *read1.fastq.gz; do
    base=$(basename "$r1" .read1.fastq.gz)
    r2="${base}.read2.fastq.gz"
    salmon quant -i $INDEX \
        -l A \
        -1 "$r1" \
        -2 "$r2" \
        --validateMappings \
        --gcBias --seqBias \
        -p 8 \
        -o "${base}_quasi.salmon_quant"
    echo "Processed: $base"
done
```

**Rationale.** Quasi mode maps raw reads directly against a transcriptome index (no genome, no junctions to discover). `--validateMappings` (salmon ≥1.0's default mapping engine) validates candidate mappings by selective alignment — an alignment-score threshold decides membership instead of a plain k-mer hash hit, which is why its mapping rate reads *lower* than STAR's: reads without a good transcriptome alignment are discarded rather than force-placed.

**What we ran (mapping rates, `num_mapped / num_processed` from `aux_info/meta_info.json`):**

| Round | Mode | HBR_Rep1 | HBR_Rep2 | HBR_Rep3 | UHR_Rep1 | UHR_Rep2 | UHR_Rep3 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| mini | direct (trimmed reads) | 36.5% | 36.5% | 36.4% | 33.8% | 31.2% | 32.2% |
| homework | quasi (raw reads, as written) | 35.2% | 35.3% | 35.1% | 32.6% | 30.3% | 32.0% |

**Pitfalls.** RNA-04 is the important one: `for r1 in *read1.fastq.gz` matches only the **raw** files (`*.read1.trimmed.fastq.gz` ends differently), so the homework quasi mode quantifies **untrimmed** reads while the mini direct mode quantifies trimmed — we ran it as written and recorded it (adapter fragments are short enough that mapping rates barely move, but the preprocessing story differs between the two documents). Also RNA-12: 30–37% here versus 100% in align mode is a **denominator difference**, not a quality collapse.

**Interpretation.** Two independent effects set the ~1/3 rate: the reads are pre-filtered to chr22+ERCC content while the index is the full 231,448-transcript GENCODE v36 (salmon discards what doesn't align well anywhere), and selective alignment is score-thresholded. The rates are remarkably stable across rounds (35.1–36.5% HBR, 30.3–33.8% UHR) — again a property of the libraries, reproducible at two depths.

**Tuning.** `-p 8` is printed in the homework batch cell — on a shared login node that is already a real thread budget. **Extensions.** `--numBootstraps 100` + tximport's `countsFromAbundance` path gives bootstrap-aware variance; `--seqBias` interacts with fragment-length distributions — check `fldLen.png` in the output when QC-ing.

## 9. 读入定量结果并过滤（tximport → edgeR）

**Original (R cells, step-by-step path):**

```r
salmon_files <- list.files(pattern = "quant\\.sf$", full.names = TRUE, recursive = TRUE)
salmon_files <- salmon_files[grepl("quasi", salmon_files)]     # homework; mini greps "_direct"
txi <- tximport(salmon_files, type = "salmon", txOut = TRUE)
dge <- DGEList(counts = txi$counts)
counts_threshold <- 10
samples_threshold <- 2
keep <- rowSums(dge$counts >= counts_threshold) >= samples_threshold
dge_filtered <- dge[keep, ]
sample_info <- data.frame(
  sample = colnames(dge_filtered),
  condition = factor(c(rep("HBR", 3), rep("UHR", 3)))  # change based on samples
)
design <- model.matrix(~condition, data = sample_info)
```

**Rationale.** `tximport` gathers per-sample `quant.sf` (Name, Length, TPM, NumReads) into one matrix; `txOut = TRUE` keeps the rows at **transcript level** (231,448 rows against the full index). The filter (`≥10 counts in ≥2 samples`) kills rows too sparse to model — at this library size that is almost everything, by construction of the dataset (chr22+ERCC-origin reads against a whole-transcriptome index).

**What we ran (and what to check before the design matrix — RNA-01):**

| Check | mini | homework |
| --- | --- | --- |
| transcripts in index | 4,989 (chr22 mini) | 231,448 (full GENCODE v36) |
| transcripts kept | 111 | 831 |
| `colnames(txi$counts)` | `NULL` → Sample1–6 | `NULL` → Sample1–6 |
| order that saved the day | HBR×3 then UHR×3 (alphabetical) | HBR×3 then UHR×3 (alphabetical) |

**Pitfalls.** RNA-01 is the structural one, in three links: (1) `list.files()` returns alphabetical order and the printed `condition` vector hard-codes HBR-first — correct **only by luck of the filenames**; (2) because the files are passed unnamed, `colnames` is `NULL` and edgeR invents Sample1–6 — sample identity is purely positional from here on; (3) nothing downstream (design matrix, contrasts, `sample_info` rownames in the heatmap) can catch a permutation. The fix that costs one line: `names(salmon_files) <- sub("...sample pattern...", ..., salmon_files)` and build `condition` from those names. Both rounds' checks printed above — order held both times, recorded, not assumed. RNA-14: `export LC_ALL=C; export LANG=C` before R (the tutorial carries this itself).

**Interpretation.** 111/4,989 and 831/231,448 look like different survival rates but are the same phenomenon at two index breadths: ~2% of rows pass in both. Depth sets *how many* transcripts clear 10 counts (111 → 831), the index sets the denominator.

**Tuning.** Raising `counts_threshold` with depth (10 is for this teaching scale; whole-transcriptome studies commonly require ≥10–20 CPM-like floors). **Extensions.** `txOut = FALSE` + a `tx2gene` frame aggregates to gene level before edgeR (the more common production path); DESeq2 with `tximport` is the drop-in alternative; `rtracklayer`/`biomaRt` (already loaded by the tutorial) can rebuild the annotation.

## 10. 差异表达、注释与火山图（glmLRT + volcano）

**Original (R cells):**

```r
dge_filtered <- estimateDisp(dge_filtered, design)
fit <- glmFit(dge_filtered, design)
results <- glmLRT(fit)
results_table <- topTags(results, n = Inf)$table
results_table$regulate <- ifelse(results_table$FDR < 0.05 & results_table$logFC > 1, "up-regulated",
                                   ifelse(results_table$FDR < 0.05 & results_table$logFC < -1, "down-regulated", "not sig"))
```

```r
gene_anno <- readRDS("/lustre1/share/RNA_seq_class/genome/gene_anno_gencode.v36.rds")
rownames_nover <- sub("\\..*$", "", rownames(results_table))
anno_transcript_nover <- sub("\\..*$", "", gene_anno$transcript_id)
match_index <- match(rownames_nover, anno_transcript_nover)
results_table$gene_name <- gene_anno$gene_name[match_index]
results_table$label <- ifelse(results_table$regulate != "not sig", results_table$gene_name, NA_character_)
```

**Rationale.** edgeR's generalized-linear-model path: `estimateDisp` estimates the negative-binomial dispersion, `glmFit` fits, `glmLRT` likelihood-ratio-tests the last term (here `conditionUHR`). The `regulate` rule is a conjunction: FDR < 0.05 **and** |logFC| > 1 — statistical significance and effect size together, which is why the volcano's dashed lines sit at ±1 and −log10(0.05). The annotation block strips GENCODE version suffixes (`ENST00000615943.1` → `ENST00000615943`) **on both sides** before `match()` — note this is exactly the discipline the heatmap block forgets (section 12).

**What we ran:**

| | mini | homework (full) |
| --- | --- | --- |
| up-regulated | 15 | 171 |
| down-regulated | 30 | 148 |
| not significant | 66 | 512 |
| Matched (gene_name) | 111 / 111 | 831 / 831 |
| top up (by FDR) | IGLC2 (logFC ≈ +10.3, FDR ≈ 1e-145), IGLC3, XBP1, MYH9 | IGLC2, IGLC3, … |
| top down | MIAT (≈ −8.5), YWHAH, MAPK8IP2, ACO2, SYNGR1, SEPTIN3, CBX7 | SYNGR1, SEPTIN3, YWHAH, RPL3, CBX6, … |

**Interpretation — the SEQC truth check.** Every direction matches the known answer: **UHR up**: immunoglobulin lambda chains IGLC2/IGLC3 (UHR's pool includes immune-cell RNA) — textbook; **HBR up (i.e., UHR-down)**: MIAT, YWHAH, MAPK8IP2, ACO2, SYNGR1, SEPTIN3 — brain-enriched transcripts. The top gene is the same transcript at both depths (IGLC2, |logFC| > 10 at FDR ~1e-145) — the signal was already saturated at 10k reads; depth bought *more* called genes (45 → 319), not a different answer. That is the core lesson of the mini/homework pairing.

**Pitfalls.** RNA-15 (`ggrepel: 9 / 280 unlabeled data points` — decluttering, benign). glmLRT takes the fitted object; passing extra arguments changes the tested term (caught and corrected during our mini-round harvest — the tutorial's own call `glmLRT(fit)` is the right one).

**Tuning.** `topTags(n = Inf)` then sorting by logFC shows strong-effect genes the FDR sort hides. **Extensions.** `decideTests` + Venn across contrasts; `topTags` with `adjust.method = "BH"` is the default FDR here.

## 11. 富集分析（GO — and the KEGG that never runs）

**Original (step-by-step cells — the executable path):**

```r
gene_list <- results_table$gene_name[results_table$regulate == "up-regulated"]
gene_list <- gene_list[!(is.na(gene_list) | gene_list == "")]
#KEGG analysis
#entrez_ids <- mapIds(...)
#kegg_enrich <- enrichKEGG(...)
#go analysis
go_enrich <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                      ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1)
go_dotplot <- dotplot(go_enrich, showCategory = 10, title = "GO Enrichment")
ggsave("GO_dotplot_up_regulate.pdf", plot = go_dotplot, width = 7, height = 5)
```

**Rationale.** Enrichment asks whether the up-regulated gene set is annotation-fold-enriched versus the background. `keyType = "SYMBOL"` feeds gene symbols directly (no manual ID conversion); `ont = "ALL"` pools BP+MF+CC; `pvalueCutoff = 1` disables the display-side cutoff so the top-10 sort is by actual significance.

**What we ran.** 15 up genes (mini) / 171 (homework). Homework top GO categories (p.adjust 0.0075–0.0088): *viral RNA genome replication*, *single-stranded viral RNA replication via double-stranded DNA intermediate*, *glycosyl compound catabolic process*, *protein targeting to mitochondrion*, *DNA dealkylation*, *cytidine catabolic process*, *cytidine deamination*, *cytidine-to-uridine editing*. The same top categories appear in the teacher's reference `GO_dotplot.pdf` (*viral RNA genome replication*, *glycosyl compound catabolic process*, p.adjust 0.0075) — cross-validated against the reference, not just self-consistent.

**Pitfalls — RNA-07/08, both about the consolidated "CODES" cell.** The step-by-step cells comment KEGG out with `#` (correct). The consolidated cell instead wraps the whole KEGG block in `''' … '''`. In R, `'''` is **not** a comment: it parses as an empty string `''` followed by a string-opener `'`, so the entire block becomes **one giant string literal** that R evaluates and auto-prints — it *accidentally* behaves like a block comment, doing nothing. The failure mode is latent: any apostrophe inside the block would close the string early and the remaining code becomes **live**. Treat `'''` blocks as non-executable prose, always. And RNA-08: the homework's consolidated cell still filters `"_direct"` (the mini naming) — with homework outputs named `_quasi` it finds zero files and dies at tximport; the step-by-step path (`grepl("quasi")`) is the executable one.

**Interpretation.** The cytidine-deamination/cytidine-to-uridine-editing cluster is not noise: immunoglobulin-expressing cells run AID/APOBEC-driven cytidine deamination for affinity maturation — the GO readout of the same IG signal the volcano shows. The "viral RNA genome replication" style categories are an annotation-side effect of the IG/ERV-region genes in the up set (read with the teacher's reference, which shows the same).

**Tuning.** `enrichKEGG` needs ENTREZ IDs (the commented block does the conversion) and network access to KEGG REST — the reason it's the step usually left off in class. **Extensions.** `compareCluster` across up/down/both sets; GSEA (`gseGO`) on the whole ranked table uses more than the top set.

## 12. 热图——沉默的全 NA 坑（the silent all-NA bug）

**Original (homework heatmap block, verbatim):**

```r
logCPM_scaled <- cpm(dge_filtered, log = TRUE, prior.count = 1)
rownames(logCPM_scaled) <- sub("\\..*$", "", rownames(logCPM_scaled))
degs <- results_table[results_table$FDR < 0.05 & abs(results_table$logFC) > 1, ]
top50 <- head(rownames(degs[order(degs$FDR), ]), 50)
top50_clean <- sub("\\..*$", "", top50)
logCPM_scaled_sub <- logCPM_scaled[top50_clean, ]
gene_map <- setNames(gene_anno$gene_name, gene_anno$transcript_id)
rownames(logCPM_scaled_sub) <- gene_map[rownames(logCPM_scaled_sub)]
rownames(sample_info) <- colnames(logCPM_scaled)
pdf("heatmap_DEG.pdf", width = 10, height = 8)
pheatmap(logCPM_scaled_sub, cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, annotation_col = sample_info,
         main = "Heatmap of Differentially Expressed Genes")
dev.off()
```

**The bug chain, link by link.** Line 2 strips versions from the rownames (`ENST00000390323.5` → `ENST00000390323`); line 6 builds `gene_map` with **versioned** keys (`setNames(gene_anno$gene_name, gene_anno$transcript_id)` — no stripping); line 7 looks the stripped names up in the versioned map → **every lookup misses**; assigning an all-`NA` vector to rownames is legal in R; pheatmap renders it — **fifty rows labeled "NA", no error, no warning**. `pdftotext heatmap_DEG.pdf` from our real run shows exactly 50 × `NA`. This is the same version-mismatch class as the annotation block in section 10, which gets it right by stripping both sides before `match()` — and the mini version's heatmap block gets it right too (a `clean_id()` helper applied to *both* the data rows and the annotation IDs, plus `make.unique()` and `labels_row =`).

**The one-line correction (executed and verified):**

```r
gene_map <- setNames(gene_anno$gene_name, sub("\\..*$", "", gene_anno$transcript_id))
```

Our corrected run (`heatmap_DEG_corrected.pdf`): **0 of 50 NA**, row labels IGLC2, IGLC3, MCM5, SERPIND1, LIF, … — the real top-DEG genes.

**Pitfalls.** RNA-09 (above). RNA-10: the teacher's reference `heatmap_DEG.pdf` has **versioned ENST IDs as row labels** — neither gene names (mini-style) nor all-NA (printed-homework-style): the reference run used yet another code variant, so it is not a valid oracle for this block (Appendix B). The homework block also drops the mini version's `scale = "row"` — rows are raw logCPM, so a highly expressed gene visually swamps the pattern scale="row" would show.

**Interpretation.** The heatmap's job here is the visual assertion of the DEG story: samples cluster by condition (HBR vs UHR columns separate), top genes split into the two direction blocks. Labels are the only thing the bug broke — the numbers, clustering, and annotation bar were always right, which is precisely why it survived into the printed handout.

**Tuning.** `scale = "row"` (z-scores per row), `cutree_rows` to force block structure, `gaps_col = 3` between the two conditions. **Extensions.** `ComplexHeatmap` for annotation-heavy figures; `label_maxsize`/fontsize when 50 labels crowd.

## 13. 双轮对比（mini vs homework/full, side by side)

| | mini (2026-09-21) | homework/full (2026-10-09) |
| --- | --- | --- |
| Reads per sample | 10,000 | 118,571–227,392 pairs |
| STAR index | chr22 mini, 381 MB | full GRCh38, ~31 GB RAM |
| STAR unique (HBR / UHR) | 53.0–53.3% / 45.1–49.9% | 51.7–51.8% / 44.6–48.5% |
| STAR too short | 11.5–12.9% | 47.1–52.7% (full index + BySJout stringency) |
| salmon align mode | 100% (aligned-fragment denominator) | 100% (same) |
| salmon direct/quasi | 31.2–36.5% (trimmed) | 30.3–35.3% (raw, as written) |
| transcripts in index | 4,989 | 231,448 |
| transcripts kept (≥10 counts, ≥2 samples) | 111 (2.2%) | 831 (0.36%) |
| DEG (FDR<0.05 and \|logFC\|>1) | 45 (15 up / 30 down) | 319 (171 up / 148 down) |
| top gene both directions | IGLC2 (+10.3) / MIAT (−8.5) | IGLC2 / SYNGR1-cluster |
| Matched annotation | 111/111 | 831/831 |
| GO top terms | (small set) | viral RNA genome replication etc., match teacher reference |
| heatmap labels | gene names (correct block) | all "NA" as printed → corrected 0/50 NA |
| wall time | ~4 min | 13 min |
| where | login node, serial | login node, serial (memory-guarded STAR) |

**What scales and what doesn't.** Wall time grew 4 → 13 min (index load dominates; the full index is ~80× the mini one). DEG count grew 45 → 319 (depth buys calls). The *answers* didn't move: mapping-rate profiles, top genes, and directions are depth-invariant — which is the whole point of running the mini first. The one thing that flipped between rounds was an artifact of the *documents*, not the data: the mini heatmap block annotates correctly, the homework block doesn't.

## 14. 检查清单、产物清单与总结（Checklist, outputs, summary）

**检查清单 (checklist — with both rounds' measured answers):**

| # | Check | mini | full |
| --- | --- | --- | --- |
| ck1 | sample order before design matrix (HBR×3 first?) | held | held |
| ck2 | `colnames(txi$counts)` NULL → Sample1–6 (positional identity) | yes | yes |
| ck3 | STAR unique rate in 44–54% band, HBR above UHR | 45.1–53.3% | 44.6–51.8% |
| ck4 | multi-loci low, UHR > HBR (IG content) | 0.5–2.1% | 1.1–3.0% |
| ck5 | quasi mapping 30–37%, align 100% (different denominators — don't compare) | 31.2–36.5% | 30.3–35.3% |
| ck6 | kept transcripts ~2% and ~0.4% of index at the two breadths | 111/4,989 | 831/231,448 |
| ck7 | DEG direction matches SEQC truth (IGL up in UHR, brain genes down) | yes | yes |
| ck8 | top gene identical across depths | IGLC2 | IGLC2 |
| ck9 | `Matched: n/n` printed by the annotation block | 111/111 | 831/831 |
| ck10 | GO top terms vs teacher reference | — | match (viral RNA genome replication, p.adjust ~0.008) |
| ck11 | heatmap labels are gene names, not "NA" | yes (correct block) | 0/50 NA **after correction** (50/50 NA as printed) |
| ck12 | `LC_ALL=C`/`LANG=C` exported before R | yes | yes |
| ck13 | ggrepel "unlabeled" warnings noted benign | 9 | 280 |

**产物清单 (outputs):** both rounds' artifacts live in the repository under `RNA-seq/results/` (analysis products only — no raw reads, no BAMs): fastqc/fastp reports, STAR `Log.final.out` ×6 each, multiqc reports, salmon `quant.sf` for all 24 quantifications, the R analysis log, the four homework-round PDFs (volcano / GO / **heatmap as printed, all-NA** / **heatmap corrected, 0-NA**) plus the three mini-round PDFs, and both rounds' full `summary` + `deviations` logs and runner scripts. The heavy intermediates (BAMs, trimmed/ raw FASTQs, full salmon directories) stay in the cluster run directories (`/lustre1/user/bjx131_pkuhpc/rna_run_mini_20260921/`, `/lustre1/user/bjx131_pkuhpc/rna_run_homework_20261009/`).

**总结.** The pipeline is sound; the pitfalls are in the hand-written seams around it — sample identity carried only by position and alphabetical luck, two salmon modes whose mapping rates share no denominator, a quasi-mode loop that quietly quantifies raw reads, and a heatmap annotation rewrite that silently renders 50 unlabeled rows. Every one of them was reproduced on real runs, and every fix was executed and verified. The SEQC dataset did its job throughout: the biological answer (IGLC2 up, brain genes down) was stable at both depths and matched the teacher's reference wherever the reference was comparable — which, as Appendix B shows, is itself something to verify, not assume.

## 附录 A. 偏差日志（both rounds, condensed but verbatim-faithful）

Full logs: `RNA-seq/results/mini/log/` and `RNA-seq/results/homework/log/` (summary + deviations + R analysis).

**mini round (2026-09-21) — selection:**

```text
DEVIATION: workdir = own run dir (tutorial: classroom share dir); no writes to the classroom share
DEVIATION: tutorial-native interactive mode: serial single-process run on the login node —
  campus-wide SLURM submission outage (all sbatch rejected; sacct empty since 9/15; the WES
  scripts that ran successfully on 9/20 also rejected → cluster-side, not ours)
DEVIATION: SLURM env-var guard + deviation logging added to the teacher's load-test script shape
DEVIATION: single consolidated runner script (tutorial assumes interactive typing)
FINDING: sample order fragility holds (list.files alphabetical, HBR×3 first, condition correct)
FINDING: tximport unnamed files → colnames NULL → edgeR Sample1-6 (positional dependence)
FINDING: KEGG block is a doc bug (step-by-step cells comment it with #; consolidated wraps in ''')
FINDING: ggrepel 9 unlabeled benign
```

**homework round (2026-10-09) — selection:**

```text
DEVIATION: workdir = own run dir (tutorial: /lustre1/share/RNA_seq_student_2026_homework/{id})
DEVIATION: tutorial-native interactive mode on the login node, per user directive 2026-10-09 —
  SLURM had RECOVERED by probe time (test-only alloc OK on cn-long) but the tutorial's own
  design is direct interactive execution; recorded, not used this round
DEVIATION: restart guards (skip steps whose outputs exist); STAR memory guard (wait for
  MemAvailable ≥ 37 GB before each invocation — index ~31 GB on a shared login node)
DEVIATION: STAR --runThreadN 4 per the homework batch script (single-sample demo cell says 1)
DEVIATION: multiqc: doc cell has placeholder path; used $PWD like the mini version
DEVIATION: salmon quasi mode runs on RAW reads (glob *read1.fastq.gz matches only untrimmed) —
  inconsistent with mini's direct mode on trimmed; run AS WRITTEN, recorded
DEVIATION: R step-by-step path used (grepl("quasi")); consolidated CODES cell still filters
  "_direct" (mini naming) AND wraps KEGG in ''' (invalid R) — unusable for homework outputs
FINDING: heatmap row labels ALL "NA" (gene_map keyed by versioned transcript_id, looked up
  with version-stripped rownames) — silent; correction pass heatmap_DEG_corrected.pdf: 0/50 NA
FINDING: quasi mapping rate 30–35% on raw reads (mini direct ~32%, consistent)
FINDING: teacher reference heatmap rows are versioned ENST IDs — teacher's actual run used a
  third code variant; third teacher-vs-doc deviation (cf. class.R salmon_files[-2], the
  mangled single-sample salmon filename)
RESULT: GO dotplot top categories match teacher reference (viral RNA genome replication,
  glycosyl compound catabolic process; p.adjust 0.0075 vs 0.0082–0.0088) — cross-validated
```

## 附录 B. 教师参考产物 vs 印刷版文档（evidence, all paths verbatim)

The teacher's own run left artifacts in `/lustre1/share/RNA_seq/`. Three independent deviations from the printed documents:

1. **The reference heatmap is a third variant.** `/lustre1/share/RNA_seq/heatmap_DEG.pdf` row labels are versioned ENST IDs (e.g. `ENST00000407418.8`). The printed homework block produces all-"NA" labels (RNA-09, reproduced); the mini block produces gene names. The reference matches neither — its code kept the versioned ENST rownames and never renamed them.
2. **`class.R` drops a sample.** `salmon_files <- salmon_files[-2]` — the 2nd quantification directory is removed before tximport. Consistent with cleanup of a duplicated directory in the teacher's own run area (the share contains both `UHR_Rep3...chr22Aligned...bam` and a doubled-name `UHR_Rep3...chr22_ERCC-Mix2...chr22Aligned...bam`), not with any method in the printed documents.
3. **The mangled single-sample salmon cell was actually run.** The doubled-name BAMs above are exactly what the homework's single-sample align cell (RNA-03, UHR name + pasted HBR `_ERCC-Mix2...chr22` suffix) produces — the reference area carries the bug's own output filename.

Where the reference *is* comparable, it agrees with our runs: `GO_dotplot.pdf` top categories and p.adjust match ours (Appendix A, RESULT). The WES edition's lesson repeats verbatim: reference artifacts are evidence of *someone's actual commands*, not a proof of the printed ones — diff against them only where the code path is known.

## 附录 C. 复现（reproduction)

Both rounds' runner scripts and the R scripts are stored verbatim in `RNA-seq/results/mini/scripts/` and `RNA-seq/results/homework/scripts/`, next to their full logs. The homework runner documents its own wrapping deviations in its header comment (SLURM header kept for provenance only — the round ran tutorial-native on the login node). Re-running the mini round takes ~4 minutes and needs only the shared environment; re-running the homework round takes ~13 minutes serial and needs ~31 GB free for each STAR invocation.
