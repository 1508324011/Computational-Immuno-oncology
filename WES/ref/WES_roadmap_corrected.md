---
title: "A Practical Roadmap for Whole-Exome Sequencing Data Analysis — Corrected & Annotated Edition"
subtitle: "WES 数据分析实战路线图 · 矫正与增强版"
lang: en
toc: true
toc-depth: 3
number-sections: false
---

# About this edition · 关于本版

This document is a **corrected and heavily annotated edition** of the course roadmap
*"A Practical Roadmap for Whole-Exome Sequencing Data"* (the original HTML/PDF remain in
this directory for cross-reference). It was produced by **actually executing the original
roadmap, command by command**, on a compute cluster — twice:

| Round | Dataset | Depth over targets | Runtime | Purpose |
| --- | --- | --- | --- | --- |
| **Demo** | `/lustre1/share/data_new` (18.2k / 25.8k read pairs, OC / PBMC) | ~0.5× | ~25 min | Exercise the full 13-step pipeline mechanically |
| **Full** | `/lustre1/share/data/OC_WES_all` (37.7M / 52.5M read pairs) | **measured 54.6× / 91.3×** usable on-target | ~4 h (8 SLURM jobs, dependency-linked) | Produce real results, compare with the teacher's `analysis_full` |

> **中文批注：** 本版路线图的做法是——先把原版命令**原样**在集群上跑两轮（一轮 0.5× 降采样教学数据、一轮全量数据），把过程中遇到/发现的**每一个坑都记录在案**，再回头写这份矫正版。所以每个步骤里你都会看到"**原版怎么写 → 实际发生了什么 → 应该怎么改 → 为什么**"的完整链条。两轮实跑本身就是最好的对照组：0.5× 数据会暴露哪些步骤在"低深度下失效"，全量数据会暴露哪些步骤在"命令写法上就有问题"。

## What each step contains · 每步包含什么

For every pipeline step, this edition adds six learning blocks around the original commands:

1. **Principle 原理** — what the algorithm actually does, in plain language.
2. **Commands** — the original roadmap commands (kept verbatim where they are correct; **corrections are shown explicitly and always justified**).
3. **What happened 实跑结果** — measured outputs from both runs (demo vs full).
4. **Pitfalls & fixes 坑与矫正** — every problem we hit, with evidence.
5. **Result interpretation 结果解读** — how to read the numbers; what "good" looks like.
6. **Tuning & beyond 调优与扩展** — parameter knobs, and what else the tool can do for you beyond this roadmap.

## Issue registry at a glance · 问题清单速览

The full detailed registry is in the **Appendix A** section at the end of this document. Summary:

| # | Severity | Where | Issue | Status in our runs |
| --- | --- | --- | --- | --- |
| 1 | **Major** | §7.3 | `MergeVcfs` of two full ApplyVQSR outputs **duplicates every variant** (each output contains *all* variants; merging them doubles the callset) | Confirmed in demo (14,782 → 29,564 records) **and in the teacher's own `analysis_full` output (204,178 → 408,356 records)** |
| 2 | **Major** | §8.3 | Rich ANNOVAR protocol references databases that do not exist in the installed humandb (`dbnsfp30a`, `exac03`, `gnomad_genome`, `avsnp147` vs installed `dbnsfp54a`, `gnomad211_exome`, `gnomad41_exome`, `avsnp151`) | Confirmed failure in both runs |
| 3 | **Major** | §7 | VQSR on a 2-sample, low-depth callset is statistically unreliable even when it runs to completion | Ran "successfully" on demo, but see the tranches analysis in §7 |
| 4 | **Major** | §10 | `cnvkit.py batch -m amplicon` on a **hybrid-capture** panel (SureSelect V6) — wrong mode; also no gene names in BED, no access file | Ran, but produces genome-wide nonsense segments on demo (cn up to 6345) |
| 5 | Medium | §9 | Mutect2 lacks a **Panel of Normals** (one exists on the cluster: `1000g_pon.hg38.vcf.gz`); germline resource uses 1000G high-conf instead of gnomAD AF-only | Runs, reduced filtering power |
| 6 | Medium | §0.2/§12 | Data path confusion: README points to `OC_WES` (deprecated downsample), new roadmap points to `data_new`; old roadmap pointed to an **empty** directory | Three inconsistent references; this edition fixes them |
| 7 | Medium | §3 | Literal roadmap writes a 33 GB SAM to disk before sorting — works, but wastes ~33 GB/sample and an extra hour | We piped `bwa | sort` in the full round (documented deviation) |
| 8 | Minor | §9 | TMB counts PASS variants from `somatic.filtered.vcf.gz` while §9.2 separately writes `somatic.pass.vcf.gz` (never used) | Cosmetic inconsistency |
| 9 | Minor | §4/§11 | No on-target depth metric anywhere (flagstat is genome-wide) — the checklist asks "is coverage sufficient?" but no command measures it | This edition adds `CollectHsMetrics` (§11) — measured: OC 54.6×, PBMC 91.3×, demo 0.3×/0.5× |
| 10 | Minor | naming | `PBMC` vs `blood` naming mixed (RG `SM:blood`, dirs `PBMC`); teacher's own output dir is spelled `2_vairant_call` | Consistent-but-confusing; noted |
| 11 | Minor | §5 | `BaseRecalibrator` runs whole-genome (no `-L`) — not wrong, ~2× slower for WES | Kept verbatim in both rounds |
| 12 | **Major**\* | reference | **The teacher's reference result itself deviates from the printed roadmap:** single-sample "joint" call (blood gVCF's `SM` tag = `OC`, so CombineGVCFs collapsed the pair) and **no `-L`** (72% of its 204,178 sites are off-target), plus issue #1's duplication (408,356 = 2×204,178). \*Not a roadmap-text bug — an execution deviation in the reference artifacts | Read from the VCF headers; documented in Appendix A12 so callset comparisons are interpreted correctly |

> **中文批注：** 最重要的是 #1（VQSR 合并导致每条变异重复两遍，连老师的参考结果都带着这个 bug）和 #4（CNVkit 模式选错）。这两个都不影响"管线能不能跑完"，但**直接影响结果文件能不能用**——这正是"跑通了"和"跑对了"的区别。

# 0. Environment, directories and data · 环境与数据

## 0.1 Required tools

Everything runs from the shared conda environment on the cluster:

```bash
source /lustre1/share/miniconda3/etc/profile.d/conda.sh
conda activate wes
```

Tool versions verified on a compute node (not the login node — see §0.4):
fastp 0.23.2, BWA 0.7.19-r1273, samtools, GATK 4.6.2.0 (Java), CNVkit 0.9.11, ANNOVAR (2020s build with hg38 humandb), bcftools.

> **中文批注：** 集群上所有重计算都必须 `sbatch` 到 `cn-long` 分区做，登录节点只用来编辑文件和提交作业。本文所有实跑结果都来自 SLURM 作业（见 §13 复现说明）。

## 0.2 Directory convention

The roadmap's environment block, kept verbatim (this is the *new* version's only
difference from the old one, which pointed `$DATA` at an **empty** directory — see issue #6):

```bash
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
export DATA=$SHARE/data_new        # demo round; full round uses $SHARE/data/OC_WES_all
export BASE=$HOME/wes_run_cluster  # we used explicit run dirs, see §13
export OUT=$BASE
```

**What lives where (verified by `ls`, not assumed):**

| Path | Content | Notes |
| --- | --- | --- |
| `/lustre1/share/references/` | `hg38.fa` + BWA index + `.dict`/`.fai`; GATK bundle (`dbsnp_138`, Mills, 1000G high-conf, hapmap, omni) **with `.idx` present**; target BED `S07604514_AllTracks_V6_60_hg38.bed`; **`1000g_pon.hg38.vcf.gz` + `.tbi`** (a ready-made Panel of Normals!); `access-5kb.hg38.bed` | read-only share |
| `/lustre1/share/data_new/` | demo FASTQs: `OC/OC_R1.fq.gz` etc. (182k pairs), `PBMC/blood_R1.fq.gz` (258k pairs) | ~0.48% subsample of the full data ≈ 0.5× over targets |
| `/lustre1/share/data/OC_WES_all/` | full FASTQs: OC 37.7M pairs (10.4 Gbp), PBMC 52.5M pairs (15.0 Gbp) | ideal ≈95×/≈135×; **measured usable 54.6× / 91.3×** (§11) |
| `/lustre1/share/data/OC_WES/` | older downsample (same size, **different content** — different md5) | deprecated; README still points here (issue #6) |
| `/lustre1/share/data/analysis_full/` | teacher's reference results (full data): alignment, joint VCF, VQSR, ANNOVAR | **only up to §8**; no somatic/TMB or CNVkit reference; `4_cnv/` is empty |
| `/lustre1/share/annovar_new/humandb/` | hg38 databases: `refGeneWithVer`, `cytoBand`, **`dbnsfp54a`, `gnomad211_exome`, `gnomad41_exome`, `avsnp151`** | the rich-protocol names in §8.3 of the roadmap are **not** here (issue #2) |

**Target region:** the BED contains 243,359 intervals totaling **60.51 Mb** (this number
matters twice: as `-L` intervals and as the TMB denominator in §9.3).

## 0.3 Reference genome

hg38 (`hg38.fa`), pre-indexed for BWA and GATK. **Check before you compute** (we verified
all of these exist — a missing `.idx` is the single most common "first job fails" cause):

```
hg38.fa  hg38.fa.fai  hg38.dict
hg38.fa.{amb,ann,bwt,pac,sa}          # BWA index
dbsnp_138.hg38.vcf{,.idx}             # also needed by HaplotypeCaller -D
Mills_and_1000G_gold_standard.indels.hg38.vcf{,.idx}
1000G_phase1.snps.high_confidence.hg38.vcf{,.idx}
hapmap_3.3.hg38.vcf{,.idx}
1000G_omni2.5.hg38.vcf{,.idx}
```

## 0.4 Roadmap (and how we ran it on SLURM)

The roadmap's steps were executed as sbatch jobs on the `cn-long` partition (20-core
nodes, 25 h wall limit). Three cluster-specific facts that cost us an hour to discover —
record them once, reuse forever:

1. **Account/QoS:** jobs must specify `--account=bjx131_g1 --qos=bjx131cnl`
   (plain `--qos=bjx131cnl` without the account is rejected: *"Invalid qos specification"*).
2. **Memory:** the `cn-long` nodes report `RealMemory=1` to SLURM, so any explicit
   `--mem=NN` is rejected (*"Memory specification can not be satisfied"*). Use
   **`--mem=0`** ("take the whole node's memory").
3. **conda + `set -u`:** if your script uses `set -u`, activate conda **before** turning
   it on — conda's `activate.d` scripts reference unbound variables and will kill the
   job in one second.

All reproduction scripts are in `scripts/` of each run directory (see §13).

**The two rounds at a glance (all jobs on `cn-long`, account `bjx131_g1`, QoS `bjx131cnl`, `--mem=0`):**

| Round | Job(s) | Steps | Wall time | Outcome |
| --- | --- | --- | --- | --- |
| demo | 3048 | 3→10 (verbatim, incl. expected-failure guards) | ~25 min | all steps executed; deviations log: 1 (rich-protocol fail) |
| full | 3049/3050 | prep per sample: fastp → `bwa\|sort` → index/flagstat → MarkDup → BQSR → HC-GVCF | 2 h 20 m / 3 h 42 m (parallel) | both clean |
| full | 3051 | 6.3→8 verbatim (joint → VQSR → MergeVcfs → ANNOVAR basic) | 11 min | duplication bug reproduced (§7.2) |
| full | 3052 | 9 verbatim (Mutect2 → orientation model → FilterMutectCalls → TMB) | 1 h 13 m | TMB = 1.78 mut/Mb |
| full | 3053 | 10 verbatim (`-m amplicon`) | 5.6 min | runs; empty antitarget baseline (§10.2) |
| full | 3054/3056→3061 | **corrected extras**: CollectHsMetrics ×4 samples | 32 min | measured depths: OC 54.6× / PBMC 91.3× (§11) |
| full | 3055 | **corrected extra**: CNVkit `-m hybrid` + access file | 11 min | arm-level events found (§10.3) |
| full | 3056 | **corrected extra**: sequential ApplyVQSR + ANNOVAR rich (installed names) | 9.6 min | zero duplication; rich annotation OK |

*Note how the corrected jobs (3055/3056) are not just "fixed" — they are the second half of the course: the roadmap-as-written is run first (3051/3053), then the corrected command on the same data, so every claim in this edition is a before/after measurement.*

# 1. FASTQ files · 原始数据

## 1.1 Inspect FASTQ

**Principle 原理.** Before any analysis, you must *look* at your data: read counts,
read length, GC profile. Two obvious-but-underestimated checks: (a) FASTQ read length
should match the sequencer mode (this library is 2×137 bp paired-end); (b) GC content
profile per cycle — a spiked or drifting GC curve in *tumor* WES is normal-ish
(capture bias + purity), but a wildly abnormal curve means contamination or adapter read-through.

**Commands** (verbatim from the roadmap): `zcat | head`, `wc -l / 4`, and
`fastp`'s reports later. On the cluster, `zcat file.fq.gz | head -8` for a peek,
`zcat file.fq.gz | awk 'NR%4==2{n++; b=length($0)} END{print n, b}'` for count × length.

**What happened 实跑结果.**

| Sample | Pairs (R1) | Read length | Bases | Depth over 60.51 Mb |
| --- | --- | --- | --- | --- |
| demo OC (`data_new`) | 182,299 | 137 bp | 0.050 Gbp | ideal ~0.82× → **measured on-target 0.3×** (usable) |
| demo PBMC | 258,112 | 137 bp | 0.071 Gbp | ideal ~1.2× → **measured on-target 0.5×** |
| full OC (`OC_WES_all`) | 37,746,701 | 137 bp | 10.41 Gbp | ideal ~172× → **measured on-target 54.6×** |
| full PBMC | 52,508,649 | 137 bp | 14.97 Gbp | ideal ~247× → **measured on-target 91.3×** |

**Interpretation 结果解读.** The "ideal" column is total bases ÷ target size — the
number you get from naive arithmetic. The **measured** column is what CollectHsMetrics
reports as `MEAN_TARGET_COVERAGE` *after* subtracting everything that cannot be used:
off-target bases (20–9%), duplicates (7–11%), overlapping mate bases (**25% — this
library's inserts peak at ~150 bp, so most read pairs overlap each other; see §2.3 and
§11). The 3× gap between ideal and usable depth is the single most common
mis-estimation in WES planning — budget your sequencing accordingly.

> **中文批注：** 这里最值得记住的对比是：全量 OC 理想深度 ~172×（10.41 Gbp ÷ 60.5 Mb），但实测可用靶区深度只有 **54.6×**——差了 3 倍多！差额去向：off-target ~20%、重复 ~7%、双端重叠 ~25%（这个文库插入片段峰值只有 ~150bp，2×137bp 的读段大量互相重叠，重叠部分只算一次）。demo 数据是全量的 0.48% 抽样，实测可用深度只有 0.3×/0.5×（中位数是 **0**）——低于任何可用的变异检测阈值，这就是后面体细胞检测为 0、CNV 全是噪声的根本原因。"理论深度" 和 "实测可用深度" 的差距本身就是 WES 设计中最重要的实务课。

# 2. Read quality control with fastp · 质控与预处理

## 2.1 Create QC output

The roadmap's fastp block (kept verbatim; `$OUT/0_fastq/` used for both raw symlinks and
clean outputs):

```bash
fastp \
  -i "$OUT/0_fastq/OC_R1.fq.gz" \
  -I "$OUT/0_fastq/OC_R2.fq.gz" \
  -o "$OUT/0_fastq/OC_R1.clean.fq.gz" \
  -O "$OUT/0_fastq/OC_R2.clean.fq.gz" \
  --detect_adapter_for_pe \
  --cut_right \
  --cut_right_window_size 4 \
  --cut_right_mean_quality 20 \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --n_base_limit 5 \
  --length_required 50 \
  --thread 8 \
  --html "$OUT/0_fastq/OC/OC.fastp.html" \
  --json "$OUT/0_fastq/OC/OC.fastp.json"
```

## 2.2 Principle 原理

- `--cut_right` is **3′ sliding-window trimming** (default window 4, mean Q20): bases at
  read ends degrade because of phasing/signal decay; trimming them prevents piles of
  false-positive mismatches at read ends later.
- `--detect_adapter_for_pe`: paired-end data lets fastp infer adapter sequences from
  read overlap — more sensitive than a fixed adapter list.
- `--length_required 50`: post-trim length floor. Short fragments contribute
  misalignments and duplicate artifacts.
- PE "overlap analysis" also gives you **insert-size estimation** for free — read it
  from the HTML report; it feeds directly into duplicate interpretation in §4.2.

## 2.3 What happened 实跑结果

Demo round (already run before this edition started):

| Sample | Before → after (reads) | Q30 (before → after) | Insert size (peak) |
| --- | --- | --- | --- |
| demo OC | 377,310 → 364,598 | 95.5% → 97.0% | **150 bp** |
| demo PBMC | 525,176 → 516,224 | 96.9% → 97.8% | **150 bp** |

Full round:

| Sample | Before → after (reads) | Q30 (before → after) | Insert size (peak) |
| --- | --- | --- | --- |
| full OC | 75,813,060 → 72,941,400 | 95.5% → 97.0% | 150 bp |
| full PBMC | 105,258,338 → 103,270,440 | 96.9% → 97.8% | 150 bp |

Two observations worth internalizing: (a) trimming barely removes anything (3.8%/1.9% of
reads) because the raw Q30 is already 95%+ — this is a healthy library; (b) the insert
peak at **150 bp** with 137 bp reads means most read pairs **overlap each other** by
~120 bp. That overlap is later discarded by coverage counters (25% of usable bases —
§11), and it inflates apparent depth if you count reads instead of non-overlapping
bases. Both rounds show the identical profile, as expected for a subsample of the same
library.

## 2.4 Pitfalls & fixes 坑与矫正

- **(p1) Output placement inconsistency:** the roadmap puts clean FASTQs at the *root* of
  `0_fastq/` but reports in `0_fastq/OC/`. Harmless, but a consistent layout
  (`0_fastq/OC/*`) would have been cleaner. We kept it verbatim.
- **(p2) `--cut_right` vs `--cut_front`:** for modern 137 bp reads, 5′ ends are usually
  clean; only 3′ trimming is standard. The roadmap's choice is correct — noted because
  many tutorials add unnecessary `--cut_front`.
- **(p3) fastp's JSON is machine-readable** — the roadmap stops at the HTML; a
  multi-sample **MultiQC** aggregation (see §2.6) is the "beyond" step.

## 2.5 Interpretation 结果解读

"Good" numbers for this library type: ≥90–95% bases ≥Q30 before trimming, a few % of
reads trimmed/adapter-detected, and — most importantly — the **duplication estimate is
NOT in fastp's scope** (it comes from alignment in §4). fastp is about read-level
hygiene; anything requiring mapping positions comes later.

## 2.6 Tuning & beyond 调优与扩展

- `--thread` scales nearly linearly; 8 is fine on a 20-core node alongside other work.
- For UMI libraries use `--umi` instead of relying on coordinate duplicates (§4).
- `--overlap_len_require`/`--overlap_diff_limit` control the PE overlap-based correction
  (`--correction`), which the roadmap does **not** enable — consider adding it for low-depth data.
- Aggregate many samples: `multiqc 0_fastq/` collects fastp JSON/HTML into one report.

> **中文批注：** fastp 三个最常用开关——`--cut_right`（滑窗去低质量尾）、`--detect_adapter_for_pe`（利用 PE 重叠自动识别接头）、`--length_required`（太短的读段宁可不要）。报告里最该看的三处：Q30 比例、接头残留比例（应为 0 或接近 0）、插入片段分布峰（应显著大于读长 137bp，否则有大量重叠读段）。
>
# 3. Alignment with BWA-MEM · 比对

## 3.1 Principle 原理

BWA-MEM seeds exact matches (MEMs) and extends them with local alignment; soft-clipped
ends handle adapter/partial overlaps. Key ideas you should be able to explain:

- **Read groups (RG)** are the provenance header of every read: `ID` (sequencing lane),
  `SM` (**sample** — the identity that GATK uses to group reads; wrong SM = wrong
  genotyping), `LB` (**library** — the unit for duplicate estimation in §4), `PL` platform.
  One BAM = one SM; multiple LBs allowed per SM.
- **`-M`**: marks split hits as *secondary* (instead of *supplementary*) for Picard
  compatibility — cosmetic today, harmless, kept from the roadmap.
- **`-t 8`**: threads. BWA-MEM scales near-linearly up to ~16 threads.

**Commands** (roadmap, verbatim — OC shown; PBMC identical with `SM:blood`, `LB:WES_PBMC`):

```bash
bwa mem -t 8 -M \
  -R "@RG\tID:OC\tPL:illumina\tSM:OC\tLB:WES_OC" \
  "$REFERENCE" \
  "$OUT/0_fastq/OC_R1.clean.fq.gz" \
  "$OUT/0_fastq/OC_R2.clean.fq.gz" \
  > "$OUT/1_alignment/OC/OC.sam"
```

## 3.2 SAM → BAM, sort, index (roadmap 3.2–3.3, verbatim)

```bash
samtools view -bS -@ 8 "$OUT/1_alignment/OC/OC.sam" > "$OUT/1_alignment/OC/OC.bam"
samtools sort -@ 8 -m 4G "$OUT/1_alignment/OC/OC.bam" -o "$OUT/1_alignment/OC/OC_sorted.bam"
samtools index -@ 8 "$OUT/1_alignment/OC/OC_sorted.bam"
```

**Pitfall & fix (issue #7).** The literal chain writes a **33 GB SAM** (full data), then a
BAM copy, then the sorted BAM — three full copies of the same information on disk, and
~⅓ extra wall time. The standard practice (and what we ran for the full round) is a pipe:

```bash
bwa mem -t 8 -M -R "@RG\tID:OC\tPL:illumina\tSM:OC\tLB:WES_OC" \
  "$REFERENCE" R1.clean.fq.gz R2.clean.fq.gz \
  | samtools sort -@ 8 -m 4G -o OC_sorted.bam -
```

The resulting `OC_sorted.bam` is **byte-for-byte equivalent in content** — same records,
same order (chromosomal, then position). Only the intermediate files never exist.
This is the *one deliberate efficiency deviation* we took in the full round; it is
logged in the run's `log/deviations.log`.

> **中文批注：** 为什么要排序？因为下游所有工具（去重、BQSR、HaplotypeCaller、CNVkit）都要求 BAM 按"染色体 → 位置"排好，这样才能顺序流式读取整个外显子组。SAM→BAM 是压缩（文本→BGZF 二进制，约省 4–6 倍空间）；管道化的意义是让"比对输出 → 排序输入"直接在内存/缓冲区里传递，省掉 33 GB 的中间 SAM 文件——结果完全一样，磁盘和时间都省一半以上。

## 3.3 Alignment metrics (roadmap 3.4, verbatim)

```bash
samtools flagstat -@ 8 "$OUT/1_alignment/OC/OC_sorted.bam" > "$OUT/1_alignment/OC/OC.flagstat.txt"
```

## 3.4 What happened 实跑结果

**Demo round (0.5×):**

| Metric | OC | PBMC |
| --- | --- | --- |
| Total reads | 397,780 | 517,164 |
| Primary mapped | 364,582 (100.0%) | 516,221 (100.0%) |
| Properly paired | 83.29% | 99.46% |
| Secondary (split) alignments | 33,182 (9.1% of primary) | 940 (0.18%) |

**Full round:**

| Metric | OC | PBMC |
| --- | --- | --- |
| Total reads | 79,607,191 | 103,455,092 |
| Primary mapped | 72,941,400 (100.00%) | 103,270,440 (100.00%) |
| Unmapped | **311** | 1,594 |
| Properly paired | **83.22%** | 99.44% |
| Mate mapped to a different chr | 11,526,594 (**15.8%**) | 444,126 (0.43%) |
| Secondary (split) alignments | 6,665,791 (9.1%) | 184,652 (0.18%) |

**Interpretation 结果解读.**

- "100% mapped" in **both** rounds — not a subsampling effect: the full round has only
  311 unmapped reads out of 79.6M. This **dataset ships with unmapped reads already
  removed** (a pre-filtered teaching dataset). In a fresh lab BAM, 96–99% mapped is
  typical, and **<90% means contamination or species mismatch**.
- OC's **83.2% properly paired — identical in both rounds — is a real property of the
  tumor library, not an artifact**: 15.8% of OC read pairs have their mate on a
  *different chromosome* (vs 0.4% in PBMC). At this scale that is the alignment
  signature of a heavily rearranged tumor genome (HGSOC typically carries massive
  structural variation). PBMC's 99.4% is what a normal genome looks like.
- 9% secondary alignments in OC vs 0.18% in PBMC: secondary/supplementary reads
  indicate repetitive content (tumor genome + exome bait in repetitive regions). Watch
  it, don't panic over it.

> **中文批注：** flagstat 三问：① 比对率——本数据集 demo 和全量都是 100%，这不是抽样巧合，而是这套教学数据在发布前就把未比对读段剔除过了（全量也只有 311/79.6M 未比对）；新鲜实验室数据应为 96–99%；② properly paired——**OC 的 83.2% 在两轮中完全一致，是肿瘤样本的真实属性而非抽样伤**：15.8% 的 OC 读段配对到不同染色体（PBMC 仅 0.43%），这是重度重排肿瘤基因组（HGSOC 典型特征）在比对层面的签名，也是下游 §10 结构性 CNV 的伏笔；③ secondary/supplementary——肿瘤样本 9% vs 正常 0.18%，重复序列含量差异，外显子组里几个百分点属正常。

## 3.5 Tuning & beyond 调优与扩展

- Threads: `-t 16` halves wall time on an idle 20-core node; `-K` (chunk size) can help huge references.
- `bwa mem -v` for verbose logs; `-Y` (soft-clip supplementary) pairs better with GATK Indel realignment era tools — modern GATK doesn't care.
- **On-target %**: `samtools view -c -L $TARGETS bam` ÷ total — the single most informative QC number for capture-based WES (feeds §11's coverage question).
- For noisy/ancient DNA: `bwa mem -B` (batch size, lower memory) and adjusted scoring.
- **MultiQC** eats flagstats too: `multiqc run_dir/` after all samples finish.

# 4. Mark PCR duplicates · 标记重复

## 4.1 Principle 原理

A "duplicate" = a read pair whose **5′ coordinates and orientation match another pair
exactly**. Random shearing makes this vanishingly unlikely for true molecules; PCR
amplification makes it common. Duplicates are **marked, not removed**: GATK counting
tools ignore them, but they stay in the BAM (auditability, and `samtools view -F 1024`
recovers them).

Why **after sorting**: MarkDuplicates needs coordinate order to see co-located mates
in a streaming fashion. This is also why the **LB tag matters** — duplicates are only
*expected* within a library; two libraries of the same sample are independent
shearatures, and marking across them would be wrong.

**Commands** (roadmap 4.1, verbatim):

```bash
gatk --java-options "-Xmx8g" MarkDuplicates \
  -I "$OUT/1_alignment/OC/OC_sorted.bam" \
  -O "$OUT/1_alignment/OC/OC_sorted.markdup.bam" \
  -M "$OUT/1_alignment/OC/OC_markdup_metrics.txt" \
  --CREATE_INDEX true
```

## 4.2 What happened 实跑结果

**Demo round:** duplication OC **0.033%**, PBMC **0.06%** — essentially zero, exactly
what a 0.48% subsample must produce (the probability that both members of a duplicate
pair survive the subsample is ~0.48%² ≈ 0.002%).

**Full round** (all values from the metrics files' per-library rows):

| Sample | Pairs examined | Dup pairs | Optical dups | PERCENT_DUPLICATION | Est. library size |
| --- | --- | --- | --- | --- | --- |
| OC | 36,468,418 | 2,583,317 | 1,853,690 (71.8% of dups) | **7.08%** | 809,511,121 |
| PBMC | 51,634,767 | 5,562,361 | 3,192,510 (57.4%) | **10.77%** | 478,822,815 |

**Teacher's own full-data run** (from `analysis_full`): OC **7.39%**, PBMC **11.01%**
— our run lands within 0.4 percentage points of both. That agreement is the strongest
single validation that our prep chain (fastp → piped bwa\|sort → MarkDup) reproduces
the reference prep faithfully.

**Interpretation 结果解读.** For hybrid-capture WES at ~100×, expect **5–20%**
duplication. Rising duplication = shrinking library complexity (bad template prep,
over-amplification); the metrics file's `ESTIMATED_LIBRARY_SIZE` is the honest
indicator (hundreds of millions to billions = healthy). One subtlety visible in the
optical-duplicate column: **71.8% of OC duplicates are optical** (adjacent-flowcell
origin), which usually indicates patterned-flowcell clustering rather than PCR
over-amplification — a sequencing-hardware property, not a library-prep failure.

## 4.3 Pitfalls & fixes 坑与矫正

- **(p1) Optical vs PCR duplicates:** the metrics split them. Optical duplicates come
  from adjacent clusters on the flowcell; if optical >> non-optical, suspect the
  patterned flowcell / loading concentration, not library prep.
- **(p2) "Unknown Library"** appears in the teacher's metrics (their RGs lacked `LB`);
  with the roadmap's RGs you get a named row. Keep `LB` — duplicate statistics per library.
- **(p3) UMI libraries** need `--barcoded-*, UMI-aware duplicate marking` (GATK
  `UmiAwareMarkDuplicatesPoN` or fgbio `CallMolecularDuplicates`); coordinate-based
  marking undercounts complexity there.

## 4.4 Tuning & beyond 调优与扩展

- `--TMP_DIR` to a scratch disk for huge BAMs (MarkDuplicates spills heavily).
- `--VALIDATION_STRINGENCY SILENT` if old BAMs trip on NM/MD tag order — better: fix the BAM.
- `samtools markdup` (with `samtools fixmate`) is a faster alternative; GATK's
  `EstimateLibraryComplexity` predicts complexity from read pairs alone.

# 5. Base quality score recalibration (BQSR) · 碱基质量校准

## 5.1 Principle 原理

Phred scores from the instrument are **estimates**, systematically off in ways that
depend on cycle, machine, context. BQSR learns an empirical error model from
"mismatches at sites we KNOW are true" (known-variant resources), then rewrites the
quality scores. The mental model: *if base calls at known-reference positions look
like Q35 in context C, treat future C-context calls as Q35, whatever the instrument said.*

**Commands** (roadmap 5.1–5.2, verbatim — the known-sites trio is dbsnp + Mills + 1000G high-conf, all present with `.idx` on this cluster):

```bash
gatk --java-options "-Xmx20g -XX:ParallelGCThreads=4" BaseRecalibrator \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/OC/OC_sorted.markdup.bam" \
  --known-sites "$DBSNP" \
  --known-sites "$INDELS" \
  --known-sites "$HIGHCONF" \
  -O "$OUT/1_alignment/OC/OC.recal_data.table"

gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" ApplyBQSR \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/OC/OC_sorted.markdup.bam" \
  --bqsr-recal-file "$OUT/1_alignment/OC/OC.recal_data.table" \
  -O "$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam"
```

## 5.2 What happened 实跑结果

Demo: recal table 108 KB / BQSR BAM 54 MB; runtime ~2.5 min per sample.
Full round: recal tables 110 KB per sample (24 covariate rows × ~400 context rows — the
same model shape as demo, just denser counts), BQSR BAMs 8.8 GB (OC) / 8.2 GB (PBMC);
runtimes from the job logs — BaseRecalibrator ~31 min (OC) / ~43 min (PBMC),
ApplyBQSR ~20 min / ~28 min. Note the BQSR BAM is ~1.9× the markdup BAM: ApplyBQSR
rewrites every base quality, and re-compressed BAMs don't always shrink — an easy
place to waste disk if you keep every intermediate.

**Interpretation 结果解读.** The recal table's `EmpiricalQuality` vs `EstimatedQ`
columns are the story: systematic deltas of 1–5 Q are normal and worth correcting;
the BQSR **report plots** (run `gatk AnalyzeCovariates`) make it visual. Two BQSR
truths for this course: (a) with dbsnp as known-sites **plus** Mills + 1000G, novel
real variants are barely affected — the model "believes" the trio; (b) on 0.5× data
BQSR still runs, but the empirical counts per covariate bin are so thin that the
correction is mostly noise — another "works mechanically, meaningless statistically" case.

## 5.3 Pitfalls & fixes 坑与矫正

- **(p1) No `-L` restriction:** BaseRecalibrator scans the *whole genome* even though
  only the exome matters downstream. Adding `-L $TARGETS` cuts wall time ~in half for
  WES. We kept it verbatim (it is not wrong, just slower).
- **(p2) BQSR needs a good chunk of data.** <10⁷ bases on-target makes the model jittery.
- **(p3) Never BQSR without the sequence dictionary match** (`hg38.dict` exists here — fine).

## 5.4 Tuning & beyond 调优与扩展

- `--native-pair-hmm-threads` doesn't apply here (that's HaplotypeCaller); for
  BaseRecalibrator the parallelism is GC threads + (in newer GATK) `--parallelism`.
- Plot covariates: `gatk AnalyzeCovariates -bqsr recal_data.table -O report.pdf` — the
  single best teaching artifact of the whole step.
- Skip-if-known-clean: for deeply-studied platforms some pipelines skip BQSR; GATK
  still recommends it whenever the mismatch model shows bias.

# 6. Germline variant calling · 胚系变异检测

## 6.1 Principle 原理

**HaplotypeCaller (HC)** is not a per-pileup SNP caller: it (1) finds *active regions*
with evidence of variation, (2) **reassembles the haplotype graph** locally, (3) runs a
pair-HMM likelihood of every read against every candidate haplotype, (4) emits the most
likely genotype. This is why it catches indels near short reads far better than pileup callers.

**GVCF mode (`-ERC GVCF`)**: emit *every* block* (variant or reference) with
confidence bands. Cost: bigger files. Benefit: **scalable joint genotyping** — add
sample 31 later by re-running `CombineGVCFs`/`GenotypeGVCFs`, never re-running HC.
This two-sample "GVCF workflow" is the exact pattern to scale to a cohort.

**Commands** (roadmap 6.2–6.4, verbatim):

```bash
gatk --java-options "-Xmx12g -XX:ParallelGCThreads=4" HaplotypeCaller \
  -ERC GVCF -R "$REFERENCE" \
  -I "$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam" \
  -D "$DBSNP" -L "$TARGETS" \
  --native-pair-hmm-threads 8 \
  -O "$OUT/1_alignment/OC/OC.g.vcf"

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" CombineGVCFs \
  -R "$REFERENCE" \
  -V "$OUT/1_alignment/OC/OC.g.vcf" \
  -V "$OUT/1_alignment/PBMC/blood.g.vcf" \
  -O "$OUT/2_variant_call/combined.g.vcf"

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" GenotypeGVCFs \
  -R "$REFERENCE" \
  -V "$OUT/2_variant_call/combined.g.vcf" \
  -D "$DBSNP" \
  -O "$OUT/2_variant_call/OC_blood_variants.vcf"
```

And the optional (roadmap 6.4) VariantAnnotator with OC's BAM — it annotates with
coverage-based fields; on GATK4 it is **not required** before filtering, and the
roadmap itself says so.

## 6.2 What happened 实跑结果

**Demo round:** joint callset **14,782 variant records** (13,632 SNPs + 1,151 indels,
plus mixed) — the tumor/normal pair genotyped jointly; PBMC's germline variants appear
in the tumor columns too.

**Full round:** joint callset **55,123 variant records** (OC+blood, both genotyped,
target-restricted: `-L $TARGETS` in HC/Combine/Genotype).

**Teacher's reference:** **204,178** records. **Beware the naive comparison — these two
numbers are NOT the same quantity** (see issue #12 in Appendix A). Measured facts about
the teacher's VCF: (a) it contains a **single sample** (OC only — their blood gVCF's
RG `SM` tag was set to `OC`, so CombineGVCFs collapsed both into one); (b) it was
called **genome-wide** — HaplotypeCaller/GenotypeGVCFs with **no `-L`** — and **72%
of its sites lie outside the capture targets** (many in low-complexity/off-bait
regions with shallow coverage). Our 55,123 is a 2-sample, on-target callset; their
204,178 is a 1-sample, genome-wide one. Depth is the dominant difference, but scope
and sample count matter too — always read the VCF header (`##GATKCommandLine`,
`#CHROM` line) before comparing callsets.

For a sanity check on our own number: 60.51 Mb of exome at ~1.2–1.4 germline variants
per kb gives ~75–85k per sample; two samples with heavy overlap (both are germline
from related-normal individuals — actually two different individuals, so overlap is
the common-variant fraction ~70%) → ~55–65k unique sites. 55,123 sits right in that
window.

## 6.3 Pitfalls & fixes 坑与矫正

- **(p1) `-D dbsnp` only adds the RS ID annotation** — it does not restrict calling.
  A common misconception is that dbsnp "tells HC what to look for"; it doesn't.
- **(p2) Never mistake this joint VCF for a somatic callset** — these are germline
  variants of BOTH individuals. The roadmap's take-home (and §9) is exactly this point.
- **(p3) PBMC as "normal"** is acceptable for paired somatic calling, but it is a
  different tissue with its own mosaic/ageing mutations; clonal hematopoiesis (CHIP)
  variants in blood WILL appear here — a hot topic in modern immuno-oncology.

## 6.4 Tuning & beyond 调优与扩展

- `--native-pair-hmm-threads 8` gave the single biggest speedup in our runs — use it always.
- For cohorts: `GenomicsDBImport` scales better than CombineGVCFs beyond ~50 samples
  (the roadmap's CombineGVCFs is fine for 2).
- `--interval-padding 50` recovers splice-region variants just outside the BED —
  nearly free, clinically useful.
- QUAL interpretation: it's Phred-scaled genotype confidence, NOT variant frequency.
  A het at QUAL 5,000 means "P(wrong) ≈ 10⁻⁵⁰⁰" — over-precise by construction.

# 7. Variant quality score recalibration (VQSR) · 变异质量再校准 — *the big correction*

## 7.1 Principle 原理

Hard filters use fixed thresholds (QD<2 → bad). **VQSR instead trains a Gaussian
mixture model** on your callset, using "truth" resources (hapmap/omni for SNPs, Mills
for indels) as positive labels and dbsnp as a "known, don't over-trust" population.
Each variant gets a **VQSLOD** score (log-odds of being real vs. artifact given its
annotation vector: QD, MQ, FS, SOR, ReadPosRankSum, MQRankSum). You keep a *truth
sensitivity tranche* (99.5% for SNPs here, 99.0% for indels), filtering the worst-scoring
fraction outside it.

**The non-negotiable precondition: enough training data.** GATK's guidance: ~30 WGS
samples (or a large cohort callset; exomes need more samples than WGS because fewer
variants each). With 2 exomes at 14k variants the model is *mathematically fittable*
but *statistically fictional* — it will happily produce a tranche file.

## 7.2 What happened 实跑结果 — the two headline findings

**Finding A: VQSR "succeeds" on the demo — and that's a trap.** Both VariantRecalibrators
ran to completion on the 14,782-variant demo callset. But look at what it actually learned:

| Evidence | Demo (0.5×, 14.8k) | Our full (55.1k, 2-sample on-target) | Teacher (204k, genome-wide) | Healthy WES expectation |
| --- | --- | --- | --- | --- |
| Accessible truth sites (hapmap+omni ∩ callset) | 7,388 | 28,324 | 85,159 | grows with cohort |
| SNP minVQSLod at 99.9% tranche | −2.56 | −11.62 | −2.47 | mild negative tail |
| SNP minVQSLod at 100% tranche | **−5889.9** | −1111.4 | −25.5 | bounded tail |
| Novel Ti/Tv at 90% tranche | **2.46** (hugging known 2.71) | 2.50 (known 2.63) | **1.69** (known 2.34) | novel ≈ 1.5–2.0, clearly below known |
| Model file size | 1.2 MB | 4.3 MB | 15.1 MB | — |
| Indel known Ti/Tv | 0.0000 (normal artifact for indels) | 0.0000 | 0.0000 | n/a |

The demo model's failure mode is visible in two numbers. First, the 100%-tranche
cutoff of −5889 means the model has **no discriminative power at the tail** — it is
separating noise from noise. Second, and more subtly: in a healthy callset the
*novel* variants (not in dbsnp) are enriched for artifacts, so **novel Ti/Tv falls
well below known Ti/Tv** — the teacher's run shows exactly that (1.69 vs 2.34).
The demo's novel Ti/Tv stays at 2.46, **hugging the known value** (2.71): the model
cannot tell its "worst" variants from its "best" ones, because at 0.5× even the
truth sites are supported by 1–2 reads. Note our own full run sits in between
(novel 2.50 vs known 2.63): with 55k on-target variants the model is usable but
still thin by GATK's own ≥30-sample guidance.

> **中文批注：** VQSR 在 14,782 个变异上"能跑完"不等于"能用"。三个铁证：① 训练用的 truth 位点只有 7,388 个（全量轮 28,324、老师 85,159）——模型是在极稀疏的标签上拟合的；② 100% tranche 的 minVQSLod = −5889——健康数据的尾部应在 −10 到 −30 之间，−5889 意味着模型在噪声里找规律；③ 最微妙的一点：健康数据里 novel Ti/Tv 应显著低于 known（老师的 1.69 vs 2.34 就是标准样子），而 demo 的 novel 2.46 几乎贴着 known 2.71——模型已无力区分"最差"和"最好"的变异，因为 0.5× 深度下连 truth 位点都只有 1–2 条读段支撑。这不是"VQSR 挂了"，而是"VQSR 拟合了一个没有统计意义的东西"——比挂了更危险，因为它不报错。

**Finding B (issue #1): the roadmap's §7.3 MergeVcfs duplicates every variant.**
`ApplyVQSR --mode SNP` does **not** output only SNPs — it outputs the **entire VCF**
with SNP-model scores applied (indels pass through unfiltered), and `--mode INDEL`
vice versa. Merging the two outputs therefore yields **every variant twice**:

| Callset | Raw joint records | After roadmap §7.3 MergeVcfs | Duplicated sites |
| --- | --- | --- | --- |
| demo | 14,782 | **29,564** | **100% (14,782/14,782)** |
| teacher `analysis_full` | 204,178 | **408,356** | **100%** |

The teacher's own reference output carries the same duplication — proving this is the
roadmap's bug, not an environment quirk. The exact damage (measured in both our runs):
every site appears **twice** with two different FILTER values — the SNP-mode record
carries the SNP-model verdict (`PASS` or `VQSRTrancheSNP…`), the INDEL-mode record
shows `.` (unfiltered passthrough) for SNPs and vice versa. Two consequences:
(a) the **raw record count doubles** (29,564; 110,246; 408,356) — so any downstream
row-counting (including ANNOVAR's avinput, §8) counts everything twice; (b) the
"PASS" count *happens* to stay correct (53,988 in our full round) because each site
has exactly one SNP-mode verdict — an accidental soft spot that masks the bug during
casual QC. One more subtlety surfaced in the full round: at sites where the SNP model
filtered the record, the *other* copy still shows `.` — so `SelectVariants
-select-expressions`-style "drop non-PASS" logic silently keeps the unfiltered twin.

## 7.3 The corrected workflow 矫正方案

Two standard fixes; **fix 1 is canonical GATK best practice** (sequential application —
each variant filtered exactly once, by its own model):

**Fix 1 — sequential ApplyVQSR (recommended):**

```bash
# SNP model → apply to full callset
gatk VariantRecalibrator ... -mode SNP -O snps.recal --tranches-file snps.tranches ...
gatk ApplyVQSR -V raw.vcf -O step1.vcf --mode SNP \
  --recal-file snps.recal --tranches-file snps.tranches --truth-sensitivity-filter-level 99.5

# INDEL model → apply to the SNP-filtered output
gatk VariantRecalibrator ... -mode INDEL -O indel.recal --tranches-file indel.tranches ...
gatk ApplyVQSR -V step1.vcf -O all.VQSR.vcf --mode INDEL \
  --recal-file indel.recal --tranches-file indel.tranches --truth-sensitivity-filter-level 99.0
```

**Fix 2 — split, filter, merge** (SelectVariants → ApplyVQSR per type → MergeVcfs),
if you want per-type VCFs on disk.

**And the honest alternative for THIS course's data (2 samples): hard filters.**
GATK's own recommendation for small callsets is `VariantFiltration`:

```bash
gatk VariantFiltration -V OC_blood_variants.vcf -O filtered.vcf \
  --filter-name "SNP_QD2"     --filter-expression "QD < 2.0 && SNP" \
  --filter-name "SNP_MQ40"   --filter-expression "MQ < 40.0 && SNP" \
  --filter-name "SNP_FS60"   --filter-expression "FS > 60.0 && SNP" \
  --filter-name "SNP_SOR3"   --filter-expression "SOR > 3.0 && SNP" \
  --filter-name "INDEL_QD2"  --filter-expression "QD < 2.0 && INDEL" \
  --filter-name "INDEL_FS200" --filter-expression "FS > 200.0 && INDEL" \
  --filter-name "INDEL_SOR10" --filter-expression "SOR > 10.0 && INDEL"
```

VQSR makes sense at cohort scale (30+ exomes); hard filters are the tool for 1–10
samples. **This is the single most transferable lesson of the whole step.**

> **中文批注：** 矫正方案记忆点：**顺序套用**（SNP 模型套完的输出，再套 INDEL 模型）——每条变异只被它自己类型的模型过滤一次；或者**先拆再合**（SelectVariants 按类型拆开 → 各自 ApplyVQSR → MergeVcfs）。而如果只有 1–10 个样本，GATK 官方建议干脆不用 VQSR，改用硬过滤（上面的 QD/MQ/FS/SOR 阈值就是 GATK 文档的标准值）。VQSR 的适用边界：约 30 个以上样本、或 callset 足够大——"能跑"≠"该用"。

## 7.4 Commands as we ran them 实跑命令说明

In both rounds we ran the roadmap's §7.1/7.2 VariantRecalibrator/ApplyVQSR verbatim
(recalibrators + mode-specific ApplyVQSR with the roadmap's 99.5/99.0 levels), and
§7.3 MergeVcfs verbatim — to *document* the duplication with measured numbers (14,782 →
29,564 demo; 55,123 → 110,246 full) rather than assert it. The corrected files for the final
deliverable are produced with Fix 1 (sequential ApplyVQSR) from the same recal models.

**Full-round numbers, as-run vs corrected:**

| Quantity | Roadmap §7.3 as-written (MergeVcfs) | Corrected (sequential ApplyVQSR, job 3056) |
| --- | --- | --- |
| Total records | **110,246** (= 2 × 55,123) | **55,123** (each site once) |
| Unique sites | 55,123 | 55,123 |
| PASS records | 53,988 | 53,988 |
| Sites carrying a `.` (unfiltered) FILTER | 53,988 | **0** |
| Filtered records | 56,258 (mostly unfiltered twins) | 1,135 (696 SNP-tranche + 427 + 12 INDEL-tranche) |
| ANNOVAR avinput lines from it | 110,246 | 55,123 |

The tranches of the full run (job 3051): SNP truth sites 28,324, minVQSLod −11.62 at
99.9% / −1111 at 100%; novel Ti/Tv 1.65 at the 100% tranche — the model behaves like a
real, if thin, callset model. The sequential application then removed exactly 1,135
of 55,123 records (2.1%) — a plausible filter-rate for a clean WES callset, in sharp
contrast to the demo's degenerate tranches.

## 7.5 Interpretation 结果解读

- The **tranches file** is the quality dashboard: `novelTiTv` per tranche should *fall*
  as you admit more variants (99 → 99.9 → 100). If it rises or wobbles, the model is
  overfitting — as it visibly does in the demo.
- `--truth-sensitivity-filter-level` semantics: "keep at least X% of truth-site variants";
  99.5% for SNPs (looser, SNP models are stronger) / 99.0% for indels (indel truth is
  scarcer; being slightly stricter compensates).
- VQSLOD distribution: `bcftools query -f '%VQSLOD\n' vcf | sort -g | awk ...` —
  a healthy run has a smooth distribution; our demo had a degenerate spike.
- **Do not compare raw QUAL between the rounds; compare tranche behavior.**

# 8. Functional annotation with ANNOVAR · 功能注释

## 8.1 Principle 原理

A VCF says *where* and *what*; annotation says *so what*. ANNOVAR's pipeline:
(1) **convert2annovar** translates VCF to a minimal tab format (`avinput`:
chr, start, end, ref, alt); (2) **annotate_variation / table_annovar** map variants
onto gene models (`refGeneWithVer` here — RefSeq with version suffixes) and optionally
onto dozens of other databases (population frequencies, deleteriousness scores, ClinVar).

The key output is `*.hg38_multianno.txt`: one row per (variant × transcript) with
`Func.refGeneWithVer` (exonic/intronic/UTR/splice/…), `ExonicFunc` (synonymous/
missense/nonsense/frameshift/…), and `AAChange.refGene` (protein change, e.g.
`BRCA1:NM_007294:exon11:c.3082G>A:p.Gly1028Arg` style coordinates).

## 8.2 Commands (roadmap 8.1–8.3, verbatim)

```bash
perl "$ANNOVAR/convert2annovar.pl" \
  -format vcf4old \
  "$OUT/2_variant_call/OC_blood.all.VQSR.vcf" \
  -includeinfo \
  -comment \
  -outfile "$OUT/3_annotation/OC_blood.all.avinput"

perl "$ANNOVAR/annotate_variation.pl" -geneanno -buildver hg38 \
  -dbtype refGeneWithVer \
  "$OUT/3_annotation/OC_blood.all.avinput" "$HUMANDB"

perl "$ANNOVAR/table_annovar.pl" \
  "$OUT/3_annotation/OC_blood.all.avinput" "$HUMANDB" \
  -buildver hg38 \
  -out "$OUT/3_annotation/OC_blood_anno" \
  -remove \
  -protocol refGeneWithVer \
  -operation g
```

## 8.3 Pitfall (issue #2): the "richer protocol" block fails — and exactly why

The roadmap's optional step, kept verbatim:

```bash
perl "$ANNOVAR/table_annovar.pl" "$OUT/3_annotation/OC_blood.all.avinput" "$HUMANDB" \
  -buildver hg38 -out "$OUT/3_annotation/OC_blood_anno_rich" -remove \
  -protocol refGeneWithVer,cytoBand,dbnsfp30a,exac03,gnomad_genome,avsnp147 \
  -operation g,r,f,f,f,f
```

This **cannot run on this cluster**: the installed humandb has
`dbnsfp54a`, `gnomad211_exome`, `gnomad41_exome`, `avsnp151` — not
`dbnsfp30a`, `exac03`, `gnomad_genome`, `avsnp147`. The error is
`Error: the database ... not found` (or "cannot open …txt"). Confirmed in both our runs;
the deviation logs record it.

**Correction** — the same intent, working protocol names:

```bash
perl "$ANNOVAR/table_annovar.pl" "$OUT/3_annotation/OC_blood.all.avinput" "$HUMANDB" \
  -buildver hg38 -out "$OUT/3_annotation/OC_blood_anno_rich" -remove \
  -protocol refGeneWithVer,cytoBand,gnomad41_exome,dbnsfp54a,avsnp151 \
  -operation g,r,f,f,f
```

(`gnomad41_exome` = gnomAD v4.1 exomes AF — the modern replacement for both `exac03`
and `gnomad_genome` for exome work; `dbnsfp54a` bundles functional predictions
including CADD/REVEL-class scores; `avsnp151` = dbSNP 151 rsIDs.)

> **中文批注：** ANNOVAR 的数据库名字**必须和 humandb 里实际存在的文件名完全一致**（`hg38_dbnsfp54a.txt` 等）。原版写的是 2020 年前后的老库名，而集群上装的是新库——这类"文档滞后于环境"的问题在生信教学环境里极其常见，遇到时先 `ls $HUMANDB | grep -iE "dbnsfp|gnomad|avsnp"` 看清楚有什么，再改 `-protocol`。另外注意每个 `-operation` 字母和 `-protocol` 里每个库一一对应（g=gene, r=region, f=filter/频率）。

## 8.4 Pitfall (issue #1 downstream): what the duplicated VCF does to ANNOVAR

Because the roadmap's §7.3 merged VCF contains every variant twice (§7.2), the
avinput and multianno outputs double-count everything:

| Callset | avinput lines | multianno rows |
| --- | --- | --- |
| demo | 29,564 (vs 14,782 real) | 29,602 |
| full | **110,246** (vs 55,123 real) | **110,284** |

Any "how many variants in gene X" count is inflated 2×. The corrected edition
annotates the **sequentially-VQSR-filtered VCF** (or the raw joint callset with hard
filters) — never the merged file from §7.3-as-written.

## 8.5 What happened 实跑结果

Demo: basic protocol produced `OC_blood_anno.hg38_multianno.txt` with columns
`Chr..AAChange.refGeneWithVer` (10 columns), 29,602 rows for 14,782 real variants.

Full round (job 3051, roadmap-verbatim): avinput **110,246** lines, basic multianno
**110,284** rows — the doubled inputs propagate exactly as predicted.

Full round **corrected** (job 3056): avinput **55,123** lines; rich protocol with the
**installed** database names (`refGeneWithVer,cytoBand,gnomad41_exome,dbnsfp54a,avsnp151`)
produced **55,165 rows × 146 columns** — including gnomAD v4.1 AFs (total + 8
sub-populations + faf95/99), SIFT4G, PolyPhen-2, MutationTaster, REVEL, CADD,
AlphaMissense, ClinPred, conservation scores (GERP++/phyloP/phastCons) and
`avsnp151` rsIDs. Of the 55,165 rows, **54,373 variants have gnomAD AF ≤ 0.01** —
i.e. 98.6% of this pair's germline variants are rare, exactly what two unrelated
individuals' exomes should look like (common variants are shared; rare ones are
private). The teacher's own annotation (`analysis_full/3_annotation`) used the older
`refGene,knownGene` protocol pair — functional/frequency interpretation was simply
not part of the reference result.

## 8.6 Tuning & beyond 调优与扩展

- `-nastring .` controls empty cells; `-polish -dot2comma`-type tidying matters for Excel users.
- Useful additional protocols available here: `clinvar_20221231`-style DBs (check humandb), `gerp++`, `spliceai` (if downloaded) — ClinVar is the one clinicians ask for first.
- `vcfanno`/`bcftools csq`/`VEP` are alternatives; ANNOVAR's niche is the flat-file simplicity and the course's muscle memory.
- **Rows ≠ variants**: one variant in 3 transcripts = 3 rows. Always `sort -u` by (chr,pos,ref,alt) before counting variants; count rows only when counting transcript-level effects.
- Convert back to VCF with `-vcfoutput` if you need annotated VCFs downstream (e.g. for reporting tools).

# 9. Somatic TMB with matched normal (Mutect2) · 体细胞突变与 TMB

## 9.1 Principle 原理

**Somatic = present in tumor, absent in the matched normal.** Mutect2's flow:
local de-novo assembly (like HC) in tumor **and** normal, then a series of
discriminative filters — population AF from a **germline resource**, cross-sample
evidence, orientation bias (FFPE/artifact), panel-of-normals recurrent artifacts —
each variant gets a posterior probability of being a somatic mutation.

**TMB (tumor mutational burden)** = (# somatic PASS variants) ÷ (target Mb). It's a
*rate*, so it only means something if (a) the denominator matches what you called
against (60.51 Mb here), and (b) the caller's sensitivity is depth-controlled
(~≥50× tumor for stable WES TMB; the measured 54.6×/91.3× pair is just above that floor —
see §11).

**Commands** (roadmap 9.1–9.3, verbatim):

```bash
gatk --java-options "-Xmx16g" Mutect2 \
  -R "$REFERENCE" \
  -I "$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam" \
  -I "$OUT/1_alignment/PBMC/blood_sorted.markdup.BQSR.bam" \
  --normal-sample blood \
  -L "$TARGETS" \
  --germline-resource "$HIGHCONF" \
  -O "$OUT/5_somatic_tmb/somatic.raw.vcf.gz" \
  --native-pair-hmm-threads 8 \
  --f1r2-tar-gz "$OUT/5_somatic_tmb/f1r2.tar.gz"

gatk LearnReadOrientationModel -I f1r2.tar.gz -O orientation_model.tar.gz

gatk FilterMutectCalls \
  -R "$REFERENCE" -V somatic.raw.vcf.gz \
  --stats somatic.raw.vcf.gz.stats \
  --ob-priors orientation_model.tar.gz \
  -O somatic.filtered.vcf.gz
```

## 9.2 Pitfalls & fixes 坑与矫正

- **(p1, issue #5) Missing Panel of Normals.** The cluster **ships one**:
  `$REF/1000g_pon.hg38.vcf.gz` (+ `.tbi`). Add `--pon $REF/1000g_pon.hg38.vcf.gz`
  and FilterMutectCalls gains a whole filter dimension (recurring capture/alignment
  artifacts). This is the single highest-value one-flag improvement to §9.
- **(p2) Germline resource semantics:** Mutect2 wants a **population AF** file to
  down-weight common germline; `1000G_phase1.snps.high_confidence` carries AF in its
  INFO and *works*, but gnomAD AF-only (`af-only-gnomad.vcf.gz`) is the modern
  recommendation — richer, newer, and exome-matched.
- **(p3, cosmetic) `somatic.pass.vcf.gz` is written but never used** by §9.3, which
  re-derives PASS from the filtered VCF (same result; keep whichever, but be consistent).
- **(p4) `--normal-sample blood`** must match the BAM's `SM` tag exactly (`SM:blood`).
  If you rename samples midway, Mutect2 fails with "FATAL error: sample not found".

## 9.3 What happened 实跑结果

**Demo round (0.5×): raw = filtered = PASS = 0, TMB = 0.00/Mb.**
This is the expected and *pedagogically perfect* outcome: at ~0.5× the tumor, there is
no statistical way to distinguish a somatic mutation from a germline variant seen by
chance — Mutect2's model correctly refuses to call anything. Zero is not a bug; it is
the honest answer.

**Full round (~54.6× usable tumor depth): raw = 2,949, filtered = 2,949, PASS = 108 → TMB = 108 / 60.5079 Mb = 1.78 mut/Mb.**

The filter breakdown reads like a textbook of why somatic calling needs the
learning filters — of the 2,841 non-PASS records (combinatorial FILTER strings):

| Dominant filter | Records | What it means |
| --- | --- | --- |
| `normal_artifact` (in 60%+ of filtered strings) | — | seen (weakly) in the normal at higher depth than this capture's noise floor — mostly germline sites with lopsided allele fractions |
| `weak_evidence` | — | tumor LOD below threshold — the 54.6× usable depth still caps sensitivity |
| `slippage` | — | homopolymer/STR-adjacent indel artifacts |
| `germline` | — | population AF says germline |
| `multiallelic` / `strand_bias` / `orientation` / `haplotype` / `clustered_events` / `map_qual` / `base_qual` | — | the full modern filter battery |

108 PASS variants with AFs from 0.056 to 0.48 (median ~0.4 — consistent with
heterozygous clonal mutations in a high-purity tumor).

**Interpretation 结果解读.**

- A realistic epithelial ovarian tumor TMB is **~1–3 mut/Mb** (HGSOC is a
  low-TMB, copy-number-driven disease). If you get 50+ mut/Mb from a quiet tumor
  sample, suspect germline leakage (bad normal matching, sample swap) or artifacts.
- Filter breakdown (`bcftools query -f '%FILTER\n' | sort | uniq -c`) tells you *why*
  variants died: a pile of `germline` filter = your "normal" isn't matching the tumor;
  piles of `orientation` = FFPE artifact; `panel_of_normals` (if you add the PoN) =
  capture artifacts.
- CHIP in the normal (blood DNMT3A/TET2/ASXL1) can masquerade as somatic in the
  tumor direction if the tumor purity is low — keep it in your differential.

> **中文批注：** TMB 的分母（60.51 Mb）必须与 `-L` 的区域一致——用全外显子组 38 Mb 还是全部捕获区域 60.5 Mb，TMB 能差 60%。第二个要点：demo 轮 TMB=0 不是流程坏了，而是"0.5× 深度下体细胞检测没有统计功效"的诚实表达。全量轮实测可用深度 54.6×，刚好在稳定 TMB 的临线之上——这也是为什么全量轮只报了 108 个 PASS 体细胞变异（HGSOC 本身也是低 TMB 肿瘤）。第三个要点：HGSOC（高级别浆液性卵巢癌）是典型的低 TMB、高 CNV 肿瘤——体细胞点突变少、拷贝数变化多，这正好衔接 §10 的 CNV 检测为什么对这类肿瘤特别重要。

## 9.4 Tuning & beyond 调优与扩展

- `--max-mnp-distance 0` for strict SNP/indel output; default 1 also emits MNPs.
- Tumor-in-normal mode (`--tumor-sample` in both inputs) when you have no normal: weaker, but possible.
- Absolute calling beyond TMB: combine with CNVkit's purity/ploidy estimate (`cnvkit.py call --purity`) → Mutect2's `--tumor-lod`-style thresholds tuned per sample.
- Report-grade filtering: `bcftools filter -i 'FILTER="PASS" && FORMAT/AF[0]>0.05'` for a confidence cut; always state VAF and depth alongside every somatic call.

# 10. CNV detection with CNVkit · 拷贝数变异

## 10.1 Principle 原理

CNVkit compares the tumor's read **depth** per bin (on- and off-target) against a
**reference** (here: the matched normal), converts to log2 ratio, segments (CBS), and
calls integer copy number. The mode flag encodes *how the panel was prepared*:

- `-m hybrid` — hybrid-capture (SureSelect/Xgen/…): **off-target reads exist** and are
  usable as a baseline → CNVkit builds antitarget bins.
- `-m amplicon` — amplicon PCR: everything is on-target by construction; no antitarget.

This panel is **S07604514 = Agilent SureSelect Human All Exon V6** — a *hybrid
capture* design. The roadmap's command uses `-m amplicon` (issue #4).

## 10.2 Commands (roadmap 10.2, verbatim) and what it produced

```bash
cnvkit.py batch \
  "$TUMOR" -n "$NORMAL" \
  -f "$REFERENCE" \
  -t "$CNVOUT/targets.bed" \
  -m amplicon \
  --output-reference "$CNVOUT/reference.cnn" \
  -d "$CNVOUT/cnvkit_out" \
  --scatter --diagram -p 8
```

**Demo round:** it *ran* (no crash) — and produced textbook-garbage: segments with
`log2` up to **+11.6** and called `cn` values of 324, 1601, 6345 (copy number cannot
exceed ~2× ploidy in a diploid sample; anything >10 is a "the reference is empty"
mathematical artifact), plus a 169-Mb "loss" spanning half of chr1 at depth ~0.5×.
Exactly what a 0.5× tumor + 0.5× single-normal reference must produce: the reference's
bin variances are pure noise, so every bin deviates "significantly".

**Full round, roadmap-verbatim (`-m amplicon`, job 3053):** no crash, plausible-looking
numbers — and the tell is still there: the `*.antitargetcoverage.cnn` files are
**31 bytes / 0 bytes** (headers only). Amplicon mode simply *has no antitarget
baseline*; everything is computed from target bins alone. The call set (290 segments:
cn=0 ×1, cn=1 ×58, cn=2 ×191, cn=3 ×39) looks reasonable at a glance — but the
weights and diagnostically important arm-level events are all muted, because without
off-target bins the log2 ratios lose their genome-wide anchor.

**Full round, corrected (`-m hybrid` + access file, job 3055):** real antitarget
coverage (1.4-Mb antitarget BED, 55-Mb of binned off-target space) → 357 segments
(cn=0 ×1, cn=1 ×69, cn=2 ×202, cn=3 ×35) with **high-weight arm-level events**:

| Event (top by weight) | log2 | cn | Weight | Reading |
| --- | --- | --- | --- | --- |
| **chr9: 78–138 Mb (9q arm) single-copy loss** | −0.76 | **1** | 7,416 | classic HGSOC 9q loss (spanning CDKN2A-adjacent region, Notch/TGF-β loci) |
| chr8: 13–124 Mb gain | +0.49 | 3 | 7,172 | whole-arm-ish chr8 gain incl. **8q24 (MYC)** — the signature HGSOC amplifier |
| chr5: 75–181 Mb gain | +0.46 | 3 | ~4,000 | chr5q+ partial gain |
| chr20 gain (two large segments) | +0.51 | 3 | ~2,000 | chr20 amp — common in OC lines |
| chr21 gain | +0.44–0.50 | 3 | ~1,000 | trisomy-21-like gain |
| **chrY loss** | −2.05 | **0** | 297 | consistent with a female patient (OC) — tumor and normal both "lost" chrY |
| chr14q / 16p / 17p focal losses | −0.33 to −0.75 | 1 | ~200–300 | sub-arm deletions |

That is a **karyotypically coherent HGSOC picture** (low-TMB, arm-level CNV-driven —
see §9's TMB=1.78 for the other half of the story). The amplicon-mode run shows
*none* of these at comparable weight — the concrete, measured cost of one wrong flag.

## 10.3 Pitfalls & fixes 坑与矫正 (issue #4)

1. **Mode:** use `-m hybrid` for this panel. Amplicon mode skips antitarget
   generation, so all off-target reads are *discarded* — losing the genome-wide
   baseline that stabilizes log2 ratios.
2. **Access file:** with hybrid mode, pass `-g $REF/access-5kb.hg38.bed` (it's on the
   cluster, ready) so antitarget bins exclude centromeres/gaps.
3. **Gene names:** the BED is 3-column (no `gene` column), so `--diagram`/scatter
   lack gene labels. `cnvkit.py import-rna`-style tricks don't apply; simplest is
   `--annotate $REF/...refGene.txt`-based annotation (`cnvkit.py batch --annotate`).
4. **Reference quality:** a single-normal reference amplifies that normal's noise.
   A pooled PoN reference (`--normal a.bam b.bam c.bam`) is the production answer;
   for the course, matched-normal is an acceptable teaching choice — *with the mode fixed*.

**Corrected command (hybrid mode + access file):**

```bash
cnvkit.py batch \
  "$TUMOR" -n "$NORMAL" \
  -f "$REFERENCE" \
  -t "$CNVOUT/targets.bed" \
  -g "$REF/access-5kb.hg38.bed" \
  -m hybrid \
  --output-reference "$CNVOUT/reference.cnn" \
  -d "$CNVOUT/cnvkit_out_hybrid" \
  --scatter --diagram -p 8
```

> **中文批注：** 为什么模式选错"能跑但结果不可信"？amplicon 模式的假设是"所有读段都来自靶区"，于是 off-target 读段全部弃用；而杂交捕获数据里 off-target 读段（通常 30–50%）恰恰是拷贝数分析的"内参基线"。丢了它们，log2 比值的分母就只剩噪声。demo 轮里 cn=1601、6345 这种荒谬值就是这么来的。矫正三件套：`-m hybrid` + `-g access文件` + （可选） BED 加基因名列。对 HGSOC 这种 CNV 驱动的肿瘤，这步修好了才算真正"上了正菜"。

## 10.4 Interpretation 结果解读

- `.cnr` = per-bin log2; `.cns` = segmented; `.call.cns` = integer CN calls
  (`cn` column: 0,1,2,3…; `weight` = bin support).
- In .call.cns, `log2` ±0.3 is the usual "gain/loss" noise floor for WES;
  arm-level events in HGSOC (chr8q gain, chr1q gain, chr16/17 loss, PTEN/RB1 loss)
  are the biologically expected pattern — and the corrected full-round run indeed found
  9q single-copy loss, chr8 gain (cn=3, spanning MYC), chr20/21 gains and chrY loss
  at high weight (§10.2).
- `--diagram` produces the genome diagram PDF; `--scatter` the per-chromosome log2
  scatter. CNVkit's `metrics` and `segmetrics` give per-segment robust stats.
- Sex chromosomes: male/female mismatch between tumor and normal inflates chrX/Y —
  CNVkit normally handles it if sex is declared (`--sex`); check the log.

## 10.5 Tuning & beyond 调优与扩展

- Purity/ploidy: `cnvkit.py call --purity 0.5` integrates with the TMB/AF interpretation.
- `cnvkit.py export seg/vcf/theta` bridges to DNAseq tools, IGV, ABSOLUTE/THetA2.
- GISTIC on a cohort of .cns/seg files → recurrent amplifications/deletions; that is
  the "cohort version" of this step.
- For immune-oncology course tie-ins: CNV of PD-L1 (CD274), copy-loss of HLA locus,
  and CCNE1 amplification are the classic HGSOC/immuno-evasion stories to grep for in
  the output.

# 11. Interpretation checklist · 结果解读检查清单

The original roadmap closes with an 11-point checklist. Here it is again — now with
**actual answers from both runs**, and the one missing metric (on-target depth) added.

| # | Question | Demo (0.5×) answer | Full (~55×/91× usable on-target) answer |
| --- | --- | --- | --- |
| 1 | Do FASTQ reports pass Q30 thresholds? | Yes: 97.0% / 97.8% (post-trim) | Yes: 97.0% / 97.8% (post-trim; raw 95.5%/96.9%) |
| 2 | Is mapping rate >95%? | 100% | 100% — but note: this **dataset ships pre-filtered** (only 311/79.6M unmapped in the full round); fresh lab data is 96–99% |
| 3 | Properly paired >90%? | OC 83.3%, PBMC 99.5% | OC **83.2%** (real tumor discordance — 15.8% mates on other chr), PBMC 99.4% |
| 4 | Duplication <20%? | ~0.03–0.06% (subsample artifact) | OC **7.08%**, PBMC **10.77%** (teacher's: 7.39/11.01) |
| 5 | **On-target depth sufficient?** *(the metric the roadmap forgot)* | **0.3×/0.5× usable, median 0 — NO** | OC **54.6×** (median 57; 87.4% ≥30×), PBMC **91.3×** (median 96; 89.6% ≥30×) — borderline-adequate tumor, good normal |
| 6 | Insert size within capture design range? | peak **150 bp** (short!) | 150 bp — 25% of usable bases lost to mate overlap (HsMetrics `PCT_EXC_OVERLAP`) |
| 7 | Are known sites recovered (dbsnp overlap)? | 7,388 truth sites | **28,324** truth sites (teacher's genome-wide: 85,159) |
| 8 | Ti/Tv of PASS variants? | joint 2.59; PASS 2.60 — novel Ti/Tv **hugs** known (2.33 vs 2.65): no discrimination | PASS **2.48** (known 2.49) — healthy range |
| 9 | Het/Hom ratio ~1.5 for exome? | **0.14** — genotypes collapse at 0.5× (single-read "homs") | **1.45** — healthy |
| 10 | Somatic count plausible for tumor type? | 0 (no power) | **108 PASS / TMB 1.78 mut/Mb** — squarely in HGSOC range |
| 11 | Do CNV calls make karyotypic sense? | No — pure noise (§10.2) | **Yes** (corrected hybrid mode): 9q cn=1 loss, chr8 gain cn=3, chrY loss — coherent HGSOC karyotype |

**The missing metric — add this to every WES run (issue #9):**

```bash
# CollectHsMetrics wants a Picard interval_list with a SAM header — convert the BED first
gatk --java-options "-Xmx4g" BedToIntervalList \
  -I "$TARGETS" -O "$OUT/qc/targets.interval_list" -SD "$REF/hg38.dict"

gatk --java-options "-Xmx8g" CollectHsMetrics \
  -I "$OUT/1_alignment/OC/OC_sorted.markdup.BQSR.bam" \
  -O "$OUT/1_alignment/OC/OC.hs_metrics.txt" \
  --TARGET_INTERVALS "$OUT/qc/targets.interval_list" \
  --BAIT_INTERVALS "$OUT/qc/targets.interval_list" \
  --REFERENCE_SEQUENCE "$REFERENCE"   # no separate bait list on the cluster; targets≈baits for teaching
```

*(Two practical notes learned the hard way in this very round: plain BEDs are rejected
— `Interval list file must contain header` — hence the BedToIntervalList step; and the
sequence dictionary lives at `hg38.dict`, not `hg38.fa.dict`.)*

The rows that matter: `MEAN_TARGET_COVERAGE`, `PCT_TARGET_BASES_30X` (or 50X/100X
columns), `PCT_SELECTED_BASES` (on-target + near-bait = capture specificity),
`FOLD_ENRICHMENT`, and — surprisingly informative — `PCT_EXC_OVERLAP` / `PCT_EXC_DUPE` /
`PCT_EXC_OFF_TARGET`: the three places your paid-for depth disappears.

Measured, both rounds (all four files under `wes_run_full_20260920/qc/`):

| Sample | MEAN_TARGET_COV | ≥30X | ≥100X | PCT_SELECTED | ZERO_CVG | EXC_DUPE | EXC_OVERLAP |
| --- | --- | --- | --- | --- | --- | --- | --- |
| demo OC | 0.3× | 0% | 0% | 79.9% | 59.4% | 0.03% | 26.8% |
| demo PBMC | 0.5× | 0% | 0% | 90.8% | 46.6% | 0.06% | 24.9% |
| full OC | **54.6×** | 87.4% | 74.2% | 79.9% | 8.4% | 7.3% | 25.6% |
| full PBMC | **91.3×** | 89.6% | 87.0% | 90.7% | 8.1% | 11.0% | 23.0% |

> **中文批注：** 检查清单的 11 个问题里，最容易被"绕过去"的是 #5（靶区深度够不够）——因为整条路线图里没有任何命令真正测量过它。上表就是实测结果：全量轮肿瘤 54.6×、正常 91.3×，80–100× 的临线刚好够用；demo 轮中位深度是 **0**。另外注意两个实测细节：① 这套文库插入片段峰值只有 150bp，双端重叠吃掉约 25% 的可用碱基——这是"付了 100× 的钱、拿到 55× 的深度"的主要原因；② 零覆盖靶区 8.4%——AllTracks 版 BED 含大量难捕获区域，这属于捕获设计的固有属性，不是实验失败。

# 12. File summary · 产物清单

What each run directory contains when everything is done (paths as produced by the
verbatim pipeline; sizes from the full round, demo sizes in parentheses where notable):

| Artifact (full round path under `wes_run_full_20260920/`) | Size | Notes |
| --- | --- | --- |
| `0_fastq/OC/OC.fastp.html/.json` (+ clean FASTQs) | ~0.5 G clean FASTQ pairs | Q30 97.0% post-trim |
| `1_alignment/OC/OC_sorted.bam` | 4.0 GB | piped `bwa\|sort` — no 33-GB SAM ever existed (issue #7) |
| `1_alignment/OC/OC_sorted.markdup.bam` | 4.7 GB | + `.bai` via `--CREATE_INDEX` |
| `1_alignment/OC/OC_sorted.markdup.BQSR.bam` | 8.8 GB | ~1.9× the markdup BAM — recompression quirk; budget for it |
| `1_alignment/OC/OC.g.vcf` | 209 MB | GVCF per sample (PBMC: 59 MB) |
| `2_variant_call/combined.g.vcf` | 264 MB | 2-sample combined GVCF |
| `2_variant_call/OC_blood_variants.vcf` | 14 MB | **55,123 records** — the raw joint callset |
| `2_variant_call/OC_blood.snps.recal` + `.tranches` (+indel) | 4.3 MB + 0.6 MB | VQSR models |
| `2_variant_call/OC_blood.all.VQSR.vcf` | 30 MB | **110,246 records — the §7.3 duplication, as-written** |
| `2_variant_call/OC_blood.all.sequentialVQSR.vcf` | 16 MB | **55,123 records — the corrected output** |
| `3_annotation/OC_blood_anno.hg38_multianno.txt` | 10 MB | basic protocol, on the duplicated VCF (110,284 rows) |
| `3_annotation/OC_blood_anno_rich_installed.hg38_multianno.txt` | 27 MB | **corrected rich protocol: 55,165 rows × 146 cols** |
| `5_somatic_tmb/somatic.raw.vcf.gz` → `somatic.filtered.vcf.gz` → `somatic.pass.vcf.gz` | 0.3 → 0.35 → 0.02 GB | 2,949 raw → 108 PASS |
| `4_cnv/cnvkit_out/` (amplicon, verbatim) | ~47 MB | empty antitarget files — the tell (§10.2) |
| `4_cnv/cnvkit_out_hybrid/` (corrected) | ~55 MB | real antitarget bins; `--scatter`/`--diagram` PDFs in both |
| `qc/*.hs_metrics.txt` | 4 files | the §11 coverage numbers, both rounds |

(Demo round: same layout, everything ~50–100× smaller — e.g. joint VCF 2.8 MB,
BQSR BAM 54 MB.)

# 13. Summary · 总结与两轮对比

**What the two rounds taught, in one table:**

| Aspect | Demo round (0.5×) | Full round (54.6×/91.3× usable on-target) |
| --- | --- | --- |
| Purpose | Pipeline mechanics | Real biology |
| Alignment | 100% mapped (dataset pre-filtered), OC properly-paired 83.3% | identical profile: OC 83.2% (real tumor discordance), PBMC 99.4% |
| Germline joint callset | 14,782 records | **55,123 records** (2-sample, on-target; teacher's genome-wide 1-sample: 204,178) |
| VQSR | ran; truth sites 7,388; minVQSLod@100% = −5889; novel Ti/Tv hugs known (2.33 vs 2.65) | truth sites 28,324; tranches behave like a real (if thin) model; sequential correction removes 2.1% |
| §7.3 duplication bug | 29,564 = 2 × 14,782 | **110,246 = 2 × 55,123** (corrected: 55,123, zero unfiltered twins) |
| ANNOVAR | basic OK, rich protocol fails (missing DBs) | verbatim rich still fails; **corrected names: 55,165 rows × 146 cols** incl. REVEL/CADD/AlphaMissense |
| Somatic / TMB | 0 — no statistical power | **108 PASS / TMB 1.78 mut/Mb** (HGSOC-plausible) |
| CNVkit | cn=6345 absurdities | verbatim amplicon: empty antitarget, muted events; **corrected hybrid: 9q loss, chr8 gain — coherent HGSOC karyotype** |
| Wall time | ~25 min (1 job) | 3 h 42 m longest chain (8 jobs, dependency-linked; prep 2:20+3:42 in parallel) |

**The four lessons to carry away:**

1. **"It ran" ≠ "it's right."** The VQSR duplication bug (§7.2) and the CNVkit mode
   error (§10.3) both run to completion. The single most important skill this course
   builds is *counting your own outputs* (variants, rows, sites) and knowing what
   numbers are plausible.
2. **Depth gates everything.** 0.5× is enough to exercise the machinery, and it is
   *honest* about its limits (somatic=0, CNV=noise). Interpretation belongs to the
   full-depth round.
3. **Tools encode assumptions** (capture mode, sample naming, database names,
   "VQSR needs 30 samples"). Reading the manual after the first error is normal;
   reading it before the first run is professional.
4. **Corrected, the pipeline is a strong teaching skeleton**: the roadmap's
   step order, RG discipline, GVCF workflow, and germline-vs-somatic separation are
   all best-practice — the fixes are surgical, not structural.

> **中文批注：** 两轮对比最大的价值是把"流程问题"和"数据问题"彻底分开：demo 轮暴露的 0、噪声、荒谬 CN 值是**数据深度问题**（换全量数据自然消失）；VQSR 重复、ANNOVAR 库名、CNVkit 模式是**流程本身的问题**（换什么数据都在，连老师的参考结果里都有）。把这两类问题分清楚，你就掌握了"排错"的第一层功力。

## Reproduction · 复现说明

Both rounds were run as SLURM jobs (partition `cn-long`, QoS `bjx131cnl`,
account `bjx131_g1`, `--mem=0`, 8–16 cores per job):

| Round | Directory | Scripts |
| --- | --- | --- |
| demo | `/lustre1/user/bjx131_pkuhpc/wes_run_cluster_20260914` | `scripts/wes_demo_steps3to10.sbatch` (fastp was already done in a prior session) |
| full (verbatim) | `/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920` | `scripts/full_prep_OC.sbatch`, `scripts/full_prep_PBMC.sbatch`, `scripts/full_joint.sbatch`, `scripts/full_somatic.sbatch`, `scripts/full_cnv.sbatch` (joint/somatic/cnv submitted with `--dependency=afterok:<prep jobs>`) |
| full (corrected extras) | same | `scripts/full_hsmetrics.sbatch` (§11), `scripts/full_cnv_hybrid.sbatch` (§10.3), `scripts/full_fix_vqsr.sbatch` (sequential ApplyVQSR + rich ANNOVAR with installed DB names) |

Each run directory contains `log/deviations.log` (every expected-failure and deviation
with timestamps) and `log/summary_*.log` (the numbers quoted throughout this edition).
The full round keeps the roadmap's literal SAM→BAM→sort chain *only in the demo round*;
the full round used the piped variant (issue #7, logged).

**Roadmap fidelity statement:** every command in steps 3–10 was executed exactly as
printed in the roadmap (same flags, same file names, same order). The only wrapper
additions are: SLURM headers, environment activation, preflight existence checks,
expected-failure `if/else` guards (which log-and-continue instead of dying), and the
summary blocks. The single algorithmic deviation is the `bwa | sort` pipe in the full
round, which changes no output content.

# Appendix A. Issue registry · 问题清单（完整版）

Each issue: **where** (roadmap section), **what** (symptom), **evidence** (measured),
**fix** (what this edition recommends), **severity**.

### A1. VQSR MergeVcfs duplication — **Major**

- Where: §7.3. What: `OC_blood.all.VQSR.vcf` contains every variant twice.
- Evidence: demo 29,564 records / 14,782 unique sites; our full round 110,246 / 55,123;
  teacher 408,356 / 204,178. Mechanism: ApplyVQSR (either mode) emits the whole callset
  (the other variant type passes through with FILTER `.`); MergeVcfs concatenates
  two whole callsets. Each site ends up as one SNP-verdict record + one `.` record —
  the PASS count *happens* to stay right, which masks the bug during casual QC.
- Fix: sequential ApplyVQSR (SNP → INDEL), or SelectVariants split → filter → merge.
  Measured result of the fix (job 3056): 55,123 records / 55,123 unique sites / 0
  unfiltered twins; filter rate 2.1%.
- Severity: Major — silently corrupts every downstream count and annotation.

### A2. ANNOVAR rich protocol names unavailable — **Major**

- Where: §8.3. What: `dbnsfp30a`, `exac03`, `gnomad_genome`, `avsnp147` absent from
  `/lustre1/share/annovar_new/humandb` (which has `dbnsfp54a`, `gnomad211_exome`,
  `gnomad41_exome`, `avsnp151`).
- Evidence: table_annovar fails in both runs; `ls humandb` listing.
- Fix: use the installed names (§8.3). Severity: Major (hard failure).

### A3. VQSR on an undersized callset — **Major (statistical)**

- Where: §7. What: with 2 exomes / ~15k variants (demo) the Gaussian mixture is
  unidentifiable; tranches become meaningless.
- Evidence (all three callsets measured, §7.2): truth sites 7,388 (demo) vs 28,324
  (our full) vs 85,159 (teacher); minVQSLod at the 100% tranche −5889.9 (demo) vs
  −1111 (full) vs −25.5 (teacher); and the decisive pattern: healthy callsets show
  novel Ti/Tv well *below* known (teacher: 1.69 vs 2.34) while the demo's novel
  Ti/Tv **hugs** known (2.46 vs 2.71) — a model that cannot separate its worst from
  its best variants.
- Fix: for ≤10 samples use VariantFiltration hard filters (§7.3); VQSR is a cohort tool
  (≈30+ samples). Severity: Major — silently wrong filtering.

### A4. CNVkit `-m amplicon` on a hybrid-capture panel — **Major**

- Where: §10.2. What: S07604514 is a SureSelect V6 hybrid capture; amplicon mode
  discards all off-target reads and builds no antitarget baseline.
- Evidence (both rounds measured): demo produced cn up to 6345 and 169-Mb "losses"
  at 0.5×; the full-round verbatim run produced **empty antitarget coverage files
  (31 bytes)** and muted arm-level events, while the corrected hybrid run on the same
  BAMs found 9q single-copy loss, chr8 gain (cn=3), chrY loss at high weight (§10.2).
  The pipeline completes either way.
- Fix: `-m hybrid` with `-g access-5kb.hg38.bed` (both on the cluster). Severity: Major
  — results unusable as-is.

### A5. Mutect2 without a Panel of Normals — Medium

- Where: §9.1. What: PoN filters recurring capture/alignment artifacts; the cluster
  ships `1000g_pon.hg38.vcf.gz` but the roadmap doesn't use it.
- Fix: `--pon $REF/1000g_pon.hg38.vcf.gz`. (Also: gnomAD AF-only is the recommended
  germline resource; 1000G high-conf works but is dated.)

### A6. Dataset path inconsistency — Medium

- Where: §0.2, README, old vs new roadmap. What: three different `$DATA` targets
  (`OC_WES` deprecated downsample / `data_new` new downsample / an **empty**
  `analysis/0_preprocess` in the old HTML).
- Fix: this edition pins every dataset by absolute path with role labels (§0.2).

### A7. 33 GB intermediate SAM per sample — Medium (efficiency)

- Where: §3.1–3.3. Fix: `bwa | sort` pipe (full round). Output identical; disk and
  ~30–40% wall time saved. (Demo round kept the literal chain for fidelity.)

### A8. TMB counts source inconsistency — Minor

- Where: §9.2/9.3: `somatic.pass.vcf.gz` written, but TMB counted from
  `somatic.filtered.vcf.gz` PASS subset. Same number; keep one convention.

### A9. No on-target depth metric — Minor

- Where: §11 checklist asks, no command measures. Fix: CollectHsMetrics (§11).

### A10. Naming inconsistencies — Minor

- `SM:blood` vs "PBMC" directories (functional but confusing); teacher's
  `2_vairant_call` directory typo (reference side); `vcf4old` format flag in ANNOVAR
  (works; `vcf4` is the modern name).

### A11. Runtime environment quirks (not roadmap bugs, but you will hit them) — Info

- SLURM on this cluster: `--account=bjx131_g1 --qos=bjx131cnl` mandatory; `--mem=0`
  required (RealMemory=1 misconfiguration); conda activate before `set -u`.
- GATK bundle `.idx` files all present — but check on any new cluster.
- CollectHsMetrics wants a Picard interval_list (run BedToIntervalList on the BED
  first) and the dictionary is `hg38.dict`, not `hg38.fa.dict`.

### A12. The teacher's reference result deviates from the printed roadmap — **Major
(for comparability; not a crash)**

- Where: `analysis_full/2_vairant_call` (the reference everyone diffs against).
- Evidence (read from the VCF headers and the files themselves):
  1. **Single-sample "joint" call:** `combined.g.vcf` has only sample `OC` — the
     blood gVCF's RG `SM` tag was set to `OC` (sample-name collision), so
     CombineGVCFs collapsed both samples into one. The printed roadmap clearly
     intends a 2-sample joint call (OC + blood).
  2. **No `-L` anywhere:** HaplotypeCaller/GenotypeGVCFs ran genome-wide — 72% of
     the 204,178 sites lie outside the capture targets, many with shallow,
     low-quality coverage.
  3. §7.3's MergeVcfs duplication is present in the reference output as well
     (408,356 = 2 × 204,178 — issue A1).
- Consequence: "204,178 vs my 55,123" is **not** a depth comparison — it's
  (genome-wide × 1 sample × doubled by the merge bug) vs (on-target × 2 samples).
  Diff callsets only after reading their `##GATKCommandLine` headers.
- Lesson for the course: the RG `SM` tag discipline the roadmap itself preaches (§2)
  is exactly what silently bit its own reference run. This is the best argument for
  step 0's "check your BAM header" habit.
- Fix: nothing to fix in the roadmap text (it prescribes `-L` and correct SM tags);
  this entry exists so students comparing against `analysis_full` know what the
  differences are and why.

# Appendix B. One-page corrected pipeline · 一页矫正版命令速查

```bash
# 0) env
source /lustre1/share/miniconda3/etc/profile.d/conda.sh && conda activate wes
export SHARE=/lustre1/share REF=$SHARE/references
export TARGETS=$REF/S07604514_AllTracks_V6_60_hg38.bed REFERENCE=$REF/hg38.fa
export DBSNP=$REF/dbsnp_138.hg38.vcf INDELS=$REF/Mills_and_1000G_gold_standard.indels.hg38.vcf
export HIGHCONF=$REF/1000G_phase1.snps.high_confidence.hg38.vcf
export PON=$REF/1000g_pon.hg38.vcf.gz ACCESS=$REF/access-5kb.hg38.bed
export ANNOVAR=$SHARE/annovar_new HUMANDB=$ANNOVAR/humandb

# 2) QC
fastp -i R1.fq.gz -I R2.fq.gz -o R1.clean.fq.gz -O R2.clean.fq.gz \
  --detect_adapter_for_pe --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
  --qualified_quality_phred 20 --unqualified_percent_limit 40 --n_base_limit 5 \
  --length_required 50 --thread 8 --html S.html --json S.json

# 3) align (piped)
bwa mem -t 8 -M -R "@RG\tID:S\tPL:illumina\tSM:S\tLB:libS" $REFERENCE R1.clean.fq.gz R2.clean.fq.gz \
  | samtools sort -@ 8 -m 4G -o S_sorted.bam -
samtools index -@ 8 S_sorted.bam
samtools flagstat -@ 8 S_sorted.bam > S.flagstat.txt
gatk BedToIntervalList -I $TARGETS -O targets.interval_list -SD $REF/hg38.dict
gatk CollectHsMetrics -I S_sorted.bam -O S.hs_metrics.txt \
  --TARGET_INTERVALS targets.interval_list --BAIT_INTERVALS targets.interval_list \
  --REFERENCE_SEQUENCE $REFERENCE   # NEW: on-target depth (§11)

# 4) markdup
gatk MarkDuplicates -I S_sorted.bam -O S_sorted.markdup.bam -M S_markdup_metrics.txt --CREATE_INDEX true

# 5) BQSR
gatk BaseRecalibrator -R $REFERENCE -I S_sorted.markdup.bam -L $TARGETS \
  --known-sites $DBSNP --known-sites $INDELS --known-sites $HIGHCONF -O S.recal.table
gatk ApplyBQSR -R $REFERENCE -I S_sorted.markdup.bam --bqsr-recal-file S.recal.table -O S.BQSR.bam

# 6) germline GVCF workflow
gatk HaplotypeCaller -ERC GVCF -R $REFERENCE -I S.BQSR.bam -D $DBSNP -L $TARGETS \
  --native-pair-hmm-threads 8 -O S.g.vcf
gatk CombineGVCFs -R $REFERENCE -V OC.g.vcf -V blood.g.vcf -O combined.g.vcf
gatk GenotypeGVCFs -R $REFERENCE -V combined.g.vcf -D $DBSNP -O joint.vcf

# 7) filtering — small callset: hard filters (VQSR needs ~30+ samples)
gatk VariantFiltration -V joint.vcf -O joint.filtered.vcf \
  --filter-name SNP_QD2 --filter-expression "QD < 2.0 && SNP" \
  --filter-name SNP_MQ40 --filter-expression "MQ < 40.0 && SNP" \
  --filter-name SNP_FS60 --filter-expression "FS > 60.0 && SNP" \
  --filter-name SNP_SOR3 --filter-expression "SOR > 3.0 && SNP" \
  --filter-name INDEL_QD2 --filter-expression "QD < 2.0 && INDEL" \
  --filter-name INDEL_FS200 --filter-expression "FS > 200.0 && INDEL" \
  --filter-name INDEL_SOR10 --filter-expression "SOR > 10.0 && INDEL"
# (cohort scale instead: VariantRecalibrator SNP → ApplyVQSR → VariantRecalibrator INDEL
#  → ApplyVQSR on the SNP output — sequential, never MergeVcfs of two whole callsets)

# 8) ANNOVAR (installed databases)
perl $ANNOVAR/convert2annovar.pl -format vcf4 joint.filtered.vcf -includeinfo -comment \
  -out ann.avinput
perl $ANNOVAR/table_annovar.pl ann.avinput $HUMANDB -buildver hg38 -out anno \
  -remove -protocol refGeneWithVer,cytoBand,gnomad41_exome,dbnsfp54a,avsnp151 \
  -operation g,r,f,f,f

# 9) somatic TMB (PoN added)
gatk Mutect2 -R $REFERENCE -I OC.BQSR.bam -I blood.BQSR.bam --normal-sample blood \
  -L $TARGETS --germline-resource $HIGHCONF --pon $PON \
  -O somatic.raw.vcf.gz --native-pair-hmm-threads 8 --f1r2-tar-gz f1r2.tar.gz
gatk LearnReadOrientationModel -I f1r2.tar.gz -O ob.tar.gz
gatk FilterMutectCalls -R $REFERENCE -V somatic.raw.vcf.gz --stats somatic.raw.vcf.gz.stats \
  --ob-priors ob.tar.gz -O somatic.filtered.vcf.gz
TMB=$(bcftools view -H -f PASS somatic.filtered.vcf.gz | wc -l)
echo "TMB = $(awk -v n=$TMB -v mb=60.5079 'BEGIN{printf "%.2f", n/mb}') mut/Mb"

# 10) CNVkit (hybrid mode)
cnvkit.py batch OC.BQSR.bam -n blood.BQSR.bam -f $REFERENCE -t $TARGETS \
  -g $ACCESS -m hybrid --output-reference reference.cnn -d cnvkit_out \
  --scatter --diagram -p 8
```

*Generated from the two-run reproduction described in §13; all quoted numbers are
measured in those runs. Chinese annotations target the course's bilingual audience.*
