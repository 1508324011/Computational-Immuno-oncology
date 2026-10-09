# RNA-seq 全流程教程 · 矫正增强版（中文版）

: 文档版本 1.0 · 2026-10-09 · 以 `RNA_seq_tutorial_mini.html` 与 `RNA_seq_homework.html`（课程原件）为底本复现，辅以两轮全链路实跑（mini 轮 2026-09-21、homework/全量轮 2026-10-09），所有偏差逐条实证、记录并矫正。命令、路径、参数一律保持原文；关键数字与英文版（`RNA_seq_tutorial_corrected.md`）逐一对拍一致。

## 0. 关于本版

本版是 WES 矫正版路线图（`WES/ref/WES_roadmap_corrected.md`）的 RNA-seq 对应篇：课程原件的命令原文照录，再补上原件没有的东西——每一步的原理讲解、两轮实跑的实测数字、带根因的坑位登记表，以及**实际执行并验证过**的矫正。

**底本（课程原件，保持原样）：**

| 文件 | 角色 |
| --- | --- |
| `RNA_seq_tutorial_mini.html` | 9 步 mini 管线（每样本 1 万条读段），课堂可执行路径 |
| `RNA_seq_homework.html` | 同一管线的全量版（全深度读段 + 完整 GRCh38 索引） |
| `rna_seq.pdf` | 14 页理论课幻灯（RNA-seq 原理、STAR/Salmon 定位、R 基础） |
| `RNA_seq_code_explaination.pdf` | 9 页逐行命令讲解 |

**两轮复现**（两轮均：命令原文、最小包装、偏差全记录）：

| 轮次 | 日期 | 数据 | STAR 索引 | 运行位置 | 耗时 |
| --- | --- | --- | --- | --- | --- |
| mini | 2026-09-21 | 12 个 fastq × 10,000 reads（20 MB） | chr22 mini（381 MB） | 登录节点，串行 | ~4 分钟 |
| homework（全量） | 2026-10-09 | 12 个 fastq，11.8万–22.7万读对（112 MB） | 完整 GRCh38（约 31 GB 内存） | 登录节点，串行 | 13 分钟 |

**为什么在登录节点跑。** 教程的原生模式就是登录节点交互式执行——两份 HTML 里没有任何 sbatch/srun/SLURM 字样，与集群相关的命令只有开头一句 `ssh`。教师自己的负载测试脚本（`02_simulate_30_students_pipeline.sh`，产物在 `/lustre1/share/RNA_seq_class/loadtest/mini30_20260928_*`）模拟的就是 **30 个学生同时在登录节点跑这条管线**——这是课程设计预期登录节点运行的直接证据。mini 轮适逢全校性 SLURM 提交瘫痪（所有 sbatch 报 `Invalid account or account/partition combination`；sacct 自 9/15 起零记录；9/20 跑成功过的 WES 脚本原样重提同样被拒——集群侧重组所致，已记入该轮偏差日志）。homework 轮按明确指令在登录节点跑——尽管开跑前探测时 SLURM 已恢复（test-only 可在 cn-long 分配）。两次选择都记录在各自的 `log/deviations.log`。

**"矫正"的含义**（与 WES 版同一纪律）：照印刷版命令原样跑、把踩到/预见的坑全部记录在案，然后用真实产物执行并验证修复——修复的证据在本仓库里，不是口头描述。本文所有数字均来自两轮实跑的测量值，不是估计值。

**每步的读法：** 每个管线步骤都带同一组小节——*原文*（印刷版命令原样）、*原理*（工具与关键参数到底在做什么）、*实跑*（两轮实测数字）、*坑*（已实证，带证据）、*解读*、*调优*、*扩展*（这步熟练后往哪走）。

## 1. 问题清单速览

课程原件中已实证的坑。完整证据与矫正见各步骤小节；两条最重要的（RNA-07/09）有专门章节。

| 编号 | 位置 | 症状 / 根因 | 矫正 |
| --- | --- | --- | --- |
| RNA-01 | R sample_info | 样本身份**全程只靠位置**：`condition` 假定 `list.files()` 恰好返回 HBR×3 在前（字母序运气），文件未命名使 `colnames = NULL`，edgeR 造出 Sample1–6 | 显式命名样本；建设计矩阵前先断言顺序 |
| RNA-02 | mini 单样本 salmon | 演示单元格缺 `--gcBias --seqBias`，批量脚本却带 → 演示与批量输出存在系统性校正差异 | 全部用批量脚本的参数 |
| RNA-03 | homework 单样本 salmon | align 模式演示单元格文件名拼接错乱（UHR 样本名 + 粘贴来的 HBR `_ERCC-Mix2...chr22` 后缀）→ 找不到文件 | 批量脚本是对的；用批量 |
| RNA-04 | homework batch_salmon_quasi | `for r1 in *read1.fastq.gz` 只匹配**原始**读段——quasi 模式定量的是**未修剪**数据（mini 的 direct 模式用的是 trimmed） | 照原文跑并记录；或把 glob 指向 `*.read1.trimmed.fastq.gz` |
| RNA-05 | homework batch_STAR.sh | shebang 印成 `!/bin/bash`（丢了 `#`） | 化妆性问题；用 `bash script.sh` 跑 |
| RNA-06 | homework multiqc 单元格 | 占位路径 `/lustre1/share/RNA_seq_student/yourfold` | 指向真实工作目录 |
| RNA-07 | 合并版 "CODES" 单元格（两份文档都有） | KEGG 块被 `''' … '''` 包住——**R 里这不是注释**：整体解析成一个巨型字符串字面量（静默空操作），块内一旦出现单引号就会提前闭串、剩余代码变成**活代码** | 分步单元格用 `#` 注释是对的；合并版单元格一律当非可执行文本 |
| RNA-08 | homework 合并版 "CODES" 单元格 | 仍在过滤 `"_direct"`（mini 命名）——homework 产物叫 `_quasi` → 找到 0 个文件，下游崩 | 走分步路径（`grepl("quasi")`） |
| RNA-09 | homework 热图块 | `gene_map` 用**带版本号**的 `transcript_id` 做键，查询用的 rownames 却已剥版本号 → 全部落空 → 50 行标签全渲染成 "NA"，静默无报错 | map 键同样剥版本号（一行修复）；已验证 0/50 NA |
| RNA-10 | 教师参考产物 | 教师 `heatmap_DEG.pdf` 的行标签是**带版本号的 ENST ID**——既不是基因名也不是印刷版代码会产出的全 NA：参考跑用的是第三种代码变体 | 参考产物不能当这块的对照基准；见附录 B |
| RNA-11 | 教师 `class.R` | `salmon_files <- salmon_files[-2]` 丢掉第 2 个样本——教师自己目录的清理动作，不是方法 | 不要照抄；与参考对拍时留意 |
| RNA-12 | 解读坑 | quasi/direct 模式映射率仅 30–37%（选择性比对偏保守）；align 模式报 100%（输入 BAM 本来就是已比对片段）——两者分母不可比 | 只在同一模式内部比较映射率 |
| RNA-13 | 两份文档 | 单样本演示 `--runThreadN 1` vs 批量脚本 `4`——1 是 30 人课堂的负载决定，不是技术上限 | 独跑可放开线程（注意共享登录节点） |
| RNA-14 | R 启动 | 本集群上 R 不先 `export LC_ALL=C`/`LANG=C` 会段错误（`address (nil) ... Segmentation fault`） | 进 R 前导出两个变量（教程自带此说明，保留） |
| RNA-15 | 火山图 | `ggrepel: N unlabeled data points`——信息性提示（mini 9 条、全量 280 条），只是去拥挤 | 良性；想多标就调大 `max.overlaps` |

## 2. 环境与数据

**原文（两份文档的设置单元格）：**

```bash
ssh .....
mkdir /lustre1/share/RNA_seq_student_2026_homework/{studentId}_{name}
cd /lustre1/share/RNA_seq_student_2026_homework/{studentId}_{name}
conda activate RNA-seq
```

**原理。** 共享环境（`/lustre1/share/miniconda3/envs/RNA-seq`）带齐了管线所需全部工具，由课程钉死版本：fastp 0.23.4、FastQC 0.12.1、STAR 2.7.11b、salmon 1.10.3、multiqc、samtools、gffread，以及 R 4.4.0（edgeR/limma/tximport/ggplot2/ggrepel/pheatmap/clusterProfiler/org.Hs.eg.db/biomaRt 全装）。用共享环境而不是每人自装，是 30 人课堂能跑起来的前提，也是教师负载测试所验证的对象。

**数据是 SEQC——唯一带真值标签的 RNA-seq 数据集。** 六个双端文库：**HBR ×3**（Human Brain Reference，脑参考 RNA）与 **UHR ×3**（Universal Human Reference RNA，通用参考 RNA），各掺 ERCC 对照（HBR 用 Mix2、UHR 用 Mix1）。这里的"差异表达"有已知答案：脑富集基因在 HBR 高；UHR 是多组织混样（含免疫细胞 RNA），免疫球蛋白转录本高。本版后面所有结果核验都靠这个真值（第 10 节）。文件名里的 `chr22`/`ErccTranscripts` 标明 FASTQ 已被预筛到 chr22（+ ERCC）内容。

**教师预建资源（不要重建——原件把建索引单元格标了 `Do NOT RUN`）：**

| 资源 | 路径 | 大小 | 谁用 |
| --- | --- | --- | --- |
| chr22 mini STAR 索引 | `/lustre1/share/RNA_seq_class/genome/mini_ref_chr22/star_index` | 381 MB | mini 轮 STAR |
| chr22 mini salmon 索引 | `.../mini_ref_chr22/salmon_index` | 小 | mini 轮 direct 模式 |
| chr22 转录本 fasta | `.../mini_ref_chr22/chr22.transcripts.fa` | 小 | mini 轮 align 模式 |
| 完整 GRCh38 STAR 索引 | `/lustre1/share/RNA_seq/genome` | SA 24.9 GB + Genome 3.2 GB + SAindex 1.6 GB（加载约 31 GB 内存） | homework 轮 STAR |
| 完整转录本 fasta | `/lustre1/share/RNA_seq/genome/GRCh38_no_alt_analysis_set_gencode.v36.transcripts.fa` | 368 MB | homework align 模式 |
| 完整 salmon 索引（k31） | `/lustre1/share/RNA_seq/genome/salmon_index` | 801 MB | homework quasi 模式 |
| 基因注释 RDS | `/lustre1/share/RNA_seq_class/genome/gene_anno_gencode.v36.rds` | — | R 注释（transcript_id → gene_name） |

**数据规模（两轮实测）：**

| 轮次 | 文件 | 每样本读对 | fastp 保留 |
| --- | --- | --- | --- |
| mini | 12 × 10,000 reads | 10,000 reads（5,000 对） | —（轻度修剪） |
| homework | 12 个 fastq.gz（112 MB） | HBR 118,571–144,826；UHR 162,373–227,392 | 97.4–98.7% |

## 3. 数据获取与原始质控

**原文（homework；mini 形状一致）：**

```bash
rsync --progress /lustre1/share/RNA_seq/RNA_seq_files/*.fastq.gz .
find . -name "*.fastq.gz" -exec fastqc {} -d . -o . \;
totalreads=$(unzip -c HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2_fastqc.zip */fastqc_data.txt | grep 'Total Sequences' | cut -f 2)
echo "Total Sequences: $totalreads"
```

**原理。** 用 `rsync` 而不是 `cp`：走 Lustre 有进度条、可断点续传。逐文件 FastQC 出标准模块组（每碱基质量、接头痕迹、重复率、GC）；homework 那句 `DO NOT UNZIP` 一语双关——读段必须保持 gzipped（后面 STAR 要 `--readFilesCommand zcat`），而 `unzip -c` 是**不解压**直接从 zip 里流式读 `fastqc_data.txt`。

**实跑。** 两轮各 12 个 fastqc zip。逐文件 `Total Sequences` 确认了上面规模表（mini 每文件 10,000；homework read1 = read2 计数，每文件 118,571…227,392）。

**坑。** 这步没有新增坑（常见坑：在错误目录跑 fastqc 会静默什么都找不到；`find -exec` 是串行——这个量级无所谓）。

**调优。** FastQC 有 `-t` 线程参数；12 个文件用不上。**扩展。** multiqc（第 6 节）会汇总这些报告；`fastq_screen` 可加污染筛查。

## 4. 接头与质量修剪

**原文（批量单元格；单样本演示是同一条命令的单次形式）：**

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

**原理。** `--detect_adapter_for_pe` 靠双端读段重叠区自动推断接头对——不需要预先知道接头序列。`-l 25` 在修剪后丢弃短于 25 bp 的读段（STAR 对超短输入会退化）。每样本的 JSON/HTML 一对，正是后面 multiqc 要吃的。

**实跑（homework 轮，六个样本）：**

| 样本 | 读入 | 读出 | 保留 |
| --- | --- | --- | --- |
| HBR_Rep1 | 237,142 | 234,120 | 98.7% |
| HBR_Rep2 | 289,652 | 285,044 | 98.4% |
| HBR_Rep3 | 259,572 | 256,166 | 98.7% |
| UHR_Rep1 | 454,784 | 444,834 | 97.8% |
| UHR_Rep2 | 324,746 | 317,040 | 97.6% |
| UHR_Rep3 | 370,884 | 361,336 | 97.4% |

**解读。** 保留率 97–99% 说明文库干净：接头少、质量高——与商品化参考 RNA 一致。UHR 修剪略多于 HBR（混样成分，短片段稍多）。

**坑。** 本步无新坑；注意 trimmed 文件名（`*.read1.trimmed.fastq.gz`）是 STAR 循环匹配的目标，而 homework quasi 模式循环匹配的是**原始**名（RNA-04）——两个 glob 编码了两种不同的预处理意图。

**调优。** `--qualified_quality_phasing`/`--cut_front` 等参数存在，但 SEQC 文库都用不上。**扩展。** cutadapt/trim-galore 是等价替代；双端重叠检测若过度修剪，收紧 `--overlap_len`。

## 5. 剪接比对（STAR）

**原文（homework 批量脚本；单样本演示除 `--runThreadN 1` 外一致）：**

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

**原理（起作用的参数）。** STAR 是剪接感知的：先精确 k-mer 播种（100 bp 读段默认 50），再延伸，并对照索引内置的剪接位点库给含内含子的比对打分。`--outFilterType BySJout`：含内含子的比对，其剪接位点必须有至少一条**其他**读段（或注释 SJ 库）支持才保留——过滤垃圾内含子的主力。`--alignSJoverhangMin 8`：**新** junction 两侧各需 ≥8 bp 悬挂。`--outFilterMultimapNmax 20`：最多 20 个位点（免疫球蛋白与重复区是真的多比对）。`--outFilterMismatchNmax 999` 加 `--outFilterMismatchNoverReadLmax 0.04`：合并起来把错配上限定在读段长度的 4%。内含子 20…1,000,000 与 mates-gap 1 Mb 是 GENCODE 兼容包络。`--quantMode TranscriptomeSAM` 是**本步的承重选择**：额外产出 `<prefix>Aligned.toTranscriptome.out.bam`（坐标换算到转录本空间）——正是 salmon align 模式（第 7 节）要吃的输入。`--outSAMattributes NH HI AS NM MD` 输出多比对层级与错配标签，下游工具要。

**实跑（unique / multi / 太短，逐样本）：**

| 轮次 | 样本 | Unique % | 多位点 % | 太短 % |
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

**解读。** 三个信号，全部跨轮跨深度稳定：(1) **HBR 比 UHR 高约 5 个点**——脑参考是单组织 RNA，免疫球蛋白/重复转录本少；UHR 混样含 IG 基因，天然多比对（UHR 多位点 1.9–3.0% vs HBR 0.5–1.1%）。(2) **1 万条读段的 mini 与 10 万+ 对的全量轮 unique 率几乎相同**（HBR 53.0–53.3 vs 51.7–51.8；UHR 45.1–49.9 vs 44.6–48.5）——比对画像由读段性质决定，与深度无关。(3) homework 轮 ~47–53% 的 "unmapped: too short"（mini 同一批读段仅 11–13%，差别在索引）是**全基因组索引 + BySJout 组合**在预筛文库上的属性：最优全基因组比对低于过滤阈值的读段被记成"太短"。mini 的 chr22 索引复现的是课程自己的数字；全量索引用覆盖换严格度。记录即可不必慌：进入定量的 unique 比对部分是一致的。

**坑。** RNA-05（印刷版批量脚本 shebang 丢了 `#`——化妆性）。RNA-13（线程 1 vs 4——课堂负载设计；我们 homework 轮用 4 并加内存守卫：每次 STAR 调用前等待 `MemAvailable` ≥ 37 GB，因为加载全量索引要 ~31 GB 且登录节点共享）。mini 索引只要 ~0.4 GB——那轮不需要守卫。

**调优。** 节点空闲时上调 `--runThreadN`；`--twopassMode Basic` 提升新 junction 灵敏度（约 2 倍耗时）；后面要跑基因组坐标工具的话 `--outSAMtype BAM SortedByCoordinate` 可省一次 samtools sort。

**扩展。** samtools `sort`/`index` + `flagstat` 做 BAM QC；RSeQC（`infer_experiment`、`inner_distance`）查链特异性与插入分布；`--quantMode GeneCounts` 在转录本 BAM 之外顺手出一份原始计数矩阵。

## 6. 汇总质控（multiqc）

**原文（homework 单元格；mini 用真实工作目录）：**

```bash
multiqc /lustre1/share/RNA_seq_student/yourfold
```

**原理。** multiqc 扫一遍目录树，认出它认识的工具产物（fastqc zip、fastp JSON/HTML、STAR `Log.final.out`、salmon `quant.sf` 等），合并成一份按样本分列的报告——一眼看全 12 个 fastq 与 6 个样本的最快方式。

**实跑。** `multiqc "$PWD" -o multiqc_report`（RNA-06：印刷路径是占位符——`yourfold` 字面就是"你的目录"）。两轮都产出了汇总报告；本文引用的 STAR 与 fastp 数字都可在其中交叉核对。

**坑。** RNA-06（占位路径）。另注意：multiqc 只汇报指定路径下找到的东西——指错目录会"成功"地产出空报告而不是报错。

**调优。** `-f` 重建、`--cl` 把命令行嵌进报告、写个 `multiqc_config.yaml` 钉死教学用的模块顺序。**扩展。** `--pdf` 导出；管线中点（STAR 后）与终点各跑一次得到前后对照视图。

## 7. 定量 I：比对模式（salmon align）

**原文（homework 批量单元格——正确版；单样本单元格的拼接错乱文件名见 RNA-03）：**

```bash
TRANSCRIPTOME="/lustre1/share/RNA_seq/genome/GRCh38_no_alt_analysis_set_gencode.v36.transcripts.fa"
for bam_file in *Aligned.toTranscriptome.out.bam; do
    sample_name=$(basename "$bam_file" Aligned.toTranscriptome.out.bam)
    salmon quant -t "$TRANSCRIPTOME" --libType A -a "$bam_file" -o "${sample_name}_align.salmon_quant" --gcBias --seqBias
    echo "Processed: $sample_name"
done
```

**原理。** 这是**基于比对**的定量：读段落在哪已由 STAR 决定；salmon 的 align 模式（`-a`）吃转录本坐标 BAM，把多比对读段按 EM 做分数分配，估计每条转录本的丰度。`--gcBias --seqBias` 从比对本身估计并校正序列内容与 GC 偏好。`-t` 是 BAM 里转录本名所指的转录本 fasta；`--libType A` 让 salmon 自己推文库类型。

**实跑。** 两轮各六样本，全部 **100.0% 片段"mapped"**——这不是质量数字：输入 BAM 里只有 STAR 已比对的片段，分母是"已比对片段"而非"文库读段"（RNA-12）。

**坑。** RNA-03（单样本演示单元格的文件名是 UHR/HBR 粘贴混排——照抄就在找不到文件处报错）。RNA-02（mini 单样本单元格缺 `--gcBias --seqBias`、批量带——一个样本用演示命令、其余用批量脚本，会造成跨样本的系统性校正差异）。

**解读。** align 模式继承 STAR 的全部决定（junction 过滤、多比对上限、"太短"丢弃）。它的计数是"STAR 眼中的转录本组"；quasi 模式（第 8 节）是 salmon 自己的、更保守的视角——两种都跑的意义就在这个对照。

**调优。** `-p` 线程；怀疑位置衰减加 `--posBias`。**扩展。** RSEM 是同类 BAM 的经典替代；`--numBootstraps` 加不确定性重采样，tximport 可接。

## 8. 定量 II：选择性比对模式（salmon quasi / `--validateMappings`）

**原文（homework 批量单元格）：**

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

**原理。** quasi 模式把原始读段直接对到转录本索引（没有基因组、无需发现 junction）。`--validateMappings`（salmon ≥1.0 的默认引擎）用选择性比对验证候选映射——由比对分数阈值决定取舍而非纯 k-mer 命中，所以映射率读起来比 STAR **低**：没有好比对的读段被丢弃而不是硬塞。

**实跑（映射率，取自各目录 `aux_info/meta_info.json` 的 `num_mapped / num_processed`）：**

| 轮次 | 模式 | HBR_Rep1 | HBR_Rep2 | HBR_Rep3 | UHR_Rep1 | UHR_Rep2 | UHR_Rep3 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| mini | direct（trimmed 读段） | 36.5% | 36.5% | 36.4% | 33.8% | 31.2% | 32.2% |
| homework | quasi（原始读段，照原文） | 35.2% | 35.3% | 35.1% | 32.6% | 30.3% | 32.0% |

**坑。** RNA-04 是重点：`for r1 in *read1.fastq.gz` 只匹配**原始**文件（`*.read1.trimmed.fastq.gz` 结尾不同），所以 homework 的 quasi 模式定量的是**未修剪**读段，而 mini 的 direct 模式定量的是 trimmed——我们照原文跑了并记录（接头片段短，映射率几乎不动，但两份文档的预处理故事不一致）。另有 RNA-12：这里的 30–37% 对 align 模式的 100% 是**分母不同**，不是质量崩塌。

**解读。** ~1/3 的映射率由两个独立因素决定：读段被预筛到 chr22+ERCC 内容而索引是 231,448 条转录本的全 GENCODE v36（salmon 丢弃无处可比的读段），加上选择性比对的分数阈值。比率跨轮惊人稳定（HBR 35.1–36.5%、UHR 30.3–33.8%）——还是文库性质，两个深度都可复现。

**调优。** `-p 8` 是 homework 批量单元格印的数——在共享登录节点上这已经是实打实的线程预算。**扩展。** `--numBootstraps 100` 配 tximport 的 `countsFromAbundance` 得到带自助方差的定量；`--seqBias` 与片段长度分布有交互——QC 时看输出里的 `fldLen.png`。

## 9. 读入定量结果并过滤（tximport → edgeR）

**原文（R 单元格，分步路径）：**

```r
salmon_files <- list.files(pattern = "quant\\.sf$", full.names = TRUE, recursive = TRUE)
salmon_files <- salmon_files[grepl("quasi", salmon_files)]     # homework；mini 过滤 "_direct"
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

**原理。** `tximport` 把各样本的 `quant.sf`（Name、Length、TPM、NumReads）收进一个矩阵；`txOut = TRUE` 保持**转录本层级**的行（全量索引 231,448 行）。过滤（`≥10 counts 且 ≥2 样本`）删掉稀疏到没法建模的行——在本文库规模下那是几乎全部，这是数据集构造使然（chr22+ERCC 来源的读段对全转录本索引）。

**实跑（建设计矩阵之前该查的东西——RNA-01）：**

| 检查 | mini | homework |
| --- | --- | --- |
| 索引内转录本数 | 4,989（chr22 mini） | 231,448（全 GENCODE v36） |
| 过滤后保留 | 111 | 831 |
| `colnames(txi$counts)` | `NULL` → Sample1–6 | `NULL` → Sample1–6 |
| 救了场的那点运气 | HBR×3 在前（字母序） | HBR×3 在前（字母序） |

**坑。** RNA-01 是结构性问题，共三环：(1) `list.files()` 返回字母序，而印刷版 `condition` 向量写死 HBR 在前——**只靠文件名的字母序运气**才对；(2) 文件未命名传入，`colnames` 是 `NULL`，edgeR 造出 Sample1–6——样本身份从此纯靠位置；(3) 下游（设计矩阵、对比、热图 `sample_info` 行名）没有任何东西能拦住一次错位。一行修复：给 `salmon_files` 命名，再由名字构造 `condition`。两轮的检查输出如上——顺序两次都成立，是记录而非假设。RNA-14：进 R 前 `export LC_ALL=C; export LANG=C`（教程自带此说明）。

**解读。** 111/4,989 与 831/231,448 看似存活率不同，实为同一现象在两种索引广度下：两轮都只有 ~2% 的行过关。深度决定**多少**转录本过 10 counts（111 → 831），索引决定分母。

**调优。** 随深度上调 `counts_threshold`（10 是本教学规模；全转录本研究常用 ≥10–20 的 CPM 类门槛）。**扩展。** `txOut = FALSE` 配 `tx2gene` 表在建 edgeR 前聚合到基因层级（更常见的生产路径）；DESeq2 配 `tximport` 是即插即用的替代；`rtracklayer`/`biomaRt`（教程已加载）可重建注释。

## 10. 差异表达、注释与火山图（glmLRT + volcano）

**原文（R 单元格）：**

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

**原理。** edgeR 的广义线性模型路径：`estimateDisp` 估计负二项离散度，`glmFit` 拟合，`glmLRT` 对最后一个模型项做似然比检验（这里是 `conditionUHR`）。`regulate` 规则是合取：FDR < 0.05 **且** |logFC| > 1——统计显著与效应量一起卡，这也是火山图虚线画在 ±1 与 −log10(0.05) 的原因。注释块把 GENCODE 版本后缀剥掉（`ENST00000615943.1` → `ENST00000615943`），**两侧都剥**再做 `match()`——注意这正是热图块忘掉的那条纪律（第 12 节）。

**实跑：**

| | mini | homework（全量） |
| --- | --- | --- |
| up-regulated | 15 | 171 |
| down-regulated | 30 | 148 |
| not sig | 66 | 512 |
| Matched（gene_name） | 111 / 111 | 831 / 831 |
| top up（按 FDR） | IGLC2（logFC ≈ +10.3，FDR ≈ 1e-145）、IGLC3、XBP1、MYH9 | IGLC2、IGLC3、… |
| top down | MIAT（≈ −8.5）、YWHAH、MAPK8IP2、ACO2、SYNGR1、SEPTIN3、CBX7 | SYNGR1、SEPTIN3、YWHAH、RPL3、CBX6、… |

**解读——SEQC 真值核验。** 每个方向都与已知答案一致：**UHR 上调**：免疫球蛋白 lambda 链 IGLC2/IGLC3（UHR 混样含免疫细胞 RNA）——教科书级；**HBR 高（即 UHR 下调）**：MIAT、YWHAH、MAPK8IP2、ACO2、SYNGR1、SEPTIN3——脑富集转录本。两个深度下 top 基因是同一条转录本（IGLC2，|logFC| > 10、FDR ~1e-145）——信号在 1 万读段时已饱和；深度买到的是**更多**被检出的基因（45 → 319），不是不同答案。这正是 mini/homework 配对的核心课程点。

**坑。** RNA-15（`ggrepel: 9 / 280 unlabeled`——去拥挤提示，良性）。glmLRT 收拟合对象；多传参数会改变被检验项（mini 轮收割时踩过并当场纠正——教程原文 `glmLRT(fit)` 才是对的）。

**调优。** `topTags(n = Inf)` 再按 logFC 排序，能看见被 FDR 排序遮住的大效应基因。**扩展。** `decideTests` 跨对比做 Venn；`topTags` 默认 BH 校正就是这里的 FDR。

## 11. 富集分析（GO——以及永远跑不起来的 KEGG）

**原文（分步单元格——可执行路径）：**

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

**原理。** 富集问的是：上调基因集相对背景在注释条目上是否过度代表。`keyType = "SYMBOL"` 直接喂基因符号（免手工转 ID）；`ont = "ALL"` 合并 BP+MF+CC；`pvalueCutoff = 1` 关掉展示侧阈值，让 top-10 按真实显著性排序。

**实跑。** mini 15 个上调基因 / homework 171 个。homework top GO 条目（p.adjust 0.0075–0.0088）：*viral RNA genome replication*、*single-stranded viral RNA replication via double-stranded DNA intermediate*、*glycosyl compound catabolic process*、*protein targeting to mitochondrion*、*DNA dealkylation*、*cytidine catabolic process*、*cytidine deamination*、*cytidine-to-uridine editing*。教师参考 `GO_dotplot.pdf` 的 top 条目相同（*viral RNA genome replication*、*glycosyl compound catabolic process*，p.adjust 0.0075）——与参考交叉验证通过，不只是自洽。

**坑——RNA-07/08，都在合并版 "CODES" 单元格上。** 分步单元格用 `#` 注释掉 KEGG（正确）。合并版单元格却把整个 KEGG 块包进 `''' … '''`。R 里 `'''` **不是**注释：它解析成空串 `''` 加一个开串符 `'`，于是整块变成**一个巨型字符串字面量**，被 R 求值后自动打印——它**碰巧**表现得像个块注释，什么都不做。失效模式是潜伏的：块内一旦出现单引号就会提前闭串，剩余代码全部变成**活代码**。一律把 `'''` 块当非可执行文本。RNA-08：homework 的合并版单元格还在过滤 `"_direct"`（mini 命名）——homework 产物叫 `_quasi`，它找到零个文件、死在 tximport；分步路径（`grepl("quasi")`）才是可执行的。

**解读。** cytidine deamination / cytidine-to-uridine editing 一簇不是噪声：表达免疫球蛋白的细胞靠 AID/APOBEC 驱动的胞嘧啶脱氨做亲和成熟——GO 读出的正是火山图上同一个 IG 信号。*viral RNA genome replication* 一类条目是上调集里 IG/ERV 区基因的注释侧效应（教师参考同样显示，可对照读）。

**调优。** `enrichKEGG` 需要 ENTREZ ID（被注释的块里做了转换）且要访问 KEGG REST 网络——这就是它常被课堂略过的原因。**扩展。** `compareCluster` 跨上调/下调/并集比较；GSEA（`gseGO`）吃整个排序表，用的信息比 top 集多。

## 12. 热图——沉默的全 NA 坑

**原文（homework 热图块，逐字）：**

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

**坑链，一环一环。** 第 2 行把 rownames 剥掉版本号（`ENST00000390323.5` → `ENST00000390323`）；第 6 行建 `gene_map` 用的是**带版本号**的键（`setNames(gene_anno$gene_name, gene_anno$transcript_id)`——没剥）；第 7 行拿剥过号的 rownames 查带号的 map → **全部落空**；把全 `NA` 向量赋给 rownames 在 R 里合法；pheatmap 照样渲染——**五十行全标 "NA"，无报错无警告**。我们实跑产出的 `heatmap_DEG.pdf` 用 `pdftotext` 抽出来正好 50 个 `NA`。这与第 10 节注释块是同一类版本号错配——那块两侧都剥所以对——mini 版的热图块也是对的（一个 `clean_id()` 帮手同时作用于**数据行与注释 ID**，外加 `make.unique()` 和 `labels_row =`）。

**一行矫正（已执行并验证）：**

```r
gene_map <- setNames(gene_anno$gene_name, sub("\\..*$", "", gene_anno$transcript_id))
```

矫正后运行（`heatmap_DEG_corrected.pdf`）：**0/50 NA**，行标签是 IGLC2、IGLC3、MCM5、SERPIND1、LIF、…——真实的 top DEG 基因。

**坑。** RNA-09（上述）。RNA-10：教师参考 `heatmap_DEG.pdf` 的行标签是**带版本号的 ENST ID**——既非基因名（mini 风格）也非全 NA（印刷版 homework 风格）：参考跑用的是又一种代码变体，不能当这块的对拍基准（附录 B）。homework 块还丢了 mini 版的 `scale = "row"`——行是裸 logCPM，高表达基因会把 scale="row" 才能显示的图案在视觉上淹没。

**解读。** 热图在这条管线里的任务是给 DEG 结论做可视化断言：样本按条件聚（HBR/UHR 列分开）、top 基因分成两个方向块。这个坑只弄坏了标签——数字、聚类、注释条一直是对的，而这恰恰是它能一路存活进印刷讲义的原因。

**调优。** `scale = "row"`（逐行 z 分）、`cutree_rows` 强制块结构、两个条件之间 `gaps_col = 3`。**扩展。** 注释多时换 `ComplexHeatmap`；50 个标签挤时调字号。

## 13. 双轮对比（mini vs homework/全量）

| | mini（2026-09-21） | homework/全量（2026-10-09） |
| --- | --- | --- |
| 每样本读段 | 10,000 | 118,571–227,392 对 |
| STAR 索引 | chr22 mini，381 MB | 全 GRCh38，约 31 GB 内存 |
| STAR unique（HBR / UHR） | 53.0–53.3% / 45.1–49.9% | 51.7–51.8% / 44.6–48.5% |
| STAR 太短 | 11.5–12.9% | 47.1–52.7%（全量索引 + BySJout 严格度） |
| salmon align 模式 | 100%（已比对片段分母） | 100%（同） |
| salmon direct/quasi | 31.2–36.5%（trimmed） | 30.3–35.3%（原始读段，照原文） |
| 索引内转录本 | 4,989 | 231,448 |
| 过滤后保留（≥10 counts 且 ≥2 样本） | 111（2.2%） | 831（0.36%） |
| DEG（FDR<0.05 且 \|logFC\|>1） | 45（15 上 / 30 下） | 319（171 上 / 148 下） |
| 两个方向的 top 基因 | IGLC2（+10.3）/ MIAT（−8.5） | IGLC2 / SYNGR1 一簇 |
| 注释 Matched | 111/111 | 831/831 |
| GO top 条目 | （集小） | viral RNA genome replication 等，与教师参考一致 |
| 热图标签 | 基因名（正确块） | 照原文全 "NA" → 矫正后 0/50 NA |
| 墙钟时间 | ~4 分钟 | 13 分钟 |
| 跑在哪 | 登录节点，串行 | 登录节点，串行（STAR 带内存守卫） |

**什么在缩放、什么不在。** 墙钟 4 → 13 分钟（索引加载主导；全量索引是 mini 的 ~80 倍）。DEG 数 45 → 319（深度买到检出）。**答案**没动：映射率画像、top 基因、方向全部随深度不变——这正是先跑 mini 的意义所在。两轮之间唯一翻转的东西来自**文档**而不是数据：mini 的热图块注释正确，homework 的不正确。

## 14. 检查清单、产物清单与总结

**检查清单（附两轮实测答案）：**

| # | 检查项 | mini | 全量 |
| --- | --- | --- | --- |
| ck1 | 建矩阵前样本顺序（HBR×3 在前？） | 成立 | 成立 |
| ck2 | `colnames(txi$counts)` NULL → Sample1–6（位置身份） | 是 | 是 |
| ck3 | STAR unique 落在 44–54% 区间，HBR 高于 UHR | 45.1–53.3% | 44.6–51.8% |
| ck4 | 多位点低，UHR > HBR（IG 成分） | 0.5–2.1% | 1.1–3.0% |
| ck5 | quasi 映射 30–37%、align 100%（分母不同——别互比） | 31.2–36.5% | 30.3–35.3% |
| ck6 | 两种索引广度下保留率 ~2% 与 ~0.4% | 111/4,989 | 831/231,448 |
| ck7 | DEG 方向符合 SEQC 真值（UHR 中 IGL 上调、脑基因下调） | 是 | 是 |
| ck8 | 两个深度 top 基因一致 | IGLC2 | IGLC2 |
| ck9 | 注释块打出 `Matched: n/n` | 111/111 | 831/831 |
| ck10 | GO top 条目对教师参考 | — | 一致（viral RNA genome replication，p.adjust ~0.008） |
| ck11 | 热图标签是基因名不是 "NA" | 是（正确的块） | 矫正后 0/50 NA（照原文 50/50 NA） |
| ck12 | 进 R 前 `LC_ALL=C`/`LANG=C` 已导出 | 是 | 是 |
| ck13 | ggrepel "unlabeled" 警告记录为良性 | 9 | 280 |

**产物清单：** 两轮的分析产物都在仓库 `RNA-seq/results/` 下（只放分析产物——无原始读段、无 BAM）：fastqc/fastp 报告、各 6 份 STAR `Log.final.out`、multiqc 报告、全部 24 份定量的 `quant.sf`、R 分析日志、homework 轮四张 PDF（火山 / GO / **照原文热图，全 NA** / **矫正热图，0 NA**）加 mini 轮三张 PDF，以及两轮完整的 `summary` + `deviations` 日志与运行脚本。重型中间产物（BAM、原始/trimmed FASTQ、完整 salmon 目录）留在集群运行目录（`/lustre1/user/bjx131_pkuhpc/rna_run_mini_20260921/`、`/lustre1/user/bjx131_pkuhpc/rna_run_homework_20261009/`）。

**总结。** 管线本身是健全的；坑全在管外围手写接缝上——样本身份只靠位置与字母序运气、两个 salmon 模式的映射率分母不可比、quasi 循环悄悄定量了原始读段、热图注释重写静默渲染出 50 行无标签。每一条都在真实运行中复现，每个修复都执行并验证过。SEQC 数据集全程尽到了职责：生物学答案（IGLC2 上调、脑基因下调）在两个深度下稳定，且在参考可对拍处与教师参考一致——而如附录 B 所示，"参考可对拍"这件事本身也是需要核验的，不是想当然。

## 附录 A. 偏差日志（两轮，节选但忠实原文）

完整日志：`RNA-seq/results/mini/log/` 与 `RNA-seq/results/homework/log/`（summary + deviations + R 分析）。

**mini 轮（2026-09-21）——节选：**

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

**homework 轮（2026-10-09）——节选：**

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

## 附录 B. 教师参考产物 vs 印刷版文档（取证，路径原文）

教师自己的运行在 `/lustre1/share/RNA_seq/` 留下了产物。与印刷版文档相对，三处独立偏离：

1. **参考热图是第三种变体。** `/lustre1/share/RNA_seq/heatmap_DEG.pdf` 的行标签是带版本号的 ENST ID（如 `ENST00000407418.8`）。印刷版 homework 块产出的是全 "NA"（RNA-09，已复现）；mini 块产出基因名。参考两者都不是——它那份代码保留了带版本号的 ENST 行名、从未改名。
2. **`class.R` 丢了一个样本。** `salmon_files <- salmon_files[-2]`——tximport 前删掉第 2 个定量目录。与教师自己运行区里重复目录的清理一致（共享目录里同时有 `UHR_Rep3...chr22Aligned...bam` 和双名 `UHR_Rep3...chr22_ERCC-Mix2...chr22Aligned...bam`），而与印刷文档里的任何方法都对不上。
3. **拼接错乱的单样本 salmon 单元格真被跑过。** 上面的双名 BAM 正是 homework 单样本 align 单元格（RNA-03，UHR 名 + 粘贴的 HBR `_ERCC-Mix2...chr22` 后缀）会产出的文件名——参考区里留着这个 bug 自己的输出文件名。

参考**可**对拍的地方与我们的结果一致：`GO_dotplot.pdf` 的 top 条目与 p.adjust 和我们吻合（附录 A 的 RESULT）。WES 版的教训原文重演：参考产物是**某人实际命令**的证据，不是印刷版命令正确的证明——只在代码路径已知处与它对拍。

## 附录 C. 复现

两轮的运行脚本与 R 脚本原样存放在 `RNA-seq/results/mini/scripts/` 与 `RNA-seq/results/homework/scripts/`，旁边是完整日志。homework 运行脚本在自己的头注释里写明了它的包装偏差（SLURM 头仅为溯源保留——该轮按教程原生方式在登录节点跑）。复跑 mini 轮约 4 分钟、只需共享环境；复跑 homework 轮串行约 13 分钟、每次 STAR 调用需 ~31 GB 空闲内存。
