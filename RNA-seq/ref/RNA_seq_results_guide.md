# RNA-seq 两轮实跑结果讲解（results guide）

: 配套 `RNA-seq/results/`（mini 轮 + homework/全量轮两轮入库产物）的逐类导读与动手练习。所有练习都能只在仓库内完成——不需要集群、不需要原始数据。每步原理与命令讲解见 `RNA_seq_tutorial_corrected.md`（中文版 `_zh.md`）。

## 1. 这批结果是什么

两轮全链路实跑（fastqc → fastp → STAR → multiqc → salmon 双模式 → R 差异表达/富集/热图），命令照教程原文、只加最小包装，全部偏差记录在 `*/log/deviations.log`。一轮是 mini（每样本 1 万读段、chr22 mini 索引、约 4 分钟），一轮是全量（每样本 11.8万–22.7万读对、完整 GRCh38 索引、约 31 GB 内存、13 分钟）。两轮都按教程原生方式在登录节点串行跑（理由与取证见矫正版第 0 节）。

| 关键数字 | mini | homework/全量 |
| --- | --- | --- |
| STAR unique（HBR/UHR） | 53.0–53.3% / 45.1–49.9% | 51.7–51.8% / 44.6–48.5% |
| 过滤后转录本 | 111 / 4,989 | 831 / 231,448 |
| DEG（FDR<0.05 且 \|logFC\|>1） | 45（15 上 / 30 下） | 319（171 上 / 148 下） |
| top 基因 | IGLC2（+10.3）/ MIAT（−8.5） | IGLC2 / SYNGR1 一簇 |
| 热图标签 | 基因名（mini 版代码正确） | 照原文全 "NA" → 矫正 0/50 NA |

## 2. 逐类导读

**02_fastqc/（fastqc zip ×12/轮）**：zip 内含 `fastqc_data.txt` 与完整 html，故未单独收散装 html。逐文件 `Total Sequences` 是规模核验的第一手出处——homework 轮 read1 与 read2 计数一致（118,571…227,392/文件）。

**04_fastp/（json+html 对 ×12/轮）**：json 是机读版修剪统计；`summary.before_filtering.total_reads` 与 `after_filtering.total_reads` 相除即保留率（97.4–98.7%）。

**05_STAR/（Log.final.out + SJ.out.tab 各 6/轮）**：`Log.final.out` 是比对画像的唯一权威出处——三行关键值：`Uniquely mapped reads %`、`% of reads mapped to multiple loci`、`% of reads unmapped: too short`。全量轮的"太短 47–53%"（mini 轮仅 11–13%）是**全基因组索引 + BySJout 过滤组合**的属性，不是跑错（矫正版第 5 节解读）。`SJ.out.tab` 是剪接位点表（0 列=染色体，第 6 列=注释状态 0=novel/1=annotated/2=genome-annotated）。

**06_multiqc/**：`multiqc_report.html.gz` 是汇总报告（fastqc+fastp+STAR+salmon 全在一个页面；gzip 无损入库——`gzip -cd` 即可浏览器打开，字节与集群原件一致已 cmp 验证）；`multiqc_data/` 是解析后的数据表（`multiqc_star.txt`、`multiqc_fastp.txt` 等，纯文本可直接 `column -t` 看）——不看 html 也能从这里核出本文所有数字。

**07_salmon_align/ 与 08_salmon_quasi|direct/**：每个 `*_salmon_quant` 目录里最有用的是 `aux_info/meta_info.json`——`num_mapped`/`num_processed` 就是映射率。**align 模式 100% 不是质量满分**：输入 BAM 只有 STAR 已比对片段，分母是"已比对片段"（矫正版 RNA-12）；quasi/direct 的 30–37% 才是"读段对索引"的口径。homework 轮 quasi 的 `quant.sf` 以 `.gz` 入库（231,448 行 × 4 列：Name/Length/EffectiveLength/TPM/NumReads——注意实际 5 列含 EffectiveLength），是 R 段计数的直接输入；align 模式的 quant.sf（下游不用）未入库、其目录保留其余全部文件。

**09_R/**：两轮的图。**homework 轮有一对热图证据**：`heatmap_DEG.pdf`（照教程原文跑出，50 行标签全部是 "NA"——沉默 bug 的直接产物）与 `heatmap_DEG_corrected.pdf`（一行修复后，50 行全是真实基因名）。mini 轮的 `heatmap_DEG_mini_final.pdf` 标签正常——mini 版热图代码本来就正确，作业版是被重写坏的。

**scripts/ 与 log/**：运行脚本原样（含 sbatch 头——homework 轮实际是教程原生登录节点跑，SLURM 头仅溯源保留）；`log/` 是 summary（逐步时间戳）、deviations（全部偏差与实锤发现）、R 分析日志（kept/DEG/Matched/top 基因的原始打印）。

## 3. 动手练习（全部仓库内可完成）

**练习 1——亲手复现热图 NA bug 证据**（本仓库最有价值的一对文件）：

```bash
cd RNA-seq/results/homework/09_R
pdftotext heatmap_DEG.pdf - | grep -c '^NA$'          # → 50（照教程原文）
pdftotext heatmap_DEG_corrected.pdf - | grep -c '^NA$' # → 0（一行矫正后）
```

对照矫正版第 12 节的坑链：为什么 R 一个报错都没报？

**练习 2——从 meta_info.json 核出两个口径的映射率**：

```bash
cd RNA-seq/results
for d in homework/08_salmon_quasi/*_quasi.salmon_quant; do
  python3 -c "import json; m=json.load(open('$d/aux_info/meta_info.json')); print('$d'.split('/')[-1][:7], round(100*m['num_mapped']/m['num_processed'],2),'%')"
done
# quasi：30.3–35.3%；再对 homework/07_salmon_align 同法 → 100%
```

两列数字放一起就是 RNA-12 的教训：分母不同，不能互比。

**练习 3——从 Log.final.out 重算比对画像表**：

```bash
cd RNA-seq/results
for f in homework/05_STAR/*Log.final.out; do
  echo "$(basename $f | cut -c1-7) $(grep -E 'Uniquely mapped reads % |mapped to multiple loci |too short ' $f | grep -oE '[0-9]+\.[0-9]+%' | tr '\n' ' ')"
done
```

对照矫正版第 5 节实跑表逐行核对（HBR 51.7–51.8 / UHR 44.6–48.5）。

**练习 4——fastp 保留率**：

```bash
cd RNA-seq/results/homework/04_fastp
python3 -c "
import json,glob
for f in sorted(glob.glob('*.fastp.json')):
    s=json.load(open(f))['summary']
    print(f.split('_ERCC')[0], s['before_filtering']['total_reads'],'->',s['after_filtering']['total_reads'])"
```

**练习 5——SEQC 真值自检**：看两轮 `09_R/volcano_plot.pdf` 的标签（或 `log/` 里 R 分析日志的 top 基因打印）——上调的是不是 IGLC2/IGLC3（UHR 的免疫球蛋白信号）、下调的是不是 SYNGR1/YWHAH/MIAT 一簇（HBR 脑富集）？两轮方向完全一致、top 基因同名——这就是"深度买检出数、不买不同答案"。

**练习 6——样本身份为什么只靠位置**：

```bash
cd RNA-seq/results/homework/08_salmon_quasi && ls -d *_quasi.salmon_quant
# 字母序恰好 HBR×3 在前、UHR×3 在后——教程的 condition 向量因此恰好对
# 若文件名换个前缀，设计矩阵静默错位（矫正版 RNA-01）
```

**练习 7——解开 quant.sf.gz 看计数矩阵的样子**：

```bash
cd RNA-seq/results/homework/08_salmon_quasi
gzip -cd *_quasi.salmon_quant/quant.sf.gz | head -5
gzip -cd *_quasi.salmon_quant/quant.sf.gz | wc -l   # → 231,449（表头+全转录本）
```

## 4. 未入库的大文件在哪

| 类别 | 位置 |
| --- | --- |
| 原始/trimmed FASTQ、STAR BAM（基因组+转录本坐标）、完整 salmon 目录 | 两轮集群运行目录（`/lustre1/user/bjx131_pkuhpc/rna_run_mini_20260921/`、`rna_run_homework_20261009/`） |
| STAR 索引（chr22 mini 381 MB / 全 GRCh38 ≈31 GB 内存） | `/lustre1/share/RNA_seq_class/genome/`、`/lustre1/share/RNA_seq/genome/`（教师预建） |
| salmon 索引、转录本 fasta、注释 RDS | 同上（教师预建） |
| 教师参考产物（GO_dotplot、带 ENST 行名的热图等） | `/lustre1/share/RNA_seq/` |

## 5. 一句话总结

管线产物本身都健全（两轮数字互相咬合、与教师参考在可对拍处一致）；真正值得学的东西在**接缝**里——位置身份、双模式映射率口径、版本号键错配、''' 伪注释——而它们的证据全部可以在本仓库内亲手复现。
