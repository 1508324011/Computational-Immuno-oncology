---
title: "WES 分析结果讲解 · Results Guide"
subtitle: "两轮实跑产物的逐文件导读——怎么看、为什么是这些数、动手验证"
lang: zh-CN
toc: true
toc-depth: 3
number-sections: false
---

# 关于本指南 · About this guide

本目录收录的是**我们自己的两轮实跑分析结果**（demo 轮 + 全量轮），不含任何原始数据——
FASTQ、BAM/SAM、GVCF 一律未入库（体积与数据政策，见 §5）。每个文件都是真实管线在集群上
跑出来的产物，**未做任何人工修改**；讲解文档（本文件）教你逐类怎么读。

三份文档的分工：

| 文档 | 回答的问题 |
| --- | --- |
| 学生版路线图（`WES/ref/WES_roadmap_to_students.*`） | 命令怎么打 |
| 矫正增强版（`WES/ref/WES_roadmap_corrected*`，中英双语） | 命令哪里对哪里错、为什么、实测证据是什么 |
| **本指南 + 本目录的文件** | 证据本体——每个结果文件长什么样、怎么读、你能亲手验证什么 |

两轮运行：

| 轮次 | 数据 | 集群运行目录（大文件仍在此） | SLURM 作业 |
| --- | --- | --- | --- |
| demo | `data_new`（0.5× 降采样） | `/lustre1/user/bjx131_pkuhpc/wes_run_cluster_20260914` | 3048 |
| full | `OC_WES_all`（全量） | `/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920` | 3049–3061（8 个，依赖链） |

# 目录结构 · What is where

目录名与路线图步骤编号一致（`0_fastq` → `5_somatic_tmb`），可直接对照：

```text
WES/results/
├── demo/                      # 0.5× demo 轮（摘录）
│   ├── 0_fastq/               #   fastp 报告（OC / PBMC 各一套 html+json）
│   ├── 1_alignment/           #   flagstat · markdup 指标 · BQSR recal 表（每样本 3 件）
│   ├── 2_variant_call/        #   联合 VCF(14,782) + SNP/INDEL 两个 VQSR 输出 + 合并产物(29,564!) + recal/tranches
│   ├── 3_annotation/          #   ANNOVAR multianno(29,602 行=翻倍传播) + avinput 日志
│   ├── 4_cnv/cnvkit_out/      #   cns/call.cns + 【31 字节的空 antitarget 文件】+ scatter/diagram
│   ├── 5_somatic_tmb/         #   raw/filtered/pass VCF——全部 0 条变异（这本身就是结果）
│   ├── scripts/               #   wes_demo_steps3to10.sbatch
│   └── log/                   #   deviations.log · summary · 3047(set -u 崩溃证据) · 3048 全程日志
├── full/                      # 全量轮（摘录）
│   ├── 0_fastq/  1_alignment/ #   同 demo 构（每样本 flagstat/指标/recal）
│   ├── 2_variant_call/        #   联合 VCF(55,123) + 顺序 VQSR 矫正产物(55,123, 0 孪生) + recal/tranches
│   ├── 3_annotation/          #   avinput 日志 + invalid_input（大 multianno 未入库，见 §5）
│   ├── 4_cnv/cnvkit_out/      #   原样 amplicon 模式：cns/call.cns + 空 antitarget + scatter
│   ├── 4_cnv/cnvkit_out_hybrid/#  矫正 hybrid 模式：cns/call.cns + 2.2MB 真实 antitarget + scatter
│   ├── 5_somatic_tmb/         #   raw(2,949) → filtered → pass(108) + filteringStats
│   ├── qc/                    #   4 份 HsMetrics + demo 联合 VCF 的 bcftools stats
│   ├── scripts/               #   8 个 sbatch 脚本（全部实跑脚本）
│   └── log/                   #   deviations · 7 份 summary · 全部作业的 .out/.err
├── README.md
└── WES_results_guide.{md,html,pdf}   # 本文档
```

# 逐类讲解 · How to read each artifact

## 1. fastp 报告（`0_fastq/`，步骤 2）

**是什么：** 每样本一份自包含 HTML（浏览器直接打开，含全部图表）+ 一份机器可读 JSON。

**三处必看：**

| 看什么 | demo 实测 | 全量实测 | 判断标准 |
| --- | --- | --- | --- |
| Q30 比例（修剪后） | OC 97.0% / PBMC 97.8% | OC 97.0% / PBMC 97.8% | ≥90–95% 健康；两轮一致=同一文库抽样 |
| 接头残留 | ~0 | ~0 | 应为 0 或接近 0 |
| **插入片段峰值** | **150 bp** | **150 bp** | 应显著大于读长 137 bp——**此处偏短** |

插入峰 150 bp 是全套数据最重要的"隐性事实"：2×137 bp 的读段大量互相重叠，重叠碱基在
覆盖度统计里只算一次（全量轮实测 **25.6%** 的可用碱基耗于此，见 HsMetrics 一节）——这是
"理论 ~172× 实测只有 54.6×"的三大去向之一。

## 2. flagstat（`1_alignment/*flagstat.txt`，步骤 3.4）

**是什么：** 比对质量的 13 行汇总。三问：

| 问题 | OC（肿瘤） | PBMC（正常） | 解读 |
| --- | --- | --- | --- |
| 比对率？ | 100%（全量轮 79.6M 读段仅 311 条未比对） | 100% | 本数据集发布前已剔除未比对读段；新鲜实验室数据 96–99% 才正常 |
| properly paired？ | demo 83.29% / 全量 83.22% | 99.4–99.5% | **OC 两轮完全一致——真实肿瘤属性**：15.8% 读段对配到不同染色体（PBMC 0.43%），重度重排基因组（HGSOC）的签名 |
| secondary 比对？ | 9.1% | 0.18% | 重复序列含量差异，几个百分点属正常 |

## 3. MarkDuplicates 指标（`1_alignment/*markdup_metrics.txt`，步骤 4）

**是什么：** 每样本的 PCR/光学重复统计表。关键行是 `LIBRARY` 数据行和 `PERCENT_DUPLICATION`。

| 样本 | 检查读段对 | 重复率 | 光学重复占比 | 文库复杂度估计 |
| --- | --- | --- | --- | --- |
| demo OC / PBMC | — | **0.033% / 0.06%** | — | 抽样必然≈0（0.48%² ≈ 0.002% 的对子双双存活） |
| full OC | 36,468,418 | **7.08%** | 71.8% | 809,511,121 |
| full PBMC | 51,634,767 | **10.77%** | 57.4% | 478,822,815 |
| 老师 `analysis_full`（对照） | — | 7.39% / 11.01% | — | — |

与老师参考差 0.4 个百分点以内——这是我们的 prep 链忠实复现参考 prep 的最强单点验证。
**71.8% 是光学重复**说明重复主要来自 patterned flowcell 相邻簇，不是 PCR 过度扩增。

## 4. BQSR recal 表（`1_alignment/*.recal_data.table`，步骤 5）

**是什么：** 机器可读的经验质量模型（ReadGroupTable / QualityScoreTable / 各协变量表）。
人读时看 `RecalTable0` 的 `EstimatedQ` vs `EmpiricalQuality`：系统性差 1–5 个 Q 属正常且
值得校准。要可视化就跑 `gatk AnalyzeCovariates`。两轮的表尺寸几乎一样（108–110 KB）——
模型形状相同，只是全量轮每格计数更密；0.5× 下这些计数薄得校准基本是噪声。

## 5. 胚系变异集 VCF（`2_variant_call/`，步骤 6–7）

**是什么：** HaplotypeCaller GVCF 流程的最终产出。本目录收录的每一份都有教学点：

| 文件 | 记录数 | 为什么收录 |
| --- | --- | --- |
| `demo/OC_blood_variants.vcf` | 14,782（13,632 SNP + 1,151 indel） | 原始联合变异集；配套 `qc/demo_joint_stats.txt`（bcftools stats 全景） |
| `demo/OC_blood.snps.VQSR.vcf` | 14,782 | ApplyVQSR **SNP 模式的输出——注意它不是只有 SNP**，indel 全部透传（FILTER=`.`） |
| `demo/OC_blood.indel.VQSR.vcf` | 14,782 | INDEL 模式同理——这四份文件是理解翻倍 bug 的全套材料 |
| `demo/OC_blood.all.VQSR.vcf` | **29,564 = 2 × 14,782** | **路线图 §7.3 的 bug 本体**——每条变异两条记录 |
| `full/OC_blood_variants.vcf` | 55,123 | 全量轮原始联合集（双样本、靶区内） |
| `full/OC_blood.all.sequentialVQSR.vcf` | 55,123（PASS 53,988） | 矫正产物：顺序套用 ApplyVQSR，零未过滤孪生，过滤率 2.1% |

**翻倍 bug 亲手验证（三条命令，§6 练习 1 会用到）：**

```bash
grep -vc '^#' demo/2_variant_call/OC_blood_variants.vcf                      # 14782
grep -vc '^#' demo/2_variant_call/OC_blood.all.VQSR.vcf                     # 29564
grep -v '^#' demo/2_variant_call/OC_blood.all.VQSR.vcf | cut -f1,2 \
  | sort | uniq -c | awk '{print $1}' | sort | uniq -c                      # "14782 2" → 每位点恰好 2 条
```

**关于 "55,123 vs 老师 204,178"：** 这两个数**不可直接比**。老师的参考 VCF 是 ①单样本
（blood 的 RG `SM` 标签误写成 `OC`，CombineGVCFs 把两样本塌缩）②全基因组（无 `-L`，72%
位点在靶区外）③被 §7.3 bug 翻倍前的原始数。比 callset 前先读 VCF 头的 `##GATKCommandLine`
与 `#CHROM` 行。

## 6. VQSR recal + tranches（`2_variant_call/*.recal`、`*.tranches`，步骤 7）

**是什么：** `*.recal` 是训练出的高斯混合模型（gzip 格式）；`*.tranches` 只有 ~550 字节，
是**质量仪表盘**——每行一个敏感度档（90 / 99 / 99.9 / 100%），列 `novelTiTv` / `minVQSLod`。

三组对照（demo / 我们全量 / 老师全量）——**这是"模型可信度"最浓缩的一张表**：

| 证据 | demo（14.8k） | 全量（55.1k） | 老师（204k） | 健康的样子 |
| --- | --- | --- | --- | --- |
| truth 位点（hapmap+omni ∩ callset） | 7,388 | 28,324 | 85,159 | 随队列增长 |
| minVQSLod（100% tranche） | **−5889.9** | −1111.4 | −25.5 | −10 ~ −30 的有界尾部 |
| novel Ti/Tv（90% tranche） | 2.46（贴着 known 2.71） | 2.50（known 2.63） | 1.69（known 2.34） | novel 明显低于 known |

**怎么读：** demo 的 −5889 意味着模型在尾部毫无判别力（噪声里分离噪声）；novel Ti/Tv 贴着
known 更致命——健康 callset 的 novel 变异富集 artifacts，Ti/Tv 必须**显著低于** known（老师
那组 1.69 vs 2.34 才是标准样子）。"跑完"≠"能用"，这是 VQSR 教学的核心一课（矫正版 §7.2）。

## 7. ANNOVAR 输出（`3_annotation/`，步骤 8）

**是什么：** `multianno.txt` 每行一个（变异 × 转录本），列含 `Func` / `ExonicFunc` /
`AAChange`。demo 收录的 `OC_blood_anno.hg38_multianno.txt` 有 **29,602 行**——对 14,782
个真实变异——因为输入正是上面那份翻倍 VCF（§7.3 bug 向下游传播的实证）。全量轮的富注释
（55,165 行 × 146 列，含 gnomAD AF / REVEL / CADD / AlphaMissense）体积 27 MB 未入库，
在集群 `wes_run_full_20260920/3_annotation/`。

`invalid_input` 是 convert2annovar 拒收的行（多为复杂重排区），`.log` 里有原因——真实数据
永远有这一项，不是错误。

## 8. 体细胞 VCF（`5_somatic_tmb/`，步骤 9）

**是什么：** Mutect2 配对检出的 raw → filtered → pass 三级，加 `filteringStats.tsv`（每个
过滤器的期望 FDR/FN——比如 `orientation` 行的 FDR 0.02 就是方向偏好模型在量化 FFPE 类
artifact）。

| 轮次 | raw | filtered | PASS | TMB |
| --- | --- | --- | --- | --- |
| demo | 0 | 0 | **0** | 0.00/Mb —— 0.5× 下无统计功效，**0 是诚实的结果不是失败** |
| full | 2,949 | 2,949 | **108** | **108 / 60.5079 Mb = 1.78 mut/Mb** —— 正中 HGSOC 的 1–3 区间 |

PASS 变异的肿瘤 AF 从 0.056 到 0.48（中位 ~0.4，与高纯度肿瘤的杂合克隆突变一致）。被滤
掉的 2,841 条里，`normal_artifact`（胚系位点等位分数倾斜）和 `weak_evidence`（54.6× 可用
深度仍是灵敏度天花板）是大头。

## 9. CNVkit 输出（`4_cnv/`，步骤 10）

**是什么：** `.cnr`=每 bin log2；`.cns`=分割后；`.call.cns`=整数化拷贝数（`cn` 列 0/1/2/3，
`weight`=bin 支持度）。本目录最大的看点是**同一批 BAM、两种模式的四组对照文件**：

| 对照 | amplicon（路线图原样） | hybrid（矫正） |
| --- | --- | --- |
| antitargetcoverage.cnn | **31 字节（只剩表头）** | **2.2 MB（真实抗靶区 bin）** |
| call.cns | 290 段，臂级事件全哑火 | 357 段，高权重臂级事件 |
| scatter.png | 事件权重低 | 事件清晰 |

hybrid 版 `call.cns` 里的高权重事件（`awk '$NF>1000'` 可直接筛出）：
**chr9:78–138 Mb cn=1（9q 单拷贝缺失，weight 7,416）**、**chr8 增益 cn=3 含 8q24/MYC
（weight 7,172）**、chr5/20/21 增益、**chrY 缺失（与女性患者一致）**——一个核型自洽的
HGSOC 图景。demo 轮的 `call.cns` 则是反面教材：cn 值 324/1601/6345——参考为空时的数学
伪象，0.5× 数据 + amplicon 模式的双重噪声。

## 10. HsMetrics（`full/qc/*.hs_metrics.txt`，矫正补充步骤）

**是什么：** CollectHsMetrics 的靶区覆盖度报告——**整条原版路线图唯一缺失的指标**，
也是"深度够不够"这个检查清单问题的唯一量化答案。每份 ~60 列；最要紧的几列：

| 样本 | MEAN_TARGET_COV | ≥30X | ZERO_CVG | EXC_DUPE | EXC_OVERLAP | PCT_SELECTED |
| --- | --- | --- | --- | --- | --- | --- |
| demo OC / PBMC | 0.3× / 0.5×（中位 **0**） | 0% | 59.4% / 46.6% | ~0.03% | 26.8% / 24.9% | 79.9% / 90.8% |
| full OC | **54.6×**（中位 57） | 87.4% | 8.4% | 7.3% | 25.6% | 79.9% |
| full PBMC | **91.3×**（中位 96） | 89.6% | 8.1% | 11.0% | 23.0% | 90.7% |

**深度去向账本**（full OC）：理论 ~172× → 脱靶（~20%）→ 重复（7.3%）→ 双端重叠（25.6%）
→ **可用 54.6×**。"付了 172× 的钱，拿到 54.6× 的深度"——测序预算的实务课。
`qc/demo_joint_stats.txt` 是 demo 联合 VCF 的 bcftools stats 全景（记录数/ SNP/indel/
multiallelic/TiTv/HetHom 一应俱全），供练习用。

## 11. 日志与脚本（`log/`、`scripts/`）

- `deviations.log`：每个预期失败与既定偏差的时间戳记录——两轮运行忠实度的直接证据。
- `summary_*.log`：每步的数字汇总（本指南引用的所有数字的原始出处）。
- `*.out/.err`：SLURM 作业全程日志。**教学彩蛋**：demo/`log/wes_demo_3047.err`（124 字节）
  是 `set -u` 撞上 conda 激活的崩溃现场；full/`log/full_hsmetrics_3058.err` 记录了
  CollectHsMetrics 的参数排障（`--INTERVALS`→`--TARGET_INTERVALS`、BED 转 interval_list）。
- `scripts/*.sbatch`：**全部实跑脚本**（demo 1 份 + full 8 份）——路线图命令原文 +
  最小 SLURM 包装，复现的入口。

# 两轮对比总表 · The two rounds at a glance

| 维度 | demo（0.5×） | full（可用 54.6×/91.3×） |
| --- | --- | --- |
| 目的 | 流程机制 | 真实生物学 |
| 联合变异集 | 14,782 | **55,123**（双样本、靶区内） |
| §7.3 翻倍 bug | 29,564 = 2×14,782（本目录可亲手验证） | 110,246 = 2×55,123（30 MB VCF 未入库） |
| VQSR 可信度 | truth 7,388，尾 −5889，novel 贴 known | truth 28,324，tranches 正常，顺序矫正滤 2.1% |
| ANNOVAR | 基础协议 29,602 行（翻倍传播）；富协议挂（库名不存在） | 矫正库名 55,165 行 × 146 列（集群） |
| 体细胞 / TMB | 0（无统计功效） | **108 PASS / TMB 1.78 mut/Mb** |
| CNVkit | cn=6345 荒谬值 | hybrid 矫正：9q 缺失、chr8 含 MYC 增益、chrY 缺失 |
| 靶区深度 | 0.3×/0.5×，中位 0 | **54.6×/91.3×** |
| 墙钟 | ~25 分钟 | 最长链 3h42m（8 作业依赖串联） |

两轮的价值：**数据问题**（0、噪声、荒谬 cn）与**流程问题**（翻倍 bug、库名、模式）被彻底
分开——前者换全量数据自然消失，后者换什么数据都在，连老师的参考产物里都有。

# 未入库的大文件 · What stayed on the cluster

原始数据与超大产物不入库（本目录只放分析结果）；集群路径与去向：

| 文件 | 尺寸 | 在哪 |
| --- | --- | --- |
| 原始/clean FASTQ | 3–24 GB/文件 | `$RUN/0_fastq/` |
| 各级 BAM（sorted/markdup/BQSR） | 4–9.4 GB/文件 | `$RUN/1_alignment/` |
| demo 的 SAM（原版逐字链产物） | 165–225 MB | demo 运行目录 |
| GVCF（OC 218 MB / combined 285 MB） | — | `$RUN/2_variant_call/`、`1_alignment/` |
| full 翻倍 VCF（30 MB）、step1/分型中间 VCF | 14–30 MB | `$RUN/2_variant_call/` |
| ANNOVAR 富注释 multianno（146 列） | 27 MB | `$RUN/3_annotation/` |
| CNVkit .cnr（每 bin）、coverage .cnn、reference.cnn | 11–20 MB | `$RUN/4_cnv/` |
| diagram PDF（两模式） | 6–7 MB | `$RUN/4_cnv/cnvkit_out*/` |
| f1r2 / orientation 模型（二进制） | 2 MB | `$RUN/5_somatic_tmb/` |
| targets.interval_list / targets.bed | 6–7 MB | `$RUN/qc/`、`4_cnv/`（参考区文件，共享区有原件） |

（`$RUN` = 上表两轮目录。）这些文件都可以用 `scripts/` 里的脚本在集群上重新生成。

# 动手练习 · Exercises you can do right here

不需要集群，只需要本目录 + bcftools/samtools/grep/awk：

1. **验证 §7.3 翻倍 bug**（§4 的三条命令）：确认 29,564 = 2 × 14,782、每个位点恰好两条、
   且两条的 FILTER 不同（一条带模型裁决、一条是 `.`）。
2. **读 tranches**：打开 demo 与 full 的 `OC_blood.snps.tranches`，各找 `novelTiTv` 随档位
   （90→99→99.9→100%）的变化方向；再对比 `minVQSLod` 列的 −5889.9 vs −1111。
3. **数 Ti/Tv 和 Het/Hom**：`qc/demo_joint_stats.txt` 的 TSTV/与 PES 段——demo 的 het/hom
   ≈ 0.14（0.5× 基因型崩坏）vs 全量 1.45（健康），数字就在文件里。
4. **找最高 AF 的体细胞变异**：`bcftools query -f '%CHROM\t%POS\t[%AF]\n' full/5_somatic_tmb/somatic.pass.vcf.gz | sort -k3 -g | tail`；
   再用 `grep -v '^#' full/5_somatic_tmb/somatic.filtered.vcf.gz | grep -c PASS` 复算 TMB 的分子。
5. **CNV 对比**：对两份 `call.cns` 各跑 `awk -F'\t' '$NF>1000' ... | head`——hybrid 有臂级
   大事件、amplicon 没有；再对 demo 的 `call.cns` 找 cn>10 的荒谬段。
6. **算深度去向账**：用 fastp JSON 里的总碱基 ÷ 60.51 Mb 得理论深度（full OC ~172×），
   对照 `qc/full_OC.hs_metrics.txt` 的 54.6×，把差额拆进 PCT_EXC_* 三个口袋。

---

*本指南与矫正增强版路线图（`WES/ref/WES_roadmap_corrected*`）配套阅读效果最佳；全部数字的
原始出处是两轮运行目录的 `log/summary_*.log`（本目录已收录副本）。*
