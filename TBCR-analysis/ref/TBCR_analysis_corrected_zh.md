# TCR/BCR 免疫组库下游分析 · 矫正增强版

**第四节课 · TBCR analysis（TRUST3 输出）**
教材来源：*TBCR analysis — Downstream analysis on TRUST output files*，徐天、曾泽贤，北京大学生命科学中心 / 定量生物学中心 / 前沿交叉学科研究院。

---

## 0. 关于本版

本文件是第四节课教材的矫正增强版，沿用本仓库 WES、RNA-seq、RNA-seq-application 三个模块的既定流程：

1. 教材先在登录节点**按原样执行**，不做优化、不在运行中修错；
2. 所有遇到或预见的问题**只记录不修改**（问题清单见 §1，具体解释放在对应 Part）；
3. 随后**另跑一份修正版**到 `results/05_corrected/`，让每条修正的影响可被测量而不是被断言（本模块用户明确要求补跑修正版）；
4. 文档、脚本、日志、图一并入库。

与前三节课不同，本模块**不是**测序流程：输入是 TRUST3（从 bulk RNA-seq 读段重建 TCR/BCR 序列的组装器）的产物，全程是 R 下游分析（`ggplot2`、`dplyr`、`ggpubr`）。没有 FASTQ、没有比对、不涉及调度器——整个模块在登录节点约一分钟跑完。

本版相对原教材增加的内容：

- 每一步为什么存在、底层文件到底装了什么；
- 本次运行的准确命令与路径，以及实测数字；
- 九条已验证的坑（`TBCR-01` 到 `TBCR-09`），每条给出修正代码；
- 375 样本 TCGA-OV 队列上"原样 vs 修正"的完整对拍；
- 教材只提了一句的扩展方向（TCR beta 链、临床整合）。

---

## 1. 问题清单速览

| 编号 | 位置 | 问题 | 修正 |
|---|---|---|---|
| TBCR-01 | Part 2 | `*_report.tsv` 首行是 `#count` 表头。教材用 `read.table()` 且不设 `header = TRUE`，能跑通只是因为 R 默认 `comment.char = "#"` 静默丢弃该行。任何地方把 `comment.char` 设为 `""`，表头就会变成一行数据。 | 保持默认，或先去掉首字段的 `#` 再用 `header = TRUE`。原样运行不改，仅记录。 |
| TBCR-02 | Part 4 | `diversity <- nrow(IGH_df)` 数的是**报告行数**，不是唯一 CDR3 数。TRUST 给同一 CDR3 分配不同 V 等位或常量区时会拆成多行：演示样本 354 条 IGH 行、只有 345 条唯一 CDR3nt（该样本内部 19 条重复 CDR3nt 行）。raw diversity 因此虚高约 2-3%。 | 先按 CDR3nt 折叠（`count` 求和）再计数。 |
| TBCR-03 | Part 3 | `IGH_df$CDR3aa != "partial"` 漏掉了第二个非序列标签 **`out_of_frame`**。演示样本留下 2 条，核苷酸长度 46 与 71——不是 3 的倍数——直接进入核苷酸长度直方图。375 个 OV 样本里这类行有 9,965 条（占 IGH 行 0.92%）。 | 用 20 种标准氨基酸的白名单过滤，而不是黑名单掉一个标签。 |
| TBCR-04 | Part 3、Part 4 | `CDR3aa` 里还有**插补**序列，用 `?` 和 `_` 作占位符（对应 `CDR3_score = 0.01`）。演示样本 8 条，OV 队列 31,914 条（2.94%）。它们能通过 `!= "partial"` 过滤，被当作正常氨基酸序列统计。 | 两级过滤：L1 去 `partial` + `out_of_frame`；L2 再去含 `?`/`_` 的插补序列。 |
| TBCR-05 | Part 4 | `clonality <- 1 - (entropy / log2(diversity))` 在样本只有 1 个克隆时除零：`log2(1) = 0`、`entropy = 0`，结果为 `NaN`。OV 队列中 1 个样本（`TCGA-23-1113`）raw diversity 为 1，得到 `NaN`；修正后 9 个样本有效克隆数不足 2、3 个完全没有。教材把这些 `NaN` 直接写进 CSV。 | 显式返回 `NA`，并统计、报告有多少样本未定义。 |
| TBCR-06 | Part 4 | 正文定义 `H(X) = -sum p_i log p_i`（自然对数，nat），代码用 `log2`（bit）。熵值比公式暗示的小约 1.44 倍。clonality 不受影响，因为分子分母同底。 | 两种都算；报告用的是哪个底。修正版同时存 nat 与 bit。 |
| TBCR-07 | Part 2 | `*` 占位符（D 或 J 未检出）被当成一个基因。演示样本 IGH 行里 `*` 在 V 出现 60 次、D 169 次、J 72 次、C 132 次；IGHJ 使用率表里 `*` 排第四。使用频率被一个非基因类别稀释。 | 从分母里剔除 `*`，或单独报告为"未分配"。 |
| TBCR-08 | Part 4 | 队列循环**不过滤** `partial` 与 `out_of_frame`，而 Part 3 过滤。OV 队列 44.5% 的 IGH 行是 `partial`（40.7%）、`out_of_frame`（0.92%）或插补（2.94%），等于 Part 3 与 Part 4 描述的是两个不同的序列总体。 | 两部分使用同一套有记录的过滤规则（L1 或 L2）。 |
| TBCR-09 | 多处 | 代码卫生：`aes(x = bcr_metrics_df$diversity)` 在 `aes()` 里引用全局环境；输出文件名拼错成 `diveristy_cor_entropy.png`；两个画图块写 `p` 而不 `print(p)`（顶层与 `Rscript` 可用，`source()` 下静默不出图）；`bcr_metrics_df <- data.frame()` 在循环里靠 `rbind` 增长；Part 0 的 `ssh`/`.bashrc` 讲的是另一台机器。 | `aes()` 用裸列名；修正版改正拼写；调用 `print()`；预分配；跳过登录说明。 |

严重性排序：**TBCR-02、03、04、08** 改变数字；**TBCR-05** 产生非法值；**TBCR-06** 只改单位；**TBCR-01、07、09** 是健壮性与规范问题。

---

## 2. 环境与数据

### 2.1 环境

教材 Part 0 是在同学自己的服务器上通过 SSH 配置 conda 并写 `.bashrc`。本机账号与环境均已就绪，Part 0 缩为一行：

```bash
source /lustre1/share/miniconda3/etc/profile.d/conda.sh
conda activate airr
```

本次运行实测版本：

| 组件 | 版本 |
|---|---|
| R | 4.3.0 (2023-04-21) |
| ggplot2 | 3.4.4 |
| dplyr | 1.1.4 |
| ggpubr | 0.6.0 |

加载包时会打印 `package 'x' was built under R version 4.3.3`，因为二进制是在略新的 R 下编译的。无害，每份日志都有。

### 2.2 输入数据

全部位于 `/lustre1/share/airr/data/`（课程只读共享）。

`TRUST_output/` 是**一个**样本 `TCGA-02-0055-01A-01R-1849-01` 的十个 TRUST3 产物（TCGA RNA-seq）：

| 文件 | 大小 | 内容 |
|---|---|---|
| `*_cdr3.out` | 115 KB | 每个 consensus assembly 一行：CDR1/2/3、V/D/J/C、CDR3 score、读段支持数 |
| `*_report.tsv` | 88 KB | 每个报告 CDR3 一行：count、frequency、V/D/J/C（兼容 VDJTools） |
| `*_final.out` | 3.3 MB | 最终 contig 及核苷酸权重 |
| `*_raw.out` | 2.5 MB | 原始 contig 及核苷酸权重 |
| `*_assembled_reads.fa` | 7.6 MB | 各 assembly 的读段（FASTA） |
| `*_full.fa` | 47 KB | 完整 TCR/BCR 序列（FASTA） |
| `*_anno.fa` | 773 KB | consensus assembly 注释（FASTA） |
| `*_cluster.out` | 86 KB | consensus assembly 的聚类归属 |
| `*_evo.out` | 1.1 KB | 最终序列的进化信息 |
| `*_simplerep.out` | 108 KB | 简化报告 |

本模块只用其中两个：Part 1 用 `cdr3.out`（assembly 层），Part 2-4 用 `report.tsv`（CDR3 层）。

`OV_results/` 是 TRUST3 从 TCGA-OV 队列抽出的 375 个 `*_report.tsv`。每个文件恰好 8 个制表符分隔字段，运行前已逐个核验。

### 2.3 列结构

`cdr3.out` **没有表头**。教材手工赋 12 个列名：

```
consensus_id, index_within_consensus, V_gene, D_gene, J_gene, C_gene,
CDR1, CDR2, CDR3, CDR3_score, read_fragment_count, CDR3_germline_similarity
```

`report.tsv` **有**表头，且以 `#` 开头：

```
#count  frequency  CDR3nt  CDR3aa  V  D  J  C
```

因为 `read.table()` 默认 `comment.char = "#"`，该行被跳过，余下 930 行是数据（TBCR-01）。教材随后用 `colnames()` 赋 8 个名字。若换成 `read.delim(..., comment.char = "")`，就会多出一行 `count = "#count"` 的伪数据。

演示 report 的两个自检：

- `sum(count) = 6341` 条读段；
- `sum(frequency) = 2`，因为 frequency 列是 **TCR 与 BCR 各自归一化**（教材明确说明）——1 + 1。这是个有用的不变量：将来 TRUST 版本若改成全表归一化，总和会是 1，按链解读频率的假设就得跟着改。

演示样本的 `CDR3_score` 有五个取值，不是教材描述的三个：

| 取值 | 教材含义 | 演示样本行数 |
|---|---|---|
| 0 | partial CDR3 | 294 |
| 0.5 | （未描述） | 1 |
| 0.67 | （未描述） | 30 |
| 0.83 | （未描述） | 137 |
| 1.00 | motif 信号最强 | 479 |

### 2.4 本仓库的产出布局

```
results/
├── 01_clonotype/            Part 1（1 张 PNG、脚本、日志）
├── 02_gene_usage/           Part 2（3 张 PNG、脚本、日志）
├── 03_cdr3_length/          Part 3（2 张 PNG、脚本、日志）
├── 04_diversity/            Part 4（3 张 PNG + OV_bcr_diversity.csv、脚本、日志）
├── 05_corrected/            修正版（9 张 PNG + 2 个 CSV + 汇总、脚本、日志）
└── ref/                     本文档与教材原文 HTML
```

各脚本 `setwd()` 进自己的步骤目录，所以教材里所有相对路径（`png("...")`、`write.csv(..., "OV_bcr_diversity.csv")`）仍然写到教材预期的位置——只是落在仓库内而不是 `~/airr_output`。

---

## 3. Part 1 —— 获取克隆型信息

### 3.1 教材做了什么

读 `*_cdr3.out`，加列名，按 V/D/J/C 任一字段的 `^TR` / `^IG` 前缀拆成 TCR 与 BCR，然后统计每个 BCR consensus assembly 出现次数并画频率直方图。

### 3.2 这一步为什么存在

TRUST3 对每个表达的受体重建一个 **consensus assembly**：若干读段在 V-D-J 重排上一致就合并成一条 contig 并赋予标识（`assemble0`、`assemble1`……）。因此 `cdr3.out` 是 *assembly 层* 的表：一个 assembly 对应多行。统计一个 assembly 出现多少次（教材表里的 `Freq`）是对克隆丰度的第一层估计，直方图回答一个基本问题：这个样本是被少数扩增克隆主导，还是长尾单例为主？

必须把 TCR 与 BCR 拆开，因为 TRUST3 从同一个 RNA-seq 文库同时调两类受体，而两者的基因命名体系永不混用。

### 3.3 命令

```r
consensus_df <- read.table("/lustre1/share/airr/data/TRUST_output/TCGA-02-0055-01A-01R-1849-01_cdr3.out", sep='\t')
colnames(consensus_df) <- c("consensus_id","index_within_consensus","V_gene","D_gene",
                            "J_gene","C_gene","CDR1","CDR2","CDR3","CDR3_score",
                            "read_fragment_count","CDR3_germline_similarity")
TCR_df <- consensus_df %>% filter(grepl("^TR", V_gene) | grepl("^TR", D_gene) |
                                  grepl("^TR", J_gene) | grepl("^TR", C_gene))
BCR_df <- consensus_df %>% filter(grepl("^IG", V_gene) | grepl("^IG", D_gene) |
                                  grepl("^IG", J_gene) | grepl("^IG", C_gene))
BCR_consensus_count <- table(factor(BCR_df$consensus_id, levels = unique(BCR_df$consensus_id)))
```

### 3.4 实测产出

| 量 | 值 |
|---|---|
| `cdr3.out` 行数 | 974 |
| 唯一 `consensus_id` | 778 |
| BCR 行数 / 唯一 BCR assembly | 941 / 745 |
| TCR 行数 / 唯一 TCR assembly | 33 / 33 |
| 既非 TCR 也非 BCR 的行 | 0 |

BCR assembly 出现次数分布（895 行落在出现 1-5 次的 assembly 里，符合长尾组库）：

| 出现次数 | 1 | 2 | 3 | 4 | 5 | 7 | 10 | 11 | 16 | 18 |
|---|---|---|---|---|---|---|---|---|---|---|
| assembly 数 | 659 | 50 | 15 | 7 | 8 | 2 | 1 | 1 | 1 | 1 |

`factor(..., levels = unique(...))` 这个构造是必要的：对普通字符向量调用 `table()` 会按字母序排序，并按观测顺序固定 levels 后，`BCR_consensus_df$Var1` 才能与 assembly 编号稳定对应。

![图 1 — `BCR_consensus_frequency_distribution.png`：demo 样本里每个 BCR consensus 被看到几次。](./figures/01_clonotype.png){width=70%}

形状本身就是这张图的全部意义：745 个 BCR assembly 里有 659 个只出现一次，尾部是少数被看到 7-18 次的 assembly。那条长尾就是克隆扩增信号——一个没有扩增的库几乎全部落在第一根柱子里。当作克隆性指标用之前有两个告诫：`Freq` 数的是**行**（report 行）而不是读段，所以一个测得很深的 assembly 与一个只见到一次的 assembly 在这张图上完全一样；x 轴是 assembly 编号的因子，柱序除"观测顺序"外没有含义。

### 3.5 坑

- **`table()` 返回因子列。** `as.data.frame(table(...))` 给出 `Var1`（assembly id，因子）与 `Freq`。任何下游排序都必须先 `as.character()`。
- **`Freq` 是行数，不是读段数。** 一个 assembly 的读段支持是 `read_fragment_count`，是另一个量；教材从未使用它。若想在这个文件上做按丰度加权的 clonality，应求和 `read_fragment_count` 而不是数行。
- **TCR 表算了却没用。** `TCR_df` 只打印，后文再未使用。演示样本只有 33 条 TCR assembly 行对 941 条 BCR 行，这在实体瘤里很典型：B 细胞浸润为主，未做 TCR 富集时 TCR 克隆性很弱。

### 3.6 参数调优

- 先按 `read_fragment_count >= 2` 过滤低支持 assembly 再聚类，可压掉组装噪声。
- 对 CDR3 做近似合并（Levenshtein 距离 1，或同 V-J 且 CDR3 同一性 90%）以合并克隆家族——教材只做到按精确标识聚类。

### 3.7 扩展

- 接上 `cluster.out` 把 assembly 折成克隆簇，在簇层面重算直方图。
- 接上 `evo.out`（本样本仅 1.1 KB）观察 BCR 的体细胞高频突变谱系。

---

## 4. Part 2 —— 获取基因使用信息

### 4.1 教材做了什么

读 `*_report.tsv`，命名 8 列，抽 TCR alpha 链（`^TRAV`/`^TRAJ`/`^TRAC`）与 BCR 重链（`^IGHV`/`^IGHD`/`^IGHJ`/`^IGHC`），再算 IGHV 使用率、IGHJ 使用率与 V-J 联合分布，配两张柱状图与一张热图。

### 4.2 这一步为什么存在

V、J 基因使用是组库最粗的描述，也是跨样本、跨研究可比性最好的描述：它由重组信号偏好、选择压力以及——在肿瘤中——抗原驱动的扩增共同决定。V-J 配对再增加一个组合维度：少数 V-J 组合常常承载大部分克隆质量。

教材刻意只看**重链**（以及原则上 beta 链），因为这两条链决定抗原结合；轻链与 alpha 链是配对伙伴，若一并计入等于把同一克隆数两遍。

### 4.3 命令

```r
air_df <- read.table(".../TCGA-02-0055-01A-01R-1849-01_report.tsv", sep='\t')
colnames(air_df) <- c("count","frequency","CDR3nt","CDR3aa","V","D","J","C")
IGH_df <- air_df %>% filter(grepl("^IGHV", V) | grepl("^IGHD", D) |
                            grepl("^IGHJ", J) | grepl("^IGHC", C))
usage_IGHV_df <- IGH_df %>% group_by(V) %>% summarise(total_count = sum(count)) %>%
  mutate(usage_frequency = total_count / sum(total_count)) %>% arrange(desc(usage_frequency))
```

### 4.4 实测产出

| 量 | 值 |
|---|---|
| report 数据行 | 930 |
| `sum(count)` 总读段 | 6341 |
| `sum(frequency)` | 2 |
| TRA 行 | 6 |
| TRB 行 | 19 |
| IGH 行 | 354 |
| IGH 读段 | 2115 |
| IGHV 条目数（含 `*`） | 80 |
| IGHJ 条目数（含 `*`） | 12 |
| 不同 V-J 配对 | 173 |

按读段数排的 IGHV 前三：

| V 基因 | 读段 | 使用率 |
|---|---|---|
| IGHV3-15\*01 | 330 | 0.1560 |
| IGHV5-51\*03 | 243 | 0.1149 |
| IGHV4-31\*03 | 167 | 0.0790 |

IGHJ 前三：

| J 基因 | 读段 | 使用率 |
|---|---|---|
| IGHJ5\*02 | 751 | 0.3551 |
| IGHJ4\*02 | 492 | 0.2326 |
| IGHJ6\*02 | 478 | 0.2260 |

最高 V-J 配对：`IGHV3-15*01` + `IGHJ5*02`，288 条读段，占 IGH 质量 0.1362。

![图 2 — `top10_IGHV_genes.png`：按读段数排的前十个 IGHV 等位。](./figures/top10_IGHV.png){width=70%}

使用频率柱状图是看"这个样本是否偏斜"最快的办法：最高那条 `IGHV3-15*01` 占 IGH 读段的 15.6%，前十名合起来主导了整个分布。注意要按**等位**而非**基因**来读——`IGHV3-15*01` 与假想的 `IGHV3-15*04` 是两根不同的柱子，这让表格变稀疏，也让"前十"由等位而非基因驱动。

![图 3 — `top10_IGHJ_genes.png`：J 段的对应版本。](./figures/top10_IGHJ.png){width=70%}

J 的使用比 V 集中得多：`IGHJ5*02` 一条就占 IGH 读段的 35.5%，前三名合起来 81%。这种不对称是预期的——这里 V 有约 80 个等位在用、J 只有约 12 个，每个 J 自然分到更多质量——也正是为什么 V-J 热图有信息量、而只画 J 的图几乎就是六个数字的排名。还要注意 `*`（未分配 J）本应排第四：它没出现在这张图里，只是因为教材的 top-10 截断恰好把它排除在外。

![图 4 — `heatmap_VJ.png`：V-J 配对计数，J 在行、V 在列，颜色 = `total_count`。](./figures/heatmap_VJ.png){width=100%}

这是配对视角，同时也是本模块"技术上正确、实际上读不动"的最佳例子。那些深色格子是真实可解释的——`IGHJ5*02` 与少数几个 V 等位、`IGHJ6*02` 与 `IGHV1-69` 类等位——但 80 个等位级的 V 标签斜着 45° 叠成一条灰带，整张网格大部分是空的（173 个 V-J 对填不满 12×80 个格子），而 `*` 那一行一列就摆在坐标轴上，提醒你未分配的片段正被当成基因处理（TBCR-07）。两行改动就能让它可用：把等位折叠到基因层（`sub("\\*.*", "", V)`）、并在作图前去掉 `*` 类别。这两处留给修正版/扩展，而不是在这里悄悄改掉。

### 4.5 坑

- **TBCR-01（表头行）。** 见 §2.3：`#count` 行被跳过纯粹依赖 `comment.char`。这在教材里看不出来，但只要有人传 `comment.char = ""`，或用 `readr::read_tsv()` 读（它默认不把 `#` 当注释，会把表头保留为一行），立刻变成数据污染。
- **TBCR-07（`*` 不是基因）。** `*` 占位符表示 D 或 J 段未检出。演示 IGH 行里它在 V 出现 60 次、D 169 次、J 72 次、C 132 次。教材按原始列分组，于是 `*` 成为一个类别：在 IGHJ 使用率表里排第四。所有报告的使用率因此被轻微稀释，V-J 热图也多出一行一列 `*`。剔除 `*` 后 IGHJ 前三名不变（修正版已核验）。
- **OR 过滤与链无关。** `grepl("^IGHV", V) | grepl("^IGHD", D) | ...` 只要四个字段中任一个匹配前缀就保留该行。对 TRUST 输出是安全的，因为一行四个字段总属同一条链，但这不是保证——混链行会漏进来。更严格的形式是要求 V 前缀匹配并检查一致性。
- **等位粒度。** 基因带等位后缀（`IGHV3-15*01`）。教材不合并等位，只差一个等位的两行就是两个类别。做队列比较时应折到基因层（`sub("\\*.*", "", V)`），否则表会非常稀疏。
- **`.groups` 不一致。** V-J 块传了 `.groups = 'drop'`，单基因块没传，于是 dplyr 会打印分组消息。无害，但这个小不一致暗示这些代码块是不同时期写的。

### 4.6 参数调优

- 折并等位：分组前 `mutate(V = sub("\\*.*", "", V))`。
- 计算使用率前先限定 productive 序列（`CDR3aa` 匹配氨基酸白名单），否则 `partial` 行会贡献 V 调用却没有有效连接区。
- 按 `count` 加权（教材即如此），而不是按行计：一个单例与一个 500 读段的克隆不是同等证据。

### 4.7 扩展

- 与参考组库（如 IMGT 参考基因分布）比较，用卡方或 Jensen-Shannon 距离。
- 分别对扩增克隆与单例克隆计算使用率，观察扩增是否受 V 限制。

---

## 5. Part 3 —— CDR3 长度分布

### 5.1 教材做了什么

去掉 `CDR3aa` 等于 `"partial"` 的行，加 `nt_length = nchar(CDR3nt)` 与 `aa_length = nchar(CDR3aa)`，两种单位各画一张直方图。

### 5.2 这一步为什么存在

CDR3 长度同时是质量指标与生物学读数。生物学上，连接区长度分布受强选择（过短缺少结合面，过长受空间位阻排斥），其众数有种属与链特异性。技术上，核苷酸长度必须是 3 的倍数才可能是框内连接，所以长度分布也是组装器的快速自检。

### 5.3 命令

```r
IGH_filter_df <- IGH_df[IGH_df$CDR3aa != "partial",]
IGH_filter_df <- IGH_filter_df %>% mutate(nt_length = nchar(CDR3nt), aa_length = nchar(CDR3aa))
nt_length_df <- IGH_filter_df %>% group_by(nt_length) %>% summarise(total_count = sum(count))
```

### 5.4 实测产出

| 量 | 值 |
|---|---|
| IGH 行 | 354 |
| 去掉的 `partial` 行 | 169（47.7%） |
| 保留行 | 185 |
| 保留读段 | 1823 / 2115 |
| 核苷酸长度范围 | 30-102 |
| 核苷酸长度众数 | 66 nt，477 读段 |
| 氨基酸长度范围 | 10-34 |
| 氨基酸长度众数 | 22 aa，477 读段 |
| 保留行中核苷酸长度**不是** 3 的倍数 | **2** |

47.7% 的 partial 比例不是样本缺陷：只要能组装出足够 CDR3 就会报告，而这个肿瘤里约一半 assembly 未覆盖完整连接区。做长度分析时剔除它们是正确做法——但注意这也剔掉了一半读段，所以 Part 4 必须对同一批序列做明确决定（TBCR-08）。

![图 5 — `CDR3_nt_length_distribution.png`：只去掉 `partial` 行之后的 CDR3 核苷酸长度。](./figures/CDR3_nt_length.png){width=70%}

分布接近一个功能库该有的样子——众数 66 nt（22 个密码子）、范围 30-102、众数附近是平滑的钟形。破坏它的两根柱子只要去找就能看见：**46 nt** 与 **71 nt**，正是躲过 `!= "partial"` 黑名单的两个 `out_of_frame` assembly（TBCR-03）。它们落在不是 3 的倍数的长度上，即物理上不可能是框内连接，这也是"过滤应该用氨基酸白名单而不是单标签黑名单"的直接证据。

![图 6 — `CDR3_aa_length_distribution.png`：同一批行的氨基酸版本。](./figures/CDR3_aa_length.png){width=70%}

这张面板从另一侧暴露同样两个问题：`nchar()` 把八个插补序列里的 `?`、`_` 占位符当成残基计数（TBCR-04），所以它们的长度多了一两位；而框外行相对于钟形落在"不该有"的位置上。生物学解读本身是标准的：10-34 aa、众数 22，即人 BCR 的经典 CDR3 环长。把图与过滤决定放在一起看——光看图无法判断那条尾巴是生物学还是记账错误。

### 5.5 坑

- **TBCR-03（`out_of_frame` 逃过过滤）。** 那两条核苷酸长度非 3 倍数的行正是 `out_of_frame` 行：46 nt 与 71 nt。它们不是 `"partial"`，所以教材的黑名单保留了它们，在核苷酸直方图上成为两根孤立柱。正确过滤是白名单：只保留匹配 `^[ACDEFGHIKLMNPQRSTVWY]+$` 的序列。
- **TBCR-04（插补序列）。** 8 条 IGH 行的氨基酸串里含 `?` 或 `_`（例如 `CATTMIRDYFYYGMDV_`、`CAKPGDYDSSSYF?QW`）。`nchar()` 把这些占位符当残基计数，于是它们以正常氨基酸的身份进入长度分布。是否保留是个判断问题——核苷酸序列完整但部分插补——但必须写明选择。
- **对因子调用 `nchar()`。** 若 `CDR3aa` 被读成因子（本模块不是，因为 R 4.0 起 `read.table` 默认 `stringsAsFactors = FALSE`），`nchar()` 会返回水平序号长度。移植到旧版 R 时要注意。
- **没有零行保护。** 样本若无一行存活，直方图为空，`ggplot` 报错而不是画空图。

### 5.6 参数调优

- 用 L2 过滤（白名单）再画同一分布，看两条 `out_of_frame` 与八条插补行实际影响多大：演示样本众数仍是 66 nt / 22 aa，读段从 477 降到 474，影响可忽略——但在别的样本里不一定。
- 按总读段归一化，以比较测序深度不同的样本。

### 5.7 扩展

- 把观测长度分布与随机 V-D-J 重排的理论分布（零模型）比较，检测选择。
- 按 V 基因家族分层：长 CDR3 在某些 V 基因中富集。

---

## 6. Part 4 —— TCR/BCR 多样性指数

### 6.1 教材做了什么

定义五个指数，先在演示样本上手工计算，再对 375 个 OV 文件循环，写 `OV_bcr_diversity.csv`，画三张 Spearman 相关图。

### 6.2 各指数衡量什么

设组库有 `N` 条唯一序列、读段计数 `c_i`、总数 `C = sum c_i`、频率 `p_i = c_i / C`。

| 指数 | 教材公式 | 衡量什么 |
|---|---|---|
| raw diversity | 唯一 CDR3 数 | 只看**丰富度**，忽略丰度 |
| CPK | `(唯一 CDR3 / 总读段) * 1000` | 按测序深度归一化的丰富度 |
| entropy | `H = -sum p_i log p_i` | 丰富度与均匀度合一 |
| clonality | `1 - H / log N` | 距"其 N 条序列均匀分布"有多远；0 = 完全均匀，1 = 单克隆主导 |
| Gini-Simpson | `1 - sum p_i^2` | 随机抽两条读段来自不同克隆的概率 |

教材关于"只分析重链（BCR）与 beta 链（TCR）"的说明是范围决定而非数学决定：这两条链决定结合，且若把配对的两条链都算进来，同一克隆会被数两遍。

### 6.3 命令

```r
diversity <- nrow(IGH_df)
CPK <- diversity/sum(IGH_df$count) *1000
freqs <- IGH_df$count/sum(IGH_df$count)
entropy <- -sum(freqs * log2(freqs))
clonality <- 1 - (entropy/log2(diversity))
gini <- 1 - sum(freqs * freqs)
```

以及 375 文件循环：

```r
for (file in files) {
  air_df <- read.table(paste0("/lustre1/share/airr/data/OV_results/",file), sep="\t")
  colnames(air_df) <- c("count","frequency","CDR3nt","CDR3aa","V","D","J","C")
  IGH_df <- air_df %>% filter(grepl("^IGHV", V) | grepl("^IGHD", D) |
                              grepl("^IGHJ", J) | grepl("^IGHC", C))
  diversity <- nrow(IGH_df); CPK <- diversity/sum(IGH_df$count) *1000
  freqs <- IGH_df$count/sum(IGH_df$count)
  entropy <- -sum(freqs * log2(freqs))
  clonality <- 1 - (entropy/log2(diversity))
  gini <- 1 - sum(freqs * freqs)
  bcr_metrics_row <- data.frame(diversity, CPK, entropy, clonality, gini)
  rownames(bcr_metrics_row) <- tools::file_path_sans_ext(basename(file))
  bcr_metrics_df <- rbind(bcr_metrics_df, bcr_metrics_row)
}
write.csv(bcr_metrics_df, "OV_bcr_diversity.csv")
```

### 6.4 实测产出

演示样本（原样）：

| 指数 | 值 |
|---|---|
| diversity | 354 |
| CPK | 167.3759 |
| entropy（log2） | 6.3314 |
| clonality | 0.252280 |
| Gini-Simpson | 0.956772 |

队列（375 样本，13.9 秒）：

| 指数 | 范围 | 中位数 |
|---|---|---|
| diversity | 1 - 26,978 | 918 |
| CPK | 15.567 - 1000 | 180.7188 |
| entropy（log2） | 0 - 13.2401 | — |
| clonality | 0.0195 - 0.7306 | — |
| Gini-Simpson | 0 - 0.9995 | — |
| `NaN` 值 | 1 个（clonality，`TCGA-23-1113`） | — |

Spearman 相关（原样）：

| 配对 | rho |
|---|---|
| diversity 对 entropy | 0.7957 |
| CPK 对 entropy | -0.2470 |
| CPK 对 clonality | -0.8790 |

符号模式符合预期：丰富度与熵同向；CPK 在深度大但未必多样的样本里可以很高（它是相对读段的比值，读段少的小样本反而可能达到 CPK 1000）；CPK 与 clonality 强负相关，因为单克隆主导的组库每读段对应的唯一 CDR3 很少。

![图 7 — `diveristy_cor_entropy.png`（文件名拼写错误是教材原文）：375 个 TCGA-OV 样本的多样性与 Shannon 熵。](./figures/diversity_cor_entropy.png){width=70%}

这条紧密的单调关系（原样跑 Spearman 0.7957）是一致性检查而非发现：两者都在量丰富度，一个用计数、一个用加权和，所以一个样本不可能"唯一 CDR3 很多而熵很低"。有意思的是线周围的散开——落在线下方的样本 singleton 多（丰富度高、均匀度低），那是"库由低支持度序列拼出来"的签名。

![图 8 — `CPK_cor_entropy.png`：每千读段克隆数对熵。](./figures/CPK_cor_entropy.png){width=70%}

关系弱得多（rho −0.2470），而且本来就该弱：CPK 按读段归一，熵按克隆频率归一，于是一个测得深但分布均匀的样本与一个测得浅但很克隆化的样本可以落在同一个 CPK 上。这里的弱相关不是分析失败，而是 CPK 与熵尽管都叫"多样性"、量的却不是同一件事的证据。

![图 9 — `CPK_cor_clonality.png`：CPK 对 clonality，三者中最强的关系。](./figures/CPK_cor_clonality.png){width=70%}

漏斗形状本身就是信息：低 CPK 处 clonality 从 0 一直散到 0.7，随着 CPK 变大云团塌到接近 0 的地板上，给出 rho −0.8790。还要注意红色拟合线在右端做了什么——它外推到了**零以下**，而 `1 − H/log(n)` 不可能为负；在有界指标上做线性模型用于描述趋势没问题，用于预测单个样本则错误。图里少了一个点（`TCGA-23-1113`），因为它的单克隆让 clonality 在原样代码里成了 `NaN`（TBCR-05）。

### 6.5 坑

- **TBCR-02（diversity 是 `nrow` 而非唯一 CDR3）。** TRUST 的 report 以完整 V/D/J/C 赋值为键，同一 CDR3 可出现在多行。演示样本：354 行、345 条唯一 CDR3nt、只有 173 条唯一 CDR3aa（氨基酸层把同义变异合并）。此处 `nrow()` 高估丰富度 2.6%；队列层面在同时修正序列质量过滤后，中位数高估约两倍（§7）。
- **TBCR-05（单克隆样本 `NaN`）。** `N = 1` 时 `entropy = 0`、`log2(N) = 0`，于是 `clonality = 1 - 0/0 = NaN`。原样运行有 1 个 OV 样本命中。把 `NaN` 写进 CSV 比 `NA` 更糟：许多工具会让 `NaN` 静默穿透均值与相关计算。
- **TBCR-06（对数底）。** 公式用 `log`（nat），代码用 `log2`（bit）。熵值相差 `ln 2 = 0.6931` 倍。clonality 不受影响，因为分子分母同底——这正是归一化形式可移植的原因。
- **TBCR-08（循环无序列过滤）。** Part 3 去了 `partial`，循环没有，于是熵是在"全长 + partial"混合序列上算的。OV 队列 40.7% 的 IGH 行是 `partial`、0.92% 是 `out_of_frame`、2.94% 是插补。这是原样与修正版差异的最大来源。
- **frequency 与 count 的关系。** 循环从 `count` 重算频率，忽略 `frequency` 列。这是对的（给定列是按链归一化的），但也意味着两列从未被交叉验证。
- **循环里 `rbind`。** 375 次没问题；10,000 次就有问题了。应预分配，或对 list 用 `dplyr::bind_rows`。
- **TBCR-09（`aes()` 里用 `$`）。** `aes(x = bcr_metrics_df$diversity)` 求值的是全局向量而不是数据框列。在数据框被过滤或改名之前它能跑；`aes(x = diversity)` 更短也更安全。

### 6.6 参数调优

- `N` 小时 clonality 不稳定（要除以 `log N`）。只在克隆数超过阈值的样本上报告它（修正版标出了有效克隆数不足 2 的 9 个样本）。
- 需要经过检验的实现与置信区间时，用 `vegan::diversity` 或 `immunarch`。
- 跨样本比较 CPK 前先把读段降采样到共同深度；CPK 按深度归一化，但当丰富度饱和时并非深度不变。

### 6.7 扩展

- **TCR beta 链。** 教材自己的理由（分析决定结合的链）同样适用于 TRB。修正版用同一套规则算了 TRB：375 个样本中 369 个至少有一个有效 TRB 克隆，diversity 中位数 38，范围 0-529（见 `OV_trb_diversity_corrected.csv`）。
- **临床整合。** 教材最后一句提到可接入临床注释、比较组间指数。TCGA-OV 临床数据可从 GDC 获取；比较是 Wilcoxon 检验，或以多样性为连续协变量的 Cox 模型。
- **稀疏化曲线。** 每样本画丰富度对测序深度的曲线，把生物学丰富度与技术抽样分开。

---

## 7. 修正版与对拍

修正版位于 `results/05_corrected/`（`scripts/corrected_analysis.R`）。它先用完全相同的代码路径重算原样指数，再施加修正，因此下面每一处差异都能归因到一条具名修正，而不是归因到不同的实现。

### 7.1 使用的定义

| 层级 | 序列过滤 | 克隆定义 |
|---|---|---|
| raw | 不过滤（原样） | 一行 report |
| L1 | 去 `partial` 与 `out_of_frame` | 唯一 CDR3nt，跨行 `count` 求和 |
| L2 | 再去含 `?`/`_` 的插补序列 | 唯一 CDR3nt，跨行 `count` 求和 |

熵同时给 nat（对应公式）与 bit（对应代码）；clonality 与底数无关。

### 7.2 演示样本

| 指数 | raw | L1 | L2 |
|---|---|---|---|
| diversity | 354 | 177 | 169 |
| CPK | 167.3759 | 97.9524 | 94.0980 |
| entropy（nat） | 4.3886 | 3.7572 | 3.7308 |
| entropy（bit） | 6.3314 | 5.4205 | 5.3824 |
| clonality | 0.2523 | 0.2741 | 0.2727 |
| Gini-Simpson | 0.9568 | 0.9393 | 0.9385 |

演示样本 354 条 IGH 行中：169 条 `partial`、2 条 `out_of_frame`、8 条插补、175 条有效；L1 保留 183 行，折成 177 条唯一 CDR3nt。

![图 10 — `corrected_top10_IGHV_genes.png`：排除 `*` 占位符后的 IGHV 使用。](./figures/corrected_top10_IGHV.png){width=70%}

与图 2 对照：前五名完全相同（`IGHV3-15*01`、`IGHV5-51*03`、`IGHV4-31*03`、`IGHV1-18*04`、`IGHV3-23*01`），表里只少了一项（80 → 79），每个百分比换成了不含"无 V 调用"读段的分母。两张图都摆出来的意义是：TBCR-07 改的是**归一化**而不是生物学——所以它是一处便宜、但仍必须显式做出的修正。

![图 11 — `corrected_top10_IGHJ_genes.png`：同一分母下的 J 版本。](./figures/corrected_top10_IGHJ.png){width=70%}

J 侧同理：`IGHJ5*02`（35.5%）、`IGHJ4*02`（23.3%）、`IGHJ6*02`（22.6%）都不变，下面几名因为 `*` 不再占位而略动。如果某个下游结论依赖"第四常用的 J 段"，修正会改变它；如果依赖前三名，就不会。

![图 12 — `corrected_CDR3_nt_length_distribution.png`：L2 过滤下的核苷酸长度（175 行）。](./figures/corrected_CDR3_nt.png){width=70%}

不可能的两根柱子没了。保留的每个长度都是 3 的倍数，众数仍是 66 nt / 474 读段，范围 30-102 nt 不变——说明修正恰好只去掉了那两个 `out_of_frame` 行，别的没动。这就是"过滤真的做到了注释里写的事"的样子：钟形不动，不可能的尖刺消失。

![图 13 — `corrected_CDR3_aa_length_distribution.png`：L2 下的氨基酸视图。](./figures/corrected_CDR3_aa.png){width=70%}

去掉 `partial`、`out_of_frame` **以及**八个插补序列之后，分布是 10-34 aa、众数 22 aa。相对图 6 可见的变化在边缘：插补序列曾把个别计数虚增了一两位，当统计量是"长度超过 20 aa 的 CDR3 占比"而不是众数时，这会真的影响结果。

![图 14 — `raw_vs_corrected_diversity.png`：原样多样性（report 行）对修正多样性（唯一 CDR3nt，L1），双对数，红色虚线是 y = x。](./figures/raw_vs_corrected_diversity.png){width=70%}

375 个样本**全部**落在等值线下方——这是本次修正最重要的后果：原样指标是系统性偏高的计数，永远不可能无偏，因为 TRUST 在给同一 CDR3 分配不同 V 等位或恒定区时会把它拆成多行。偏差在低多样性处按比例更大（原样多样性 <10 个克隆的样本约 2 倍），到 10,000 以上压缩到约 1.2 倍——也就是说它恰好扭曲了那些被用来下"克隆性"结论的样本。贴着底边（修正值 = 1、原样 6-10）的点，就是折叠后只剩一个克隆的样本，也正是 §7.4 里 clonality 未定义的那几个。

![图 15 — `corrected_diversity_cor_entropy.png`：修正后的丰富度-熵关系。](./figures/corrected_diversity_cor_entropy.png){width=70%}

形状与图 7 相同，略紧一些，x 范围下移，因为计数轴减半。相关系数几乎不动（0.7957 → 0.7862），这正是"尺度修正"应有的样子：绝对数全变了，但样本的**排序**被保留（原样与 L1 多样性的 Spearman 0.9975）。反过来，如果某个修正改变了排序，那是关于修正本身的危险信号，而不是关于数据的。

![图 16 — `corrected_CPK_cor_entropy.png`：L1 下 CPK 对熵。](./figures/corrected_CPK_cor_entropy.png){width=70%}

三者中最弱的关系保持在 −0.2421（原 −0.2470）。把它读成"CPK 按读段归一这件事本身就让它与熵脱钩"的证据——去掉 partial 与重复行之后，两个量按不同倍数变化，却依旧同样弱地相关，说明这种弱关系是结构性的，不是计数 bug 的产物。

![图 17 — `corrected_CPK_cor_clonality.png`：L1 下 CPK 对 clonality，本组最强的关系。](./figures/corrected_CPK_cor_clonality.png){width=70%}

rho −0.8729（原 −0.8790）：相关是稳健的，云团并不相同。修正后 clonality 平均值上移（中位 0.1843 → 0.2196，即去掉重复行后每个样本看起来**更**克隆化，因为折叠把质量集中了），低 CPK 处漏斗略窄。三个没有有效克隆的样本是缺席而不是被画在 0 上——修正版返回 `NA` 而不是编造一个值（TBCR-05）。

### 7.3 队列（375 样本）

| 指数 | raw 中位数 | L1 中位数 | L1 / raw | Spearman(raw, L1) |
|---|---|---|---|---|
| diversity | 918 | 446 | 0.484 | 0.9975 |
| CPK | 180.72 | 113.35 | 0.640 | 0.9877 |
| clonality | 0.1843 | 0.2196 | 1.172 | 0.9690 |
| Gini-Simpson | 0.9820 | 0.9645 | 0.985 | 0.9672 |

范围，raw 到 L1：

| 指数 | raw | L1 |
|---|---|---|
| diversity | 1 - 26,978 | 0 - 14,979 |
| CPK | 15.567 - 1000 | 8.408 - 833.333 |
| clonality | 0.0195 - 0.7306 | 0.0290 - 0.8168 |
| Gini-Simpson | 0 - 0.9995 | 0 - 0.9991 |

解读：**绝对**数值变化很大（diversity 约减半，因为一半的行是 partial，剩下的还要折叠），但样本**排序**基本保持（Spearman 0.9975）。因此"A 样本比 B 样本更多样"这类结论仍然成立；"本队列多样性中位数为 918"这类结论不成立，因为那个数字是数 partial 行的产物。这就是采用修正定义的实际理由。

### 7.4 未定义值

| 层级 | clonality 的 `NA` | entropy / CPK / Gini 的 `NA` | 克隆数不足 2 的样本 | 无有效克隆的样本 |
|---|---|---|---|---|
| raw | 1 | 0 | 1 | 0 |
| L1 | 9 | 3 | 9 | 3 |
| L2 | 9 | 3 | 9 | 3 |

三个无有效 IGH 克隆的样本（`TCGA-13-0800`、`TCGA-23-1113`、`TCGA-24-1419`）所有指数均未定义。它们应从队列统计中显式排除，而不是静默插补。

### 7.5 相关结构

| 配对 | raw | L1 | L2 |
|---|---|---|---|
| diversity 对 entropy | 0.7957 | 0.7862 | 0.7875 |
| CPK 对 entropy | -0.2470 | -0.2421 | -0.2393 |
| CPK 对 clonality | -0.8790 | -0.8729 | -0.8722 |

相关结构对修正稳健——再次说明修正版改变的是尺度而不是生物学。

### 7.6 同一次修正里的其他内容

- **基因使用率**剔除 `*`：IGHV 表从 80 项降到 79 项，前五个 V 基因不变（`IGHV3-15*01`、`IGHV5-51*03`、`IGHV4-31*03`、`IGHV1-18*04`、`IGHV3-23*01`）；前三个 IGHJ 基因不变（`IGHJ5*02`、`IGHJ4*02`、`IGHJ6*02`）。
- **CDR3 长度**在 L2 上：175 行，核苷酸长度全部是 3 的倍数（两条 `out_of_frame` 已去掉），范围 30-102 nt，众数 66 nt / 22 aa，474 读段。
- **画图**在修正版里用裸列名、调用 `print()`、修正文件名拼写；原样的拼写错误 `diveristy_cor_entropy.png` 为保真保留。

### 7.7 修正版**没有**做什么

- 没有重跑 TRUST3，也没有质疑组装器的调用。
- 没有在使用率表里折并 V 等位（留作扩展）。
- 没有加入临床协变量，因为教材没有提供。
- 没有把 `NaN` 静默替换成数字：未定义就保持 `NA`。

---

## 8. 检查清单、产物与总结

### 8.1 检查清单

| 检查 | 预期 | 实测 |
|---|---|---|
| ck1 | `airr` 环境可激活并加载 ggplot2/dplyr/ggpubr | 是（R 4.3.0，3.4.4 / 1.1.4 / 0.6.0） |
| ck2 | 演示 `cdr3.out` 12 列、无表头 | 是，974 行 |
| ck3 | TCR 与 BCR 行恰好划分全表 | 是，33 + 941 = 974，0 剩余 |
| ck4 | 演示 `report.tsv` 8 列 | 是，跳过 `#` 表头后 930 数据行 |
| ck5 | report frequency 合计为 2（TCR、BCR 各自归一化） | 是 |
| ck6 | IGH 使用率合计为 1 | 是 |
| ck7 | CDR3 长度直方图众数符合预期 | 是，66 nt / 22 aa |
| ck8 | 375 个 OV 样本全部处理 | 是，13.9 秒 |
| ck9 | diversity CSV 每样本一行 | 是，375 行 |
| ck10 | 三张相关图产出 | 是 |
| ck11 | `NaN` 被记录而非隐藏 | 是，1 个样本；修正后 9 个未定义 |
| ck12 | 修正版能精确复现原样指数 | 是，演示 raw 行与 Part 4 六位小数一致 |
| ck13 | 修正前后排序保持 | 是，diversity 的 Spearman 0.9975 |

### 8.2 产物清单

```
01_clonotype/    BCR_consensus_frequency_distribution.png、log_part1.txt、scripts/part1_clonotype.R
02_gene_usage/   top10_IGHV_genes.png、top10_IGHJ_genes.png、heatmap_VJ.png、
                 log_part2.txt、scripts/part2_gene_usage.R
03_cdr3_length/  CDR3_nt_length_distribution.png、CDR3_aa_length_distribution.png、
                 log_part3.txt、scripts/part3_cdr3_length.R
04_diversity/    OV_bcr_diversity.csv、diveristy_cor_entropy.png、CPK_cor_entropy.png、
                 CPK_cor_clonality.png、log_part4.txt、scripts/part4_diversity.R
05_corrected/    OV_bcr_diversity_corrected.csv、OV_trb_diversity_corrected.csv、
                 corrected_summary.txt、raw_vs_corrected_diversity.png、
                 corrected_diversity_cor_entropy.png、corrected_CPK_cor_entropy.png、
                 corrected_CPK_cor_clonality.png、
                 corrected_CDR3_nt_length_distribution.png、corrected_CDR3_aa_length_distribution.png、
                 corrected_top10_IGHV_genes.png、corrected_top10_IGHJ_genes.png、
                 log_corrected.txt、scripts/corrected_analysis.R
```

### 8.3 总结

本模块是一分钟级的 R 分析，输入 11 个文件加一个 375 样本队列；它也是本课程第一个"教材自己的定义、而不是命令"成为薄弱点的模块。九条问题已核验并记录。改变数字的有三条：把 `nrow()` 当作丰富度度量而 TRUST report 并不以 CDR3 为键（TBCR-02）；黑名单只去 `partial` 却留下 `out_of_frame`（TBCR-03）；队列循环跳过了 Part 3 施加的质量过滤（TBCR-08）。修正版把多样性中位数减半（918 到 446），同时保持样本排序（Spearman 0.9975），因此相对比较仍然成立，绝对队列层面的陈述不成立。原样运行有 1 个样本得到 `NaN` clonality，修正后 3 个样本完全没有有效 IGH 克隆；全部按未定义报告而不插补。

---

## 附录 A —— 偏差日志

| # | 教材 | 本次运行 | 原因 |
|---|---|---|---|
| 1 | 交互式 `R`，Part 之间状态延续 | 每个 Part 一个 `Rscript`，状态重新推导 | 可复现、可留日志 |
| 2 | 输出到 `~/airr_output` | 通过 `setwd()` 输出到 `results/NN_*/` | 仓库布局 |
| 3 | Part 0 的 SSH + `.bashrc` + `conda init` | 只 `conda activate airr` | 已在该集群 |
| 4 | 输入绝对路径 | 不变 | 只读课程共享 |
| 5 | 无诊断打印 | 增加 `[diag]` 打印 | 为本文档留证 |
| 6 | 三个画图块用 `p` 不 `print()` | 原样运行不变，修正版加 `print()` | 原样运行保真 |
| 7 | 教材停在三张相关图 | 增加修正版、TRB 扩展、对拍 | 用户明确要求 |
| 8 | 文件名 `diveristy_cor_entropy.png` | 原样保留；修正版改正拼写 | 保真 |

## 附录 B —— 复现

```bash
source /lustre1/share/miniconda3/etc/profile.d/conda.sh
conda activate airr
export LC_ALL=C LANG=C

cd /gpfs1/home/bjx131_pkuhpc/class/TBCR-analysis/results
Rscript 01_clonotype/scripts/part1_clonotype.R     2>&1 | tee 01_clonotype/log_part1.txt
Rscript 02_gene_usage/scripts/part2_gene_usage.R   2>&1 | tee 02_gene_usage/log_part2.txt
Rscript 03_cdr3_length/scripts/part3_cdr3_length.R 2>&1 | tee 03_cdr3_length/log_part3.txt
Rscript 04_diversity/scripts/part4_diversity.R     2>&1 | tee 04_diversity/log_part4.txt
Rscript 05_corrected/scripts/corrected_analysis.R  2>&1 | tee 05_corrected/log_corrected.txt
```

总墙钟时间约 45 秒（队列循环 14 秒，修正版 24 秒）。不使用调度器；分析单线程，登录节点完全够用。

用于核验的关键数字：`cdr3.out` 974 行；745 个 BCR assembly；report 930 行；IGHV 首位 `IGHV3-15*01` 使用率 0.1560；354 条 IGH 行中 169 条 partial；演示样本 diversity 354；队列 375 行；raw diversity 中位数 918；修正后中位数 446；原样 1 个 `NaN` clonality，修正后 9 个未定义。
