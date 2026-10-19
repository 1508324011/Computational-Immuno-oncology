# TBCR 分析结果讲解

本文件导读 `TBCR-analysis/results/` 下的全部产物，并给出可以亲手复现的练习。配套的矫正增强版文档在 `ref/TBCR_analysis_corrected.md`（英文）与 `ref/TBCR_analysis_corrected_zh.md`（中文）。

## 1. 怎么读这个目录

```
results/
├── 01_clonotype/   Part 1 原样产出：BCR consensus assembly 频率直方图
├── 02_gene_usage/  Part 2 原样产出：IGHV / IGHJ top10 柱状图、V-J 热图
├── 03_cdr3_length/ Part 3 原样产出：CDR3 核苷酸 / 氨基酸长度分布
├── 04_diversity/   Part 4 原样产出：375 样本多样性 CSV + 三张相关图
└── 05_corrected/   修正版产出：同四项指标在 L1/L2 过滤下的结果 + 对拍图
```

每个步骤目录里都有 `scripts/`（教材原文代码，文件头登记偏差）与 `log_*.txt`（完整运行日志，含 `[diag]` 诊断行）。所有脚本都会 `setwd()` 到自己的目录，所以直接 `Rscript scripts/xxx.R` 即可重跑。

## 2. 逐类导读

### 2.1 `01_clonotype/`

- `BCR_consensus_frequency_distribution.png`：横轴是某个 consensus assembly 在样本里出现的次数（`Freq`），纵轴是有多少个 assembly 具有该次数。演示样本 745 个 BCR assembly 中 659 个只出现一次，最高的一个出现 18 次——长尾单例为主、少数扩增克隆。
- 读这张图要记住：`Freq` 是**行数**（报告条目数），不是读段数；若要按读段加权，用 `cdr3.out` 的 `read_fragment_count` 列。

### 2.2 `02_gene_usage/`

- `top10_IGHV_genes.png`：IGHV 使用率前 10。演示样本首位 `IGHV3-15*01`，占 IGH 读段 0.1560。
- `top10_IGHJ_genes.png`：IGHJ 使用率前 10。首位 `IGHJ5*02`，占 0.3551。
- `heatmap_VJ.png`：V（横）× J（纵）配对计数热图，颜色越深读段越多。173 个不同配对，最强的是 `IGHV3-15*01` + `IGHJ5*02`（288 读段，占 0.1362）。
- 注意热图里会出现 `*` 行/列：那是"未检出 D/J"的占位符，不是一个基因（详见矫正文档 TBCR-07）。修正版已剔除。

### 2.3 `03_cdr3_length/`

- `CDR3_nt_length_distribution.png` / `CDR3_aa_length_distribution.png`：按读段数加权的长度分布。演示样本众数 66 nt / 22 aa（各 477 读段）。
- 核苷酸长度图上会出现两根不在 3 倍数上的孤立柱（46 与 71 nt）：那是 `out_of_frame` 序列，教材的 `!= "partial"` 过滤没拦住它们（TBCR-03）。修正版的 `corrected_CDR3_nt_length_distribution.png` 里这两根柱消失。

### 2.4 `04_diversity/`

- `OV_bcr_diversity.csv`：375 个 TCGA-OV 样本 × 5 个指数（`diversity`、`CPK`、`entropy`、`clonality`、`gini`），第一列是样本名（去掉 `_report.tsv` 后缀）。`TCGA-23-1113` 的 `clonality` 是 `NA`——该样本只有 1 条 IGH 行，公式在此处 0/0。
- `diveristy_cor_entropy.png`（文件名拼写照教材保留）：diversity 对 entropy，Spearman rho 0.7957。
- `CPK_cor_entropy.png`：CPK 对 entropy，rho -0.2470。
- `CPK_cor_clonality.png`：CPK 对 clonality，rho -0.8790（预期中的强负相关）。
- 三张图里都会看到一条提示"Removed 1 rows containing non-finite values"——就是那个 `NaN` 样本，属预期。

### 2.5 `05_corrected/`

- `corrected_summary.txt`：修正版全部关键数字的纯文本汇总，写文档时直接引用它。
- `OV_bcr_diversity_corrected.csv`：逐样本同时给出 raw / L1 / L2 三套指标，外加 `rows_IGH`、`reads`、`n_partial`、`n_oof`、`n_imputed`、`n_valid` 六个审计列。想核对任何差异，从这张表开始。
- `OV_trb_diversity_corrected.csv`：TCR beta 链的同一套指标（教材只做 BCR 重链，这是扩展）。
- `raw_vs_corrected_diversity.png`：raw 对 L1 的 diversity 散点（双对数），虚线是 y = x；点整体落在虚线下方，说明 raw 系统性高估，但斜率接近 1，说明排序基本不变。
- 三张 `corrected_*_cor_*.png`：修正后的相关图（entropy 用 nat，与公式一致）。
- 两张 `corrected_CDR3_*_length_distribution.png`：L2 过滤后的长度分布（无非 3 倍数柱）。
- 两张 `corrected_top10_IGH*.png`：剔除 `*` 后的使用率前 10。

## 3. 手把手练习

**练习 1 —— 复现"表头行靠注释符跳过"这件事**

```bash
cd /gpfs1/home/bjx131_pkuhpc/class/TBCR-analysis/results
source /lustre1/share/miniconda3/etc/profile.d/conda.sh && conda activate airr
Rscript -e '
f <- "/lustre1/share/airr/data/TRUST_output/TCGA-02-0055-01A-01R-1849-01_report.tsv"
a <- read.table(f, sep="\t")                       # 默认 comment.char="#"
b <- read.table(f, sep="\t", comment.char="")      # 关掉注释
cat("默认:", nrow(a), "行；关注释:", nrow(b), "行\n")
cat("关注释后首行 count =", as.character(b$V1[1]), "\n")'
```

预期：默认 930 行，关注释 931 行，首行 `count` 是 `#count`。这就是 TBCR-01。

**练习 2 —— 亲手验证 diversity 的两种口径**

```bash
Rscript -e '
suppressMessages(library(dplyr))
a <- read.table("/lustre1/share/airr/data/TRUST_output/TCGA-02-0055-01A-01R-1849-01_report.tsv", sep="\t")
colnames(a) <- c("count","frequency","CDR3nt","CDR3aa","V","D","J","C")
ig <- a %>% filter(grepl("^IGHV",V)|grepl("^IGHD",D)|grepl("^IGHJ",J)|grepl("^IGHC",C))
cat("nrow(IGH)          =", nrow(ig), "\n")
cat("唯一 CDR3nt        =", length(unique(ig$CDR3nt)), "\n")
cat("唯一 CDR3aa        =", length(unique(ig$CDR3aa)), "\n")
cat("partial 行         =", sum(ig$CDR3aa=="partial"), "\n")
cat("out_of_frame 行    =", sum(ig$CDR3aa=="out_of_frame"), "\n")
cat("含 ?/_ 插补行      =", sum(grepl("[?_]", ig$CDR3aa)), "\n")
cat("保持 3 倍数?       =", all(nchar(ig$CDR3nt[ig$CDR3aa!="partial"]) %% 3 == 0), "\n")'
```

预期：354 / 345 / 173 / 169 / 2 / 8 / FALSE。最后一行 FALSE 正是 TBCR-03：`!= "partial"` 之后仍有非 3 倍数序列。

**练习 3 —— 从 CSV 里找出所有未定义样本**

```bash
cd /gpfs1/home/bjx131_pkuhpc/class/TBCR-analysis/results/05_corrected
Rscript -e '
m <- read.csv("OV_bcr_diversity_corrected.csv")
cat("样本数:", nrow(m), "\n")
cat("L1 clonality 为 NA:", sum(is.na(m$L1_clonality)), "\n")
cat("无有效克隆(L1 diversity==0):", sum(m$L1_diversity==0), "\n")
print(m[is.na(m$L1_clonality), c("sample","rows_IGH","n_partial","n_oof","n_imputed","n_valid","L1_diversity")])'
```

预期：375 样本、9 个 NA、3 个零有效克隆，名单与文档 §7.4 一致。

**练习 4 —— 量化 raw 与修正版的关系**

```bash
Rscript -e '
m <- read.csv("OV_bcr_diversity_corrected.csv")
cat("diversity  spearman =", round(cor(m$raw_diversity, m$L1_diversity, method="spearman"),4), "\n")
cat("L1/raw 中位比       =", round(median(m$L1_diversity/m$raw_diversity, na.rm=TRUE),4), "\n")
cat("clonality spearman  =", round(cor(m$raw_clonality, m$L1_clonality, method="spearman", use="complete.obs"),4), "\n")'
```

预期：0.9975 / 0.484 / 0.969。解读：绝对值缩了一半，排序几乎不变——这是本文档最重要的结论。

**练习 5 —— 看 `*` 占位符如何进入使用率表**

```bash
Rscript -e '
suppressMessages(library(dplyr))
a <- read.table("/lustre1/share/airr/data/TRUST_output/TCGA-02-0055-01A-01R-1849-01_report.tsv", sep="\t")
colnames(a) <- c("count","frequency","CDR3nt","CDR3aa","V","D","J","C")
ig <- a %>% filter(grepl("^IGHV",V)|grepl("^IGHD",D)|grepl("^IGHJ",J)|grepl("^IGHC",C))
u <- ig %>% group_by(J) %>% summarise(n=sum(count)) %>% mutate(f=n/sum(n)) %>% arrange(desc(f))
print(as.data.frame(u))'
```

预期：`*` 排第四（111 读段，占 0.0525；前三是 751 / 492 / 478），证明它不是基因却被当作基因计数。剔除后前三名不变。

**练习 6 —— 复核训练集的长度分布众数**

```bash
Rscript -e '
m <- read.table("/lustre1/share/airr/data/TRUST_output/TCGA-02-0055-01A-01R-1849-01_report.tsv", sep="\t")
colnames(m) <- c("count","frequency","CDR3nt","CDR3aa","V","D","J","C")
ig <- m[grepl("^IGHV",m$V)|grepl("^IGHD",m$D)|grepl("^IGHJ",m$J)|grepl("^IGHC",m$C),]
ig <- ig[grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", ig$CDR3aa),]
nt <- nchar(ig$CDR3nt)
cat("L2 行数:", nrow(ig), " 众数 nt:", as.integer(names(which.max(table(nt)))), "\n")
cat("全部 3 的倍数:", all(nt %% 3 == 0), "\n")'
```

预期：175 行、众数 66、全部为 3 的倍数。与修正版 `corrected_summary.txt` 的 "L2 rows: 175 / nt mode: 66" 对上。

**练习 7 —— 用 TRB 扩展数据看 TCR 侧**

```bash
cd /gpfs1/home/bjx131_pkuhpc/class/TBCR-analysis/results/05_corrected
Rscript -e '
t <- read.csv("OV_trb_diversity_corrected.csv")
cat("有有效 TRB 克隆的样本:", sum(t$diversity>0), "/", nrow(t), "\n")
cat("TRB diversity 中位数(非零):", median(t$diversity[t$diversity>0]), "\n")
cat("TRB diversity 范围:", range(t$diversity), "\n")'
```

预期：369/375、38、0-529。

## 4. 一句话总结

第四节课用约一分钟的纯 R 分析，把 TRUST3 的组装产物变成克隆型频率、基因使用、CDR3 长度和五个多样性指数；教材的薄弱点不在命令而在定义——`nrow()` 冒充唯一 CDR3、`partial` 黑名单漏掉 `out_of_frame`、队列循环跳过质量过滤——修正后多样性中位数从 918 降到 446，而样本排序保持（Spearman 0.9975），因此相对比较可信、绝对数值不可直接引用。
