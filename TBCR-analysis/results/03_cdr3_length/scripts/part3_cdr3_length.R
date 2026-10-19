# part3_cdr3_length.R — 教材 Part 3 原文（as-written）
# =============================================================================
# 教材：TBCR analysis（Tian Xu / Zexian Zeng, PKU）
# 偏差登记：
#   1. 交互式 R → Rscript
#   2. 教材在 ~/airr_output 下跑，本脚本 setwd 到本步骤目录
#   3. 教材 Part 3 延续 Part 2 的 R 会话（IGH_df 已存在）；脚本形式在此
#      重新读表并重建 IGH_df（与 Part 2 完全相同的过滤表达式）
#   4. LC_ALL=C；加诊断打印（原文 print 保留）
# 注意（原文行为，未改）：两个分布图块用 `p` 而非 print(p)，Rscript 顶层
#   自动打印；source() 执行则不出图。
# =============================================================================
setwd('/gpfs1/home/bjx131_pkuhpc/class/TBCR-analysis/results/03_cdr3_length')
cat('=== Part 3: Calculate the distribution of CDR3 length ===\n')
cat('R:', R.version.string, '\n')
cat('time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')

library(ggplot2)
library(dplyr)
library(ggpubr)

## 重建 Part 2 的 air_df / IGH_df（教材为同一会话延续）
air_df <- read.table("/lustre1/share/airr/data/TRUST_output/TCGA-02-0055-01A-01R-1849-01_report.tsv", sep='\t')
colnames(air_df) <- c("count","frequency","CDR3nt","CDR3aa","V","D","J","C")
IGH_df <- air_df %>%
  filter(grepl("^IGHV", V) |
         grepl("^IGHD", D) |
         grepl("^IGHJ", J) |
         grepl("^IGHC", C))

IGH_filter_df <- IGH_df[IGH_df$CDR3aa != "partial",]
print(head(IGH_filter_df))

cat('\n[diag] IGH rows:', nrow(IGH_df),
    ' partial rows:', sum(IGH_df$CDR3aa == "partial"),
    ' kept:', nrow(IGH_filter_df),
    ' kept fraction:', round(nrow(IGH_filter_df)/nrow(IGH_df), 4), '\n')
cat('[diag] kept read count:', sum(IGH_filter_df$count),
    ' (of', sum(IGH_df$count), ')\n')

## count the character lengths
IGH_filter_df <- IGH_filter_df %>%
  mutate(
    nt_length = nchar(CDR3nt),
    aa_length = nchar(CDR3aa)
  )

print(head(IGH_filter_df))

cat('\n[diag] nt length range:', min(IGH_filter_df$nt_length), '-', max(IGH_filter_df$nt_length),
    ' multiples of 3:', all(IGH_filter_df$nt_length %% 3 == 0), '\n')
cat('[diag] aa length range:', min(IGH_filter_df$aa_length), '-', max(IGH_filter_df$aa_length), '\n')

## nt_length
nt_length_df <- IGH_filter_df %>%
  group_by(nt_length) %>%
  summarise(total_count = sum(count))

nt_length_df <- as.data.frame(nt_length_df)

print(head(nt_length_df, 10))
cat('\n[diag] nt length mode:', nt_length_df$nt_length[which.max(nt_length_df$total_count)],
    ' count:', max(nt_length_df$total_count),
    ' total:', sum(nt_length_df$total_count), '\n')
cat('[diag] top5 nt lengths by count:\n')
print(head(nt_length_df[order(-nt_length_df$total_count), ], 5))

p <- ggplot(nt_length_df, aes(x = nt_length, y = total_count)) +
  geom_col(fill = "skyblue", color = "black", alpha = 0.8) +
  labs(title = "Histogram of Nucleotide Length of CDR3 Sequence",
       x = "Nucleotide Length",
       y = "Total Count") +
  theme_minimal()

png("CDR3_nt_length_distribution.png", width = 8, height = 6, units = "in", res = 300)
p
dev.off()

## aa_length
aa_length_df <- IGH_filter_df %>%
  group_by(aa_length) %>%
  summarise(total_count = sum(count))

aa_length_df <- as.data.frame(aa_length_df)

print(head(aa_length_df, 10))
cat('\n[diag] aa length mode:', aa_length_df$aa_length[which.max(aa_length_df$total_count)],
    ' count:', max(aa_length_df$total_count), '\n')
cat('[diag] top5 aa lengths by count:\n')
print(head(aa_length_df[order(-aa_length_df$total_count), ], 5))

p <- ggplot(aa_length_df, aes(x = aa_length, y = total_count)) +
  geom_col(fill = "skyblue", color = "black", alpha = 0.8) +
  labs(title = "Histogram of Amino Acid Length of CDR3 Sequence",
       x = "Amino Acid Length",
       y = "Total Count") +
  theme_minimal()

png("CDR3_aa_length_distribution.png", width = 8, height = 6, units = "in", res = 300)
p
dev.off()

cat('\n=== Part 3 done ===\n')
