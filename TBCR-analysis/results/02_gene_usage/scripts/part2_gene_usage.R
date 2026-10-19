# part2_gene_usage.R — 教材 Part 2 原文（as-written）
# =============================================================================
# 教材：TBCR analysis（Tian Xu / Zexian Zeng, PKU）
# 偏差登记：
#   1. 交互式 R → Rscript
#   2. 教材在 ~/airr_output 下跑，本脚本 setwd 到本步骤目录
#   3. 教材 Part 2 延续 Part 1 的 R 会话；脚本形式在此重新读表（等价）
#   4. LC_ALL=C 固定排序行为
#   5. 加了诊断打印（原文 print 保留）
# 注意（原文行为，未改）：
#   - report.tsv 首行是 "#count\tfrequency\t..."，靠 read.table 默认
#     comment.char="#" 被跳过，故 read.table 不设 header=TRUE 也不会多一行
#   - 热图块写的是 `p` 而非 print(p)，在 Rscript 顶层会自动打印；若用
#     source() 执行则不会出图
# =============================================================================
setwd('/gpfs1/home/bjx131_pkuhpc/class/TBCR-analysis/results/02_gene_usage')
cat('=== Part 2: Obtain gene usage information ===\n')
cat('R:', R.version.string, '\n')
cat('time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')

library(ggplot2)
library(dplyr)
library(ggpubr)

## A. Load CDR3 information
air_df <- read.table("/lustre1/share/airr/data/TRUST_output/TCGA-02-0055-01A-01R-1849-01_report.tsv", sep='\t')
print(head(air_df))

colnames(air_df) <- c("count","frequency","CDR3nt","CDR3aa","V","D","J","C")
print(head(air_df))
cat('\n[diag] report rows:', nrow(air_df), ' sum(count):', sum(air_df$count),
    ' sum(frequency):', sum(air_df$frequency), '\n')

## B. Extract TCR and BCR information individually
## Extract rows of TCR alpha chain
TRA_df <- air_df %>%
  filter(grepl("^TRAV", V) |
         grepl("^TRAJ", J) |
         grepl("^TRAC", C))

print(head(TRA_df))

## Extract rows of BCR heavy chain
IGH_df <- air_df %>%
  filter(grepl("^IGHV", V) |
         grepl("^IGHD", D) |
         grepl("^IGHJ", J) |
         grepl("^IGHC", C))

print(head(IGH_df))

cat('\n[diag] TRA rows:', nrow(TRA_df), ' TRB rows:',
    sum(grepl("^TRBV", air_df$V)), ' IGH rows:', nrow(IGH_df), '\n')
cat('[diag] IGH sum(count):', sum(IGH_df$count),
    ' IGH sum(frequency):', sum(IGH_df$frequency), '\n')

## C. Calculate gene usage frequency
##  Group by V columns, calculate the sum of the counts for each group, and calculate the usage frequency
usage_IGHV_df <- IGH_df %>%
  group_by(V) %>%
  summarise(total_count = sum(count)) %>%
  mutate(usage_frequency = total_count / sum(total_count))

##  For ease of reading, we sort the usage_frequency column in descending order
usage_IGHV_df <- usage_IGHV_df %>%
  arrange(desc(usage_frequency))

usage_IGHV_df <- as.data.frame(usage_IGHV_df)

print(head(usage_IGHV_df))
cat('\n[diag] distinct IGHV genes:', nrow(usage_IGHV_df),
    ' usage sum:', sum(usage_IGHV_df$usage_frequency), '\n')
cat('[diag] top5 IGHV:', paste(usage_IGHV_df$V[1:5], collapse=', '), '\n')

# select the top 10 genes
top_10_IGHV_usage_df <- usage_IGHV_df[1:10,]

p <- ggplot(top_10_IGHV_usage_df, aes(x = reorder(V, -usage_frequency), y = usage_frequency)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Top 10 IGHV Gene Usage Frequency",
       x = "V Gene",
       y = "Usage Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png('top10_IGHV_genes.png')
print(p)
dev.off()

##  Group by J columns, calculate the sum of the counts for each group, and calculate the usage frequency
usage_IGHJ_df <- IGH_df %>%
  group_by(J) %>%
  summarise(total_count = sum(count)) %>%
  mutate(usage_frequency = total_count / sum(total_count))

##  For ease of reading, we sort the usage_frequency column in descending order
usage_IGHJ_df <- usage_IGHJ_df %>%
  arrange(desc(usage_frequency))

usage_IGHJ_df <- as.data.frame(usage_IGHJ_df)

print(head(usage_IGHJ_df))
cat('\n[diag] distinct IGHJ genes:', nrow(usage_IGHJ_df),
    ' usage sum:', sum(usage_IGHJ_df$usage_frequency), '\n')
cat('[diag] top5 IGHJ:', paste(usage_IGHJ_df$J[1:5], collapse=', '), '\n')

# select the top 10 genes
top_10_IGHJ_usage_df <- usage_IGHJ_df[1:10,]

p <- ggplot(top_10_IGHJ_usage_df, aes(x = reorder(J, -usage_frequency), y = usage_frequency)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Top 10 IGHJ Gene Usage Frequency",
       x = "J Gene",
       y = "Usage Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png('top10_IGHJ_genes.png')
print(p)
dev.off()

##  Group by V and J columns, calculate the sum of the counts for each group, and calculate the usage frequency
usage_VJ_df <- IGH_df %>%
  group_by(V, J) %>%
  summarise(total_count = sum(count), .groups = 'drop') %>%   ## remove the combined group column
  mutate(usage_frequency = total_count / sum(total_count)) %>%
  arrange(desc(usage_frequency))

usage_VJ_df <- as.data.frame(usage_VJ_df)

print(head(usage_VJ_df))
cat('\n[diag] distinct V-J pairs:', nrow(usage_VJ_df),
    ' top pair:', usage_VJ_df$V[1], '+', usage_VJ_df$J[1],
    ' count:', usage_VJ_df$total_count[1], '\n')

p <- ggplot(usage_VJ_df, aes(x = V, y = J, fill = total_count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "Heatmap of V and J Gene Counts", x = "V Gene", y = "J Gene") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

png("heatmap_VJ.png", width = 12, height = 6, units = "in", res = 300)
p
dev.off()

cat('\n=== Part 2 done ===\n')
