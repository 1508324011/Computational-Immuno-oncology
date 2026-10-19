# part1_clonotype.R — 教材 Part 1 原文（as-written）
# =============================================================================
# 教材：TBCR analysis（Tian Xu / Zexian Zeng, PKU）
# 偏差登记（仅这些，逻辑逐字照抄）：
#   1. 交互式 R → Rscript（教材为逐块粘贴，状态在同一会话延续）
#   2. 教材在 ~/airr_output 下跑，本脚本 setwd 到本步骤目录
#   3. LC_ALL=C 固定排序/正则行为
#   4. 加了诊断打印（原文 print(head(...)) 保留）
# =============================================================================
setwd('/gpfs1/home/bjx131_pkuhpc/class/TBCR-analysis/results/01_clonotype')
cat('=== Part 1: Obtain clonetype information ===\n')
cat('R:', R.version.string, '\n')
cat('time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')

library(ggplot2)
library(dplyr)
library(ggpubr)

## A. Load consensus assemblies information
consensus_df <- read.table("/lustre1/share/airr/data/TRUST_output/TCGA-02-0055-01A-01R-1849-01_cdr3.out", sep='\t')
print(head(consensus_df))

colnames(consensus_df) <- c("consensus_id","index_within_consensus","V_gene","D_gene","J_gene","C_gene","CDR1","CDR2","CDR3","CDR3_score","read_fragment_count","CDR3_germline_similarity")
print(head(consensus_df))
cat('\n[diag] cdr3.out rows:', nrow(consensus_df),
    ' unique consensus_id:', length(unique(consensus_df$consensus_id)), '\n')

## B. Extract TCR and BCR information individually
## Extract TCR information
TCR_df <- consensus_df %>%
  filter(grepl("^TR", V_gene) |
         grepl("^TR", D_gene) |
         grepl("^TR", J_gene) |
         grepl("^TR", C_gene))

# Extract BCR information
BCR_df <- consensus_df %>%
  filter(grepl("^IG", V_gene) |
         grepl("^IG", D_gene) |
         grepl("^IG", J_gene) |
         grepl("^IG", C_gene))

print(head(TCR_df))
cat('\n[diag] TCR rows:', nrow(TCR_df), ' unique consensus_id:', length(unique(TCR_df$consensus_id)), '\n')
cat('[diag] BCR rows:', nrow(BCR_df), ' unique consensus_id:', length(unique(BCR_df$consensus_id)), '\n')
cat('[diag] rows neither TR nor IG:', nrow(consensus_df) - nrow(TCR_df) - nrow(BCR_df), '\n')

## C. Show the consensus assemblies information of the sample
BCR_consensus_count <- table(factor(BCR_df$consensus_id, levels = unique(BCR_df$consensus_id)))
BCR_consensus_df <- as.data.frame(BCR_consensus_count)
print(head(BCR_consensus_df))

## count the number of different assemblies in the sample
print(nrow(BCR_consensus_df))
cat('\n[diag] BCR consensus assemblies:', nrow(BCR_consensus_df),
    ' total rows:', sum(BCR_consensus_df$Freq), '\n')
cat('[diag] Freq table:\n'); print(table(BCR_consensus_df$Freq))

## count the number of different assemblies in the sample
png('BCR_consensus_frequency_distribution.png')
p <- ggplot(BCR_consensus_df, aes(x = Freq)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Frequency Distribution of Consensus Assemblies",
       x = "Frequency (Freq)",
       y = "Count") +
  theme_minimal()
print(p)
dev.off()

cat('\n=== Part 1 done ===\n')
