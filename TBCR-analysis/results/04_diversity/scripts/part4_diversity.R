# part4_diversity.R — 教材 Part 4 原文（as-written）
# =============================================================================
# 教材：TBCR analysis（Tian Xu / Zexian Zeng, PKU）
# 偏差登记：
#   1. 交互式 R → Rscript
#   2. 教材在 ~/airr_output 下跑，本脚本 setwd 到本步骤目录
#   3. LC_ALL=C；加诊断打印（原文 print 保留）
# 注意（原文行为，未改）：
#   - 公式正文写 log（自然对数），代码用 log2 → entropy 单位为 bit
#   - 循环里不剔除 partial / out_of_frame 序列（与 Part 3 口径不同）
#   - diversity = nrow(IGH_df) 即 report 行数，不是唯一 CDR3 数
#   - diversity == 1 的样本 clonality = 1 - 0/0 = NaN
# =============================================================================
setwd('/gpfs1/home/bjx131_pkuhpc/class/TBCR-analysis/results/04_diversity')
cat('=== Part 4: Calculate TCR/BCR diversity index of each sample ===\n')
cat('R:', R.version.string, '\n')
cat('time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')

library(ggplot2)
library(dplyr)
library(ggpubr)

## A. Test on previous sample
air_df <- read.table("/lustre1/share/airr/data/TRUST_output/TCGA-02-0055-01A-01R-1849-01_report.tsv", sep='\t')
colnames(air_df) <- c("count","frequency","CDR3nt","CDR3aa","V","D","J","C")
print(head(air_df))

IGH_df <- air_df %>%
  filter(grepl("^IGHV", V) |
         grepl("^IGHD", D) |
         grepl("^IGHJ", J) |
         grepl("^IGHC", C))

print(head(IGH_df))

diversity <- nrow(IGH_df)
print(diversity)

CPK <- diversity/sum(IGH_df$count) *1000
print(CPK)

freqs <- IGH_df$count/sum(IGH_df$count)
entropy <- -sum(freqs * log2(freqs))
print(entropy)

clonality <- 1 - (entropy/log2(diversity))
print(clonality)

gini <- 1 - sum(freqs * freqs)
print(gini)

cat('\n[diag] single-sample indices: diversity=', diversity,
    ' unique CDR3nt=', length(unique(IGH_df$CDR3nt)),
    ' unique CDR3aa=', length(unique(IGH_df$CDR3aa)), '\n')
cat('[diag] partial rows in IGH:', sum(IGH_df$CDR3aa == "partial"),
    ' out_of_frame:', sum(IGH_df$CDR3aa == "out_of_frame"), '\n')
cat('[diag] total reads:', sum(IGH_df$count), '\n')

## B. Test on cohort samples
## create a data.frame to storage the indices value of each sample
bcr_metrics_df <- data.frame()

## list the files under OV_results folder
files <- list.files("/lustre1/share/airr/data/OV_results")

cat('\n[diag] OV files:', length(files), '\n')
t0 <- Sys.time()

## Following the code above, the indices for each sample are
## calculated in a loop and stored in bcr_metrics_df
for (file in files) {
  air_df <- read.table(paste0("/lustre1/share/airr/data/OV_results/",file), sep="\t")
  colnames(air_df) <- c("count","frequency","CDR3nt","CDR3aa","V","D","J","C")
  IGH_df <- air_df %>%
  filter(grepl("^IGHV", V) |
         grepl("^IGHD", D) |
         grepl("^IGHJ", J) |
         grepl("^IGHC", C))

  diversity <- nrow(IGH_df)
  CPK <- diversity/sum(IGH_df$count) *1000

  freqs <- IGH_df$count/sum(IGH_df$count)
  entropy <- -sum(freqs * log2(freqs))

  clonality <- 1 - (entropy/log2(diversity))
  gini <- 1 - sum(freqs * freqs)

  bcr_metrics_row <- data.frame(diversity, CPK, entropy, clonality, gini)
  rownames(bcr_metrics_row) <- tools::file_path_sans_ext(basename(file))  # remove the file suffix

  bcr_metrics_df <- rbind(bcr_metrics_df, bcr_metrics_row)
}

print(head(bcr_metrics_df))
write.csv(bcr_metrics_df, "OV_bcr_diversity.csv")

cat('\n[diag] loop elapsed:', round(as.numeric(difftime(Sys.time(), t0, units='secs')), 1), 's\n')
cat('[diag] samples:', nrow(bcr_metrics_df), '\n')
cat('[diag] diversity range:', min(bcr_metrics_df$diversity), '-', max(bcr_metrics_df$diversity),
    ' median:', median(bcr_metrics_df$diversity), '\n')
cat('[diag] CPK range:', round(min(bcr_metrics_df$CPK),4), '-', round(max(bcr_metrics_df$CPK),4),
    ' median:', round(median(bcr_metrics_df$CPK),4), '\n')
cat('[diag] entropy range:', round(min(bcr_metrics_df$entropy),4), '-', round(max(bcr_metrics_df$entropy),4), '\n')
cat('[diag] clonality range:', round(min(bcr_metrics_df$clonality, na.rm=TRUE),4), '-',
    round(max(bcr_metrics_df$clonality, na.rm=TRUE),4), '\n')
cat('[diag] gini range:', round(min(bcr_metrics_df$gini),4), '-', round(max(bcr_metrics_df$gini),4), '\n')
cat('[diag] NaN per column:\n'); print(colSums(is.na(bcr_metrics_df)))
cat('[diag] rows with NaN:\n'); print(rownames(bcr_metrics_df)[!complete.cases(bcr_metrics_df)])

## C. Indices correlation
p <- ggplot(bcr_metrics_df,aes(x = bcr_metrics_df$diversity,
                y = bcr_metrics_df$entropy))+
  geom_point()+
  geom_smooth(method = "lm",
              se = F,
              color = "red")+
  theme_classic()+
  stat_cor(data = bcr_metrics_df,
           method = "spearman")

png('diveristy_cor_entropy.png')
print(p)
dev.off()

p <- ggplot(bcr_metrics_df,aes(x = bcr_metrics_df$CPK,
                y = bcr_metrics_df$entropy))+
  geom_point()+
  geom_smooth(method = "lm",
              se = F,
              color = "red")+
  theme_classic()+
  stat_cor(data = bcr_metrics_df,
           method = "spearman")

png('CPK_cor_entropy.png')
print(p)
dev.off()

p <- ggplot(bcr_metrics_df,aes(x = bcr_metrics_df$CPK,
                y = bcr_metrics_df$clonality))+
  geom_point()+
  geom_smooth(method = "lm",
              se = F,
              color = "red")+
  theme_classic()+
  stat_cor(data = bcr_metrics_df,
           method = "spearman")

png('CPK_cor_clonality.png')
print(p)
dev.off()

# 诊断：独立复算 Spearman（与图内 stat_cor 对拍）
cat('\n[diag] spearman rho (diversity~entropy):',
    round(cor(bcr_metrics_df$diversity, bcr_metrics_df$entropy, method='spearman', use='complete.obs'), 4), '\n')
cat('[diag] spearman rho (CPK~entropy):',
    round(cor(bcr_metrics_df$CPK, bcr_metrics_df$entropy, method='spearman', use='complete.obs'), 4), '\n')
cat('[diag] spearman rho (CPK~clonality):',
    round(cor(bcr_metrics_df$CPK, bcr_metrics_df$clonality, method='spearman', use='complete.obs'), 4), '\n')

cat('\n=== Part 4 done ===\n')
