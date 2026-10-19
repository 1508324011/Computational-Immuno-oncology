# corrected_analysis.R — TBCR 教材修正版（对拍原样结果）
# =============================================================================
# 本脚本不放宽教材逻辑，只修正四类已实测的口径问题，并把原样口径一并复算用于对拍：
#   C1 diversity 口径：教材 nrow(IGH_df) 是 report 行数（同一 CDR3nt 会因
#      V/C 等位或常量区赋值不同被拆成多行）→ 修正为按 CDR3nt 折叠求和后计数
#   C2 序列有效性：教材 Part 3 只剔 CDR3aa=='partial'，实测还有 'out_of_frame'
#      与含 '?'/'_' 的插补序列 → 修正为两级过滤
#        L1(宽松)：剔 partial + out_of_frame
#        L2(严格)：再剔含 '?'/'_' 的插补序列，只留 20 种标准氨基酸
#   C3 NaN：diversity==1 时教材 clonality = 1 - 0/log2(1) = NaN →
#      修正为显式 NA 并计数报告（不编造数值）
#   C4 对数底数：正文公式写 log（自然对数），代码用 log2 →
#      修正版同时给 nat（log）与 bit（log2），clonality 两者恒等
#   另：C5 基因使用率的 '*' 占位符（未检出 D/J）不作为基因计入分母
# 输出：
#   OV_bcr_diversity_corrected.csv   逐样本 raw + L1 + L2 三套指标
#   corrected_summary.txt            全部关键数字（供文档引用）
#   修正版图：相关性 3 张、CDR3 长度 2 张、gene usage 2 张、raw-vs-corrected 1 张
# =============================================================================
setwd('/gpfs1/home/bjx131_pkuhpc/class/TBCR-analysis/results/05_corrected')
options(stringsAsFactors = FALSE)
suppressMessages({library(ggplot2); library(dplyr); library(ggpubr)})
cat('=== TBCR corrected analysis ===\n')
cat('R:', R.version.string, '\n')
cat('time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')

AA_RE   <- '^[ACDEFGHIKLMNPQRSTVWY]+$'   # 完整、无歧义的氨基酸序列
DATA    <- '/lustre1/share/airr/data'
DEMO    <- file.path(DATA, 'TRUST_output', 'TCGA-02-0055-01A-01R-1849-01_report.tsv')
OVDIR   <- file.path(DATA, 'OV_results')
log_con <- file('corrected_summary.txt', open = 'wt')
say <- function(...) { s <- paste0(...); cat(s, '\n'); writeLines(s, log_con) }

read_report <- function(path) {
  x <- read.table(path, sep = '\t')
  colnames(x) <- c('count','frequency','CDR3nt','CDR3aa','V','D','J','C')
  x
}
is_IGH <- function(x) grepl('^IGHV', x$V) | grepl('^IGHD', x$D) |
                        grepl('^IGHJ', x$J) | grepl('^IGHC', x$C)
seq_class <- function(aa) {
  ifelse(aa == 'partial', 'partial',
  ifelse(aa == 'out_of_frame', 'out_of_frame',
  ifelse(grepl(AA_RE, aa), 'valid',
  ifelse(grepl('[?_]', aa), 'imputed', 'other'))))
}

## ---- 指标函数 ------------------------------------------------------------
# 输入：一个克隆层面的数据框（count 为克隆读数），返回 5 指标
indices_of <- function(counts) {
  counts <- counts[counts > 0]
  N <- length(counts)                      # diversity
  if (N == 0) return(c(diversity=0, CPK=NA, entropy_nat=NA, entropy_bit=NA,
                       clonality=NA, gini=NA))
  tot <- sum(counts)
  p <- counts / tot
  H_nat <- -sum(p * log(p))
  H_bit <- -sum(p * log2(p))
  cpl <- if (N > 1) 1 - H_bit / log2(N) else NA_real_   # N==1 时未定义
  c(diversity = N, CPK = N / tot * 1000, entropy_nat = H_nat, entropy_bit = H_bit,
    clonality = cpl, gini = 1 - sum(p^2))
}
# 原样口径：直接用 report 行（行即克隆，含低质量序列）
metrics_raw <- function(ig) indices_of(as.numeric(ig$count))
# 折叠口径：按 CDR3nt 折叠求和（同一 CDR3 的等位/常量区分行合并）
collapse_cdr3 <- function(ig) {
  if (nrow(ig) == 0) return(numeric(0))
  d <- aggregate(count ~ CDR3nt, data = ig, FUN = sum)
  setNames(d$count, d$CDR3nt)
}

## ---- C1/C2/C3/C4：单样本演示 ---------------------------------------------
ig <- read_report(DEMO) %>% filter(is_IGH(.))
cls <- seq_class(ig$CDR3aa)
say('## demo sample (TCGA-02-0055)')
say('IGH report rows: ', nrow(ig), '  unique CDR3nt: ', length(unique(ig$CDR3nt)),
    '  unique CDR3aa: ', length(unique(ig$CDR3aa)))
for (k in sort(unique(cls))) say('  seq_class ', k, ': ', sum(cls == k))
say('total reads: ', sum(ig$count))

ig_L1 <- ig[cls %in% c('valid','imputed'), ]
ig_L2 <- ig[cls == 'valid', ]
say('L1 kept rows: ', nrow(ig_L1), '  L2 kept rows: ', nrow(ig_L2))

demo <- rbind(
  raw = metrics_raw(ig),
  L1  = indices_of(collapse_cdr3(ig_L1)),
  L2  = indices_of(collapse_cdr3(ig_L2))
)
say('demo indices (rows = raw / L1 / L2):')
say(paste(sprintf('%-6s', rownames(demo)), collapse = ''))
for (j in colnames(demo)) say(sprintf('%-12s ', j), paste(sprintf('%-10.6f', demo[, j]), collapse = ''))

## ---- 375 样本批量 --------------------------------------------------------
files <- list.files(OVDIR)
say('')
say('## cohort (OV_results)')
say('samples: ', length(files))
t0 <- Sys.time()
res <- vector('list', length(files))
for (i in seq_along(files)) {
  x  <- read_report(file.path(OVDIR, files[i]))
  ig_i <- x[is_IGH(x), ]
  cl <- seq_class(ig_i$CDR3aa)
  r0 <- metrics_raw(ig_i)
  r1 <- indices_of(collapse_cdr3(ig_i[cl %in% c('valid','imputed'), ]))
  r2 <- indices_of(collapse_cdr3(ig_i[cl == 'valid', ]))
  res[[i]] <- data.frame(
    sample = tools::file_path_sans_ext(basename(files[i])),
    rows_IGH = nrow(ig_i), reads = sum(ig_i$count),
    n_partial = sum(cl == 'partial'), n_oof = sum(cl == 'out_of_frame'),
    n_imputed = sum(cl == 'imputed'), n_valid = sum(cl == 'valid'),
    raw_diversity = r0['diversity'], raw_CPK = r0['CPK'], raw_entropy_bit = r0['entropy_bit'],
    raw_clonality = r0['clonality'], raw_gini = r0['gini'],
    L1_diversity = r1['diversity'], L1_CPK = r1['CPK'], L1_entropy_nat = r1['entropy_nat'],
    L1_entropy_bit = r1['entropy_bit'], L1_clonality = r1['clonality'], L1_gini = r1['gini'],
    L2_diversity = r2['diversity'], L2_CPK = r2['CPK'], L2_entropy_nat = r2['entropy_nat'],
    L2_entropy_bit = r2['entropy_bit'], L2_clonality = r2['clonality'], L2_gini = r2['gini'],
    check.names = FALSE)
}
m <- do.call(rbind, res)
say('loop elapsed: ', round(as.numeric(difftime(Sys.time(), t0, units = 'secs')), 1), ' s')

## ---- 对拍统计 ------------------------------------------------------------
cmp <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  c(n_finite = sum(ok),
    rho = suppressWarnings(cor(a[ok], b[ok], method = 'spearman')),
    med_ratio = median(b[ok] / a[ok]),
    n_changed = sum(a[ok] != b[ok]))
}
say('')
say('## as-written vs corrected')
for (nm in c('diversity','CPK','clonality','gini')) {
  a <- m[[paste0('raw_', nm)]]; b <- m[[paste0('L1_', nm)]]
  cc <- cmp(a, b)
  say(sprintf('%-11s raw=median %.4f  L1=median %.4f  finite=%d  spearman=%.4f  median(L1/raw)=%.4f  changed=%d',
              nm, median(a, na.rm = TRUE), median(b, na.rm = TRUE),
              cc['n_finite'], cc['rho'], cc['med_ratio'], cc['n_changed']))
}
say(sprintf('entropy bit->nat factor (log 2 = %.6f), raw bit median %.4f -> L1 nat median %.4f',
            log(2), median(m$raw_entropy_bit, na.rm = TRUE), median(m$L1_entropy_nat, na.rm = TRUE)))
say('')
say('NA counts:')
for (nm in c('raw_clonality','L1_clonality','L2_clonality','L1_diversity','L2_diversity',
             'L1_entropy_nat','L1_CPK','L1_gini')) {
  say(sprintf('  %-16s NA = %d', nm, sum(is.na(m[[nm]]))))
}
say('samples that are degenerate (L1 diversity < 2): ',
    sum(m$L1_diversity < 2, na.rm = TRUE),
    ' -> ', paste(m$sample[m$L1_diversity < 2], collapse = ', '))
say('samples with L1 diversity == 0 (no valid IGH): ', sum(m$L1_diversity == 0, na.rm = TRUE),
    ' -> ', paste(m$sample[m$L1_diversity == 0], collapse = ', '))
say('')
say('ranges (raw -> L1):')
for (nm in c('diversity','CPK','clonality','gini')) {
  a <- m[[paste0('raw_', nm)]]; b <- m[[paste0('L1_', nm)]]
  say(sprintf('  %-11s %10.4f..%10.4f  ->  %10.4f..%10.4f',
              nm, min(a, na.rm = TRUE), max(a, na.rm = TRUE),
              min(b, na.rm = TRUE), max(b, na.rm = TRUE)))
}
say(sprintf('  %-11s %10.4f..%10.4f  (L2)', 'diversity',
            min(m$L2_diversity, na.rm = TRUE), max(m$L2_diversity, na.rm = TRUE)))

## Spearman 相关性：原样 vs 修正
say('')
say('## spearman (as-written)')
say(sprintf('  diversity~entropy  %.4f', suppressWarnings(cor(m$raw_diversity, m$raw_entropy_bit, method='spearman', use='complete.obs'))))
say(sprintf('  CPK~entropy        %.4f', suppressWarnings(cor(m$raw_CPK, m$raw_entropy_bit, method='spearman', use='complete.obs'))))
say(sprintf('  CPK~clonality      %.4f', suppressWarnings(cor(m$raw_CPK, m$raw_clonality, method='spearman', use='complete.obs'))))
say('## spearman (corrected L1, entropy in nat)')
say(sprintf('  diversity~entropy  %.4f', suppressWarnings(cor(m$L1_diversity, m$L1_entropy_nat, method='spearman', use='complete.obs'))))
say(sprintf('  CPK~entropy        %.4f', suppressWarnings(cor(m$L1_CPK, m$L1_entropy_nat, method='spearman', use='complete.obs'))))
say(sprintf('  CPK~clonality      %.4f', suppressWarnings(cor(m$L1_CPK, m$L1_clonality, method='spearman', use='complete.obs'))))
say('## spearman (corrected L2, entropy in nat)')
say(sprintf('  diversity~entropy  %.4f', suppressWarnings(cor(m$L2_diversity, m$L2_entropy_nat, method='spearman', use='complete.obs'))))
say(sprintf('  CPK~entropy        %.4f', suppressWarnings(cor(m$L2_CPK, m$L2_entropy_nat, method='spearman', use='complete.obs'))))
say(sprintf('  CPK~clonality      %.4f', suppressWarnings(cor(m$L2_CPK, m$L2_clonality, method='spearman', use='complete.obs'))))

## ---- 落盘 ----------------------------------------------------------------
write.csv(m, 'OV_bcr_diversity_corrected.csv', row.names = FALSE)

## ---- 修正版图 ------------------------------------------------------------
# 1. raw vs corrected diversity（log-log）
p <- ggplot(m, aes(x = raw_diversity, y = L1_diversity)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed') +
  scale_x_log10() + scale_y_log10() +
  labs(title = 'Raw (report rows) vs Corrected (unique CDR3nt, L1) diversity',
       x = 'Raw diversity = nrow(IGH_df)', y = 'Corrected diversity') +
  theme_classic()
png('raw_vs_corrected_diversity.png', width = 7, height = 6, units = 'in', res = 300)
print(p); dev.off()

# 2. 三张相关性图（修正 L1，entropy = nat，避免 aes 里引用 $ 列）
mk_cor <- function(xv, yv, xl, yl, fname) {
  p <- ggplot(m, aes(x = .data[[xv]], y = .data[[yv]])) +
    geom_point(size = 1.2, alpha = 0.6) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') +
    theme_classic() +
    labs(x = xl, y = yl) +
    stat_cor(method = 'spearman')
  png(fname, width = 7, height = 6, units = 'in', res = 300)
  print(p); dev.off()
}
mk_cor('L1_diversity','L1_entropy_nat','Diversity (unique CDR3nt, L1)','Entropy (nat)',
       'corrected_diversity_cor_entropy.png')
mk_cor('L1_CPK','L1_entropy_nat','CPK (L1)','Entropy (nat)','corrected_CPK_cor_entropy.png')
mk_cor('L1_CPK','L1_clonality','CPK (L1)','Clonality (L1)','corrected_CPK_cor_clonality.png')

## ---- C2 修正后 CDR3 长度分布（demo 样本） --------------------------------
say('')
say('## CDR3 length, demo sample')
len_df <- ig_L2 %>%
  mutate(nt_length = nchar(CDR3nt), aa_length = nchar(CDR3aa))
say('L2 rows: ', nrow(len_df), '  nt multiple-of-3: ', all(len_df$nt_length %% 3 == 0))
say('nt length range: ', min(len_df$nt_length), '-', max(len_df$nt_length))
say('aa length range: ', min(len_df$aa_length), '-', max(len_df$aa_length))
ntd <- len_df %>% group_by(nt_length) %>% summarise(total_count = sum(count), .groups = 'drop')
aad <- len_df %>% group_by(aa_length) %>% summarise(total_count = sum(count), .groups = 'drop')
say('nt mode: ', ntd$nt_length[which.max(ntd$total_count)], ' (', max(ntd$total_count), ' reads)')
say('aa mode: ', aad$aa_length[which.max(aad$total_count)], ' (', max(aad$total_count), ' reads)')
p <- ggplot(ntd, aes(x = nt_length, y = total_count)) +
  geom_col(fill = 'skyblue', colour = 'black', alpha = 0.8) +
  labs(title = 'CDR3 Nucleotide Length (corrected: valid sequences only)',
       x = 'Nucleotide Length', y = 'Total Count') + theme_minimal()
png('corrected_CDR3_nt_length_distribution.png', width = 8, height = 6, units = 'in', res = 300)
print(p); dev.off()
p <- ggplot(aad, aes(x = aa_length, y = total_count)) +
  geom_col(fill = 'skyblue', colour = 'black', alpha = 0.8) +
  labs(title = 'CDR3 Amino Acid Length (corrected: valid sequences only)',
       x = 'Amino Acid Length', y = 'Total Count') + theme_minimal()
png('corrected_CDR3_aa_length_distribution.png', width = 8, height = 6, units = 'in', res = 300)
print(p); dev.off()

## ---- C5 基因使用率排除 '*' -----------------------------------------------
say('')
say('## gene usage, demo sample (C5: exclude * placeholder)')
say('IGH rows with J == "*": ', sum(ig$J == '*'),
    '  V == "*": ', sum(ig$V == '*'),
    '  D == "*": ', sum(ig$D == '*'),
    '  C == "*": ', sum(ig$C == '*'))
usage <- function(v) {
  u <- data.frame(gene = v, count = ig$count) %>% filter(gene != '*') %>%
    group_by(gene) %>% summarise(total_count = sum(count), .groups = 'drop') %>%
    mutate(usage_frequency = total_count / sum(total_count)) %>% arrange(desc(usage_frequency))
  as.data.frame(u)
}
v_all <- data.frame(gene = ig$V, count = ig$count) %>% group_by(gene) %>%
  summarise(total_count = sum(count), .groups = 'drop') %>% mutate(f = total_count / sum(total_count)) %>% arrange(desc(f))
say('IGHV distinct incl "*": ', nrow(v_all), '  excluding "*": ', nrow(usage(ig$V)))
say('top5 IGHV (excl *): ', paste(head(usage(ig$V)$gene, 5), collapse = ', '))
say('top5 IGHJ (excl *): ', paste(head(usage(ig$J)$gene, 5), collapse = ', '))
uV <- usage(ig$V); uJ <- usage(ig$J)
p <- ggplot(head(uV, 10), aes(x = reorder(gene, -usage_frequency), y = usage_frequency)) +
  geom_bar(stat = 'identity', fill = 'skyblue', colour = 'black', alpha = 0.7) +
  labs(title = 'Top 10 IGHV Usage (corrected: "*" excluded)', x = 'V Gene', y = 'Usage Frequency') +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
png('corrected_top10_IGHV_genes.png', width = 8, height = 6, units = 'in', res = 300)
print(p); dev.off()
p <- ggplot(head(uJ, 10), aes(x = reorder(gene, -usage_frequency), y = usage_frequency)) +
  geom_bar(stat = 'identity', fill = 'skyblue', colour = 'black', alpha = 0.7) +
  labs(title = 'Top 10 IGHJ Usage (corrected: "*" excluded)', x = 'J Gene', y = 'Usage Frequency') +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
png('corrected_top10_IGHJ_genes.png', width = 8, height = 6, units = 'in', res = 300)
print(p); dev.off()

## ---- 扩展：TRB（仅统计，不画图）------------------------------------------
say('')
say('## extension: TRB (beta chain) diversity, same C1/C2/C3 rules')
trb_rows <- 0; trb_ok <- 0
trb <- lapply(files, function(f) {
  x <- read_report(file.path(OVDIR, f))
  tb <- x[grepl('^TRBV', x$V) | grepl('^TRBD', x$D) | grepl('^TRBJ', x$J) | grepl('^TRBC', x$C), ]
  cl <- seq_class(tb$CDR3aa)
  r <- indices_of(collapse_cdr3(tb[cl == 'valid', ]))
  c(diversity = unname(r['diversity']), CPK = unname(r['CPK']),
    entropy_nat = unname(r['entropy_nat']), clonality = unname(r['clonality']),
    gini = unname(r['gini']))
})
trb <- do.call(rbind, trb)
say('samples with >=1 valid TRB clone: ', sum(trb[, 'diversity'] > 0, na.rm = TRUE), ' / ', nrow(trb))
say('TRB diversity median (nonzero): ', median(trb[trb[, 'diversity'] > 0, 'diversity'], na.rm = TRUE))
say('TRB diversity range: ', min(trb[, 'diversity'], na.rm = TRUE), '-', max(trb[, 'diversity'], na.rm = TRUE))
write.csv(data.frame(sample = tools::file_path_sans_ext(basename(files)), trb, check.names = FALSE),
          'OV_trb_diversity_corrected.csv', row.names = FALSE)

say('')
say('=== corrected analysis done ===')
close(log_con)
