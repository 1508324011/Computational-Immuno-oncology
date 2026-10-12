## =============================================================================
## Lesson 3 — RNA-seq application: Part 3 Immune infiltration analysis
## Source: RNA-seq-application.html (code blocks 17-26, as written)
## Deviations (pre-registered, additive only):
##   - run via Rscript (tutorial: interactive R) under LC_ALL=C
##   - setwd to this step's directory so png()/CIBERSORT outputs land here
##   - re-derive exprsn & metadata$Response here (separate R process)
##   - write each deconvolution method's result to deconv_methods/*.tsv
##     (tutorial only printed res_qts; saving all is additive, non-destructive)
##   - MCPcounter (1/7) wrapped in tryCatch: conda 1.1.0 hardcodes markers as
##     http://raw.githubusercontent.com fetch (login node resets http);
##     per as-written instruction it is recorded as blocked, not bypassed.
##     Six remaining methods run exactly as written.
##   - diagnostic prints: deconvolution_methods, each method's dim
## =============================================================================
setwd("/gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/04_immune_infiltration")
options(stringsAsFactors = FALSE)

## ---- code block 17: load libs, set CIBERSORT ----
library(tidyverse)
library(pheatmap)
library(immunedeconv)
set_cibersort_binary("/lustre1/share/rnaimmune/CIBERSORT.R")
set_cibersort_mat("/lustre1/share/rnaimmune/LM22.txt")

## ---- re-derive exprsn (Part 1 blocks 6,7) ----
exprsn <- read.table('/lustre1/share/rnaimmune/data/expression/VanAllen2015_CTLA4_Melanoma_RNASeq', sep='\t', header=TRUE, row.names=1)
exprsn <- t(t(exprsn)/colSums(exprsn)*1e6)

## ---- code block 18: list methods ----
cat("== deconvolution_methods ==\n"); print(deconvolution_methods)

## ---- code block 19: run 7 methods (MCPcounter skipped — see note) ----
## NOTE: MCPcounter 1.1.0 (conda) hardcodes markers as
##   read.table(curl:::curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"))
## i.e. an *http* (not https) runtime fetch of signature genes from GitHub.
## Login node resets http to raw.githubusercontent.com (Connection reset by peer);
## https reaches HTTP 200. Per instruction to run the tutorial code as written
## (no local-copy override of the genes=/probesets= defaults), MCPcounter is a
## hard blocker and is recorded as such; the remaining six methods run as written.
res_mcp <- tryCatch({
  cat("\n== [1/7] mcp_counter (as written, expected to fail on http fetch) ==\n")
  deconvolute(exprsn, "mcp_counter")
}, error=function(e){
  cat(sprintf("   -> BLOCKED (recorded): %s\n", conditionMessage(e)))
  writeLines(sprintf("[%s] MCPcounter blocked: %s\n  conda MCPcounter 1.1.0 hardcodes http://raw.githubusercontent.com markers; login node resets http. https reachable but tutorial code as-written uses the http default.", Sys.time(), conditionMessage(e)), "mcp_counter_blocked.txt")
  NULL
})
if(!is.null(res_mcp)){ cat("dim:"); print(dim(res_mcp)); write_tsv(as_tibble(res_mcp), "deconv_methods/res_mcp_counter.tsv") } else {
  writeLines("NULL — see mcp_counter_blocked.txt", "deconv_methods/res_mcp_counter.tsv")
}

cat("\n== [2/7] epic ==\n")
res_epc <- deconvolute(exprsn, "epic")
cat("dim:"); print(dim(res_epc))
write_tsv(as_tibble(res_epc), "deconv_methods/res_epic.tsv")

cat("\n== [3/7] quantiseq ==\n")
res_qts <- deconvolute(exprsn, "quantiseq")
cat("dim:"); print(dim(res_qts))
write_tsv(as_tibble(res_qts), "deconv_methods/res_quantiseq.tsv")

cat("\n== [4/7] xcell ==\n")
res_xcl <- deconvolute(exprsn, "xcell")
cat("dim:"); print(dim(res_xcl))
write_tsv(as_tibble(res_xcl), "deconv_methods/res_xcell.tsv")

cat("\n== [5/7] cibersort ==\n")
res_cbs <- deconvolute(exprsn, "cibersort") # will save results in the work dir
cat("dim:"); print(dim(res_cbs))
write_tsv(as_tibble(res_cbs), "deconv_methods/res_cibersort.tsv")

cat("\n== [6/7] cibersort_abs ==\n")
res_cba <- deconvolute(exprsn, "cibersort_abs")
cat("dim:"); print(dim(res_cba))
write_tsv(as_tibble(res_cba), "deconv_methods/res_cibersort_abs.tsv")

cat("\n== [7/7] timer ==\n")
res_tmr <- deconvolute(exprsn, "timer", indications = rep('skcm', ncol(exprsn)))
cat("dim:"); print(dim(res_tmr))
write_tsv(as_tibble(res_tmr), "deconv_methods/res_timer.tsv")

## ---- code block 20: quanTIseq with rmgenes ----
cat("\n== quanTIseq rmgenes (overrides res_qts) ==\n")
res_qts <- deconvolute(
    exprsn, method='quantiseq',
    rmgenes = c('TCL1A', 'TCF4', 'CD37', 'SPIB', 'IRF8', 'BCL11A'))
cat("dim:"); print(dim(res_qts))
write_tsv(as_tibble(res_qts), "deconv_methods/res_quantiseq_rmgenes.tsv")

## ---- code block 21: timer_available_cancers ----
cat("\n== timer_available_cancers ==\n"); print(timer_available_cancers)

## ---- code block 22: print res_qts ----
cat("\n== res_qts (rmgenes version) ==\n"); print(res_qts)

## ---- code block 23: quanTIseq bar plot (stacked, fractions) ----
df_plot <- res_qts[,1:10] %>% gather(sample, fraction, -cell_type) # only show 10 samples
print(df_plot)

# bar plot
p <- ggplot(aes(x=sample, y=fraction, fill=cell_type), data=df_plot)
p <- p + geom_bar(stat='identity')
p <- p + coord_flip()
p <- p + scale_fill_brewer(palette="Paired")
png('qts_barplot.png')
print(p)
dev.off()

## ---- code block 24: quanTIseq scatter/dot plot (relative scores) ----
df_plot <- res_qts[,1:10] %>% gather(sample, score, -cell_type) # only show 9 samples

# scatter plot
p <- ggplot(aes(x=sample, y=score, color=cell_type), data=df_plot)
p <- p + geom_point()
p <- p + facet_wrap(~cell_type, scales='free_y', nrow=2)
p <- p + scale_color_brewer(palette="Paired", guide='none')
p <- p + theme(axis.text = element_text(angle=90, vjust=0.5, hjust=1))
png('qts_scatterplot.png', width = 960)
print(p)
dev.off()

## ---- re-derive metadata + Response factor (Part 2 blocks 11,12,14) ----
metadata <- read.table('/lustre1/share/rnaimmune/data/follow_up/VanAllen2015_CTLA4_Melanoma_RNASeq', sep='\t', header=TRUE, row.names=1)
metadata <- metadata[colnames(exprsn),]
metadata$CYT <- colMeans(exprsn[c('GZMA','PRF1'),])
metadata <- metadata %>% mutate(CYT_level = ifelse(CYT > median(CYT), "high", "low"))
metadata$CTLA4 <- exprsn['CTLA4',]
metadata$Gender <- factor(metadata$Gender, levels=c(0,1)); levels(metadata$Gender) <- c('Female','Male')

## ---- code block 25: Response factor + heatmap (unscaled) ----
print(metadata$Response)
metadata$Response <- factor(metadata$Response, levels=c(0,1))
levels(metadata$Response) <- c('Non-responder', 'Responder')
print(metadata$Response)

head(res_qts)
res_qts_heatmap <- data.frame(res_qts)
rownames(res_qts_heatmap) <- res_qts$cell_type
res_qts_heatmap <- res_qts_heatmap[,-1]
head(res_qts_heatmap)

# heatmap plot
png('qts_heatmap.png', width = 800, height=1000)
pheatmap(res_qts_heatmap,
         scale = 'none',
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         annotation_col = metadata
)
dev.off()

## ---- code block 26: heatmap row-scaled ----
# heatmap plot
# scale each row for easy comparison between samples
png('qts_heatmap_row_scaled.png', width = 800, height=1000)
pheatmap(res_qts_heatmap,
         scale = 'row',
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         annotation_col = metadata
)
dev.off()

cat("\n== Part 3 done. outputs: qts_barplot/scatterplot/heatmap/heatmap_row_scaled.png + deconv_methods/*.tsv ==\n")
