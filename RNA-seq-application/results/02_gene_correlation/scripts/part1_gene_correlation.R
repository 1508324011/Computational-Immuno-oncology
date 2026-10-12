## =============================================================================
## Lesson 3 — RNA-seq application: Part 1 Gene correlation analysis
## Source: RNA-seq-application.html (code blocks 6-10, as written)
## Deviations (pre-registered, additive only):
##   - run via Rscript (tutorial: interactive R) under LC_ALL=C
##   - setwd to this step's directory so png() relative paths land here
##   - diagnostic prints: matrix dims, colSums before/after CPM
## =============================================================================
setwd("/gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/02_gene_correlation")

## ---- code block 6: load expression matrix ----
exprsn <- read.table('/lustre1/share/rnaimmune/data/expression/VanAllen2015_CTLA4_Melanoma_RNASeq', sep='\t', header=TRUE, row.names=1)
cat("== exprsn dims ==\n"); print(dim(exprsn))
print(exprsn[1:5,1:5])
cat("== colSums before CPM (head) ==\n"); print(head(colSums(exprsn)))

## ---- code block 7: CPM normalization ----
exprsn <- t(t(exprsn)/colSums(exprsn)*1e6)
cat("== colSums after CPM (head, should be 1e6) ==\n"); print(head(colSums(exprsn)))

## ---- code block 8: GZMA vs PRF1 correlation ----
library(Hmisc)
library(ggplot2)

gene1 <- 'GZMA'
gene2 <- 'PRF1'

pearson_corr <- rcorr(exprsn[gene1,], exprsn[gene2,], type = 'pearson')
spearman_corr <- rcorr(exprsn[gene1,], exprsn[gene2,], type = 'spearman')

cat("== Pearson (GZMA vs PRF1) ==\n"); print(pearson_corr)
cat("== Spearman (GZMA vs PRF1) ==\n"); print(spearman_corr)

## ---- code block 9: GZMA-PRF1 scatter plot ----
df <- data.frame(t(rbind(exprsn[gene1,], exprsn[gene2,])))
colnames(df) <- c('gene1', 'gene2')

pearson_label <- paste0('Pearson r = ', round(pearson_corr$r[1,2], 2), ', p = ', round(pearson_corr$P[1,2], 4))
spearman_label <- paste0('Spearman ρ = ', round(spearman_corr$r[1,2], 2), ', p = ', round(spearman_corr$P[1,2], 4))

# scatter plot
p <- ggplot(df, aes(gene1, gene2))
p <- p + geom_point()
p <- p + annotate("text", x = 50, y = 13, label = pearson_label, size = 5)
p <- p + annotate("text", x = 50, y = 10, label = spearman_label, size = 5)
p <- p + xlab(gene1)
p <- p + ylab(gene2)
p <- p + ggtitle(paste("Pearson & Spearman Correlation of", gene1, 'and', gene2))
png('GZMA_cor_PRF1.png')
print(p)
dev.off()

## ---- code block 10: CTLA4 vs CYT signature correlation ----
geneset <- c('GZMA', 'PRF1')
sig <- colMeans(exprsn[geneset,])
sig_name <- 'CYT'
gene <- 'CTLA4'

pearson_corr <- rcorr(sig, exprsn[gene,], type = 'pearson')
spearman_corr <- rcorr(sig, exprsn[gene,], type = 'spearman')

cat("== Pearson (CTLA4 vs CYT) ==\n"); print(pearson_corr)
cat("== Spearman (CTLA4 vs CYT) ==\n"); print(spearman_corr)

df <- data.frame(t(rbind(sig, exprsn[gene,])))
colnames(df) <- c('sig', 'gene')

pearson_label <- paste0('Pearson r = ', round(pearson_corr$r[1,2], 2), ', p = ', round(pearson_corr$P[1,2], 4))
spearman_label <- paste0('Spearman ρ = ', round(spearman_corr$r[1,2], 2), ', p = ', round(spearman_corr$P[1,2], 4))

# scatter plot
p <- ggplot(df, aes(gene, sig))
p <- p + geom_point()
p <- p + annotate("text", x = 100, y = 20, label = pearson_label, size = 5)
p <- p + annotate("text", x = 100, y = 16, label = spearman_label, size = 5)
p <- p + xlab(gene)
p <- p + ylab(sig_name)
p <- p + ggtitle(paste("Pearson & Spearman Correlation of", gene, 'and', sig_name))
png('CTLA4_cor_CYT.png')
print(p)
dev.off()

cat("== Part 1 done. outputs: GZMA_cor_PRF1.png, CTLA4_cor_CYT.png ==\n")
