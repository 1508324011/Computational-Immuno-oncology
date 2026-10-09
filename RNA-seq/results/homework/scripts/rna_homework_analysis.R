## =============================================================================
## RNA-seq HOMEWORK — R analysis (tutorial step-by-step blocks 21-31, as written)
## Source: RNA_seq_homework.html code cells (quasi-mode quant dirs)
## Deviations (pre-registered, additive only):
##   - run via Rscript (tutorial: interactive R) under LC_ALL=C
##   - diagnostic prints: sample order, kept transcripts, DEG counts
##   - KEGG block kept commented exactly as written;
##     the consolidated "CODES" cell wraps it in ''' (invalid R — doc bug) and
##     still filters "_direct" (mini naming) — step-by-step path used instead
## =============================================================================

## 1.Import the packages
library(readr)     # read data
library(edgeR)     # DEG analysis
library(limma)     # Limma
library(ggplot2)   # ggplot
library(tximport)
library(ggrepel)   # Volcano plot
library(biomaRt)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)  # 用于注释
library(rtracklayer)

## List all quant.sf files — homework filters on "quasi"
salmon_files <- list.files(pattern = "quant\\.sf$", full.names = TRUE, recursive = TRUE)
salmon_files <- salmon_files[grepl("quasi", salmon_files)]
cat("== salmon_files (order feeds the design matrix) ==\n")
print(salmon_files)

## 2.Import data using tximport
txi <- tximport(salmon_files, type = "salmon", txOut = TRUE)
cat("== sample columns (must be HBR,HBR,HBR,UHR,UHR,UHR) ==\n")
print(colnames(txi$counts))

## 3.Create DGEList and filter low express genes
dge <- DGEList(counts = txi$counts)
counts_threshold <- 10
samples_threshold <- 2
keep <- rowSums(dge$counts >= counts_threshold) >= samples_threshold
dge_filtered <- dge[keep, ]
cat("== transcripts kept:", nrow(dge_filtered), "of", nrow(dge), "==\n")

## 4.make sample information
sample_info <- data.frame(
  sample = colnames(dge_filtered),
  condition = factor(c(rep("HBR", 3), rep("UHR", 3)))  # change based on samples
)
print(sample_info)

## 5.design matrix
design <- model.matrix(~condition, data = sample_info)

## 6.DEG analysis edgeR
dge_filtered <- estimateDisp(dge_filtered, design)
fit <- glmFit(dge_filtered, design)
results <- glmLRT(fit)

## 7.extract results
results_table <- topTags(results, n = Inf)$table
results_table$regulate <- ifelse(results_table$FDR < 0.05 & results_table$logFC > 1, "up-regulated",
                                   ifelse(results_table$FDR < 0.05 & results_table$logFC < -1, "down-regulated", "not sig"))
cat("== regulation summary ==\n")
print(table(results_table$regulate, useNA = "ifany"))

## And for gene_mapping ,you can do
gene_anno <- readRDS(
  "/lustre1/share/RNA_seq_class/genome/gene_anno_gencode.v36.rds"
)

rownames_nover <- sub("\\..*$", "", rownames(results_table))
anno_transcript_nover <- sub("\\..*$", "", gene_anno$transcript_id)

match_index <- match(
  rownames_nover,
  anno_transcript_nover
)

results_table$gene_name <- gene_anno$gene_name[match_index]

results_table$label <- ifelse(
  results_table$regulate != "not sig",
  results_table$gene_name,
  NA_character_
)

head(results_table$gene_name)

cat(
  "Matched:",
  sum(!is.na(results_table$gene_name)),
  "/",
  nrow(results_table),
  "\n"
)

## 8.volcano_plot
pdf("volcano_plot.pdf", width = 7, height = 7)
ggplot(results_table, aes(x = logFC, y = -log10(FDR))) +
  geom_point(alpha = 0.5, size = 2, aes(color = regulate)) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_text_repel(data = subset(results_table, regulate != "not sig"),
                  aes(label = label),
                  size = 2,
                  segment.color = "black", show.legend = FALSE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(FDR)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## 9.Function analysis(KEGG and GO)
gene_list <- results_table$gene_name[results_table$regulate == "up-regulated"]
gene_list <- gene_list[!(is.na(gene_list) | gene_list == "")]
cat("== up-regulated genes for enrichment:", length(gene_list), "==\n")

#KEGG analysis
#entrez_ids <- mapIds(org.Hs.eg.db,keys = gene_list,column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
#entrez_ids <- na.omit(entrez_ids)
#kegg_enrich <- enrichKEGG(gene= entrez_ids,organism     = "hsa",keyType      = "kegg",pvalueCutoff = 1,qvalueCutoff = 1)
#kegg_dotplot <- dotplot(kegg_enrich, showCategory = 10, title = "KEGG Pathways")
#ggsave("KEGG_dotplot.pdf", plot = kegg_dotplot, width = 7, height = 5)

#GO analysis
go_enrich <- enrichGO(gene         = gene_list,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "SYMBOL",
                      ont          = "ALL",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

go_dotplot <- dotplot(go_enrich, showCategory = 10, title = "GO Enrichment")
ggsave("GO_dotplot_up_regulate.pdf", plot = go_dotplot, width = 7, height = 5)

## 10.heatmap
logCPM_scaled <- cpm(dge_filtered, log = TRUE, prior.count = 1)
rownames(logCPM_scaled) <- sub("\\..*$", "", rownames(logCPM_scaled))
degs <- results_table[results_table$FDR < 0.05 & abs(results_table$logFC) > 1, ]
top50 <- head(rownames(degs[order(degs$FDR), ]), 50)
top50_clean <- sub("\\..*$", "", top50)
logCPM_scaled_sub <- logCPM_scaled[top50_clean, ]
gene_map <- setNames(gene_anno$gene_name, gene_anno$transcript_id)
rownames(logCPM_scaled_sub) <- gene_map[rownames(logCPM_scaled_sub)]
rownames(sample_info) <- colnames(logCPM_scaled)
pdf("heatmap_DEG.pdf", width = 10, height = 8)
pheatmap(
  logCPM_scaled_sub,
  cluster_rows = TRUE,        # Cluster the rows (genes)
  cluster_cols = TRUE,        # Cluster the columns (samples)
  show_rownames = TRUE,       # Display gene names
  annotation_col = sample_info,  # Add sample condition as annotation
  main = "Heatmap of Differentially Expressed Genes"
)
dev.off()

cat("== R analysis done: volcano_plot.pdf, GO_dotplot_up_regulate.pdf, heatmap_DEG.pdf ==\n")
