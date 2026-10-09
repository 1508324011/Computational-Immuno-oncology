## Correction pass — heatmap only. Same as scripts/rna_homework_analysis.R step 10
## EXCEPT one line (marked FIX): gene_map keys must be version-stripped to match
## the version-stripped rownames (homework block uses versioned keys -> all-NA labels,
## recorded in log/deviations.log; output: heatmap_DEG_corrected.pdf)
library(edgeR)
library(tximport)
library(pheatmap)

salmon_files <- list.files(pattern = "quant\\.sf$", full.names = TRUE, recursive = TRUE)
salmon_files <- salmon_files[grepl("quasi", salmon_files)]
txi <- tximport(salmon_files, type = "salmon", txOut = TRUE)
dge <- DGEList(counts = txi$counts)
keep <- rowSums(dge$counts >= 10) >= 2
dge_filtered <- dge[keep, ]
sample_info <- data.frame(
  sample = colnames(dge_filtered),
  condition = factor(c(rep("HBR", 3), rep("UHR", 3)))
)
rownames(sample_info) <- colnames(dge_filtered)

logCPM_scaled <- cpm(dge_filtered, log = TRUE, prior.count = 1)
rownames(logCPM_scaled) <- sub("\\..*$", "", rownames(logCPM_scaled))

## recompute DEG order (same definition as homework)
dge_filtered <- estimateDisp(dge_filtered, model.matrix(~condition, data = sample_info))
results <- glmLRT(glmFit(dge_filtered, model.matrix(~condition, data = sample_info)))
results_table <- topTags(results, n = Inf)$table
degs <- results_table[results_table$FDR < 0.05 & abs(results_table$logFC) > 1, ]
top50 <- head(rownames(degs[order(degs$FDR), ]), 50)
top50_clean <- sub("\\..*$", "", top50)
logCPM_scaled_sub <- logCPM_scaled[top50_clean, ]

## FIX: version-stripped keys (homework: setNames(gene_anno$gene_name, gene_anno$transcript_id) -> all NA)
gene_anno <- readRDS("/lustre1/share/RNA_seq_class/genome/gene_anno_gencode.v36.rds")
gene_map <- setNames(gene_anno$gene_name, sub("\\..*$", "", gene_anno$transcript_id))
rownames(logCPM_scaled_sub) <- gene_map[rownames(logCPM_scaled_sub)]
cat("== corrected heatmap row labels (should be gene names, no NA) ==\n")
print(head(rownames(logCPM_scaled_sub)))
cat("NA labels:", sum(is.na(rownames(logCPM_scaled_sub))), "of", nrow(logCPM_scaled_sub), "\n")

pdf("heatmap_DEG_corrected.pdf", width = 10, height = 8)
pheatmap(
  logCPM_scaled_sub,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  annotation_col = sample_info,
  main = "Heatmap of Differentially Expressed Genes (corrected labels)"
)
dev.off()
cat("== done: heatmap_DEG_corrected.pdf ==\n")
