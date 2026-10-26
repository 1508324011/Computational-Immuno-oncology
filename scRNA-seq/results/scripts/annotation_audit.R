#!/usr/bin/env Rscript
# annotation_audit.R -- post-hoc audit of the tutorial's two hard-coded annotation lists.
#
# The tutorial gives the cluster -> cell-type mapping twice, in two places that do not
# agree with each other:
#   * cell 39 (markdown-ish table): 0 Naive CD4 T, 1 CD14+ Mono, 2 Memory CD4 T, ...
#   * cell 43 (new.cluster.ids):    0 Naive CD4 T, 1 Memory CD4 T, 2 CD14+ Mono, ...
# Which one is right depends on how Louvain happened to number the communities in that
# particular run, so we list both next to the observed top markers and let the data decide.
#
# usage: Rscript annotation_audit.R <pbmc3k_final.rds> <markers_all_clusters.csv> <out_csv>
args <- commandArgs(trailingOnly = TRUE)
suppressMessages({library(Seurat); library(dplyr)})

seu <- readRDS(args[[1]])
mk  <- read.csv(args[[2]], stringsAsFactors = FALSE)
mk$cluster <- as.character(mk$cluster)

table_label <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
                 "FCGR3A+ Mono", "NK", "DC", "Platelet")          # tutorial cell 39
rename_label <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T",
                  "FCGR3A+ Mono", "NK", "DC", "Platelet")         # tutorial cell 43

clusters <- sort(unique(as.character(seu$seurat_clusters)), method = "radix")
sig <- mk %>% dplyr::filter(p_val_adj < 0.05, avg_log2FC > 1)

out <- data.frame(
  cluster        = clusters,
  n_cells        = as.integer(table(factor(seu$seurat_clusters, levels = clusters))),
  table_label    = table_label[match(clusters, as.character(seq_along(table_label) - 1))],
  rename_label   = rename_label[match(clusters, as.character(seq_along(rename_label) - 1))],
  top5_markers   = sapply(clusters, function(cl)
                     paste(head(sig$gene[sig$cluster == cl], 5), collapse = ", ")),
  stringsAsFactors = FALSE)

# A one-word verdict per cluster from the panel we can check unambiguously: which of the
# unambiguous programmes (B, monocyte, DC, platelet, NK, CD8, CD4) dominates.
panel <- list(
  CD4    = c("IL7R", "CCR7", "LEF1", "LTB", "IL32"),
  CD8    = c("CD8A", "CD8B", "GZMK", "CCL5", "CTSW"),
  NK     = c("GNLY", "GZMB", "KLRF1", "FGFBP2", "SPON2"),
  B      = c("MS4A1", "CD79A", "CD79B", "TCL1A"),
  Mono   = c("CD14", "LYZ", "S100A8", "S100A9", "FCN1"),
  FcMono = c("FCGR3A", "MS4A7", "LST1"),
  DC     = c("FCER1A", "CLEC10A", "SERPINF1", "CST3"),
  Plat   = c("PPBP", "PF4", "GP9"))
data <- LayerData(seu, assay = "RNA", layer = "data")
programme <- sapply(clusters, function(cl) {
  cells <- which(as.character(seu$seurat_clusters) == cl)
  score <- sapply(panel, function(g) {
    g <- intersect(g, rownames(data))
    mean(as.matrix(data[g, cells, drop = FALSE]))
  })
  names(sort(score, decreasing = TRUE))[1]
})
out$observed_programme <- programme

out$verdict <- with(out, ifelse(
  grepl("Mono", rename_label) & observed_programme %in% c("Mono", "FcMono"), "rename vector consistent",
  ifelse(grepl("CD4", rename_label) & observed_programme == "CD4", "rename vector consistent",
  ifelse(grepl("CD8", rename_label) & observed_programme == "CD8", "rename vector consistent",
  ifelse(grepl("NK", rename_label) & observed_programme == "NK", "rename vector consistent",
  ifelse(grepl("B", rename_label) & observed_programme == "B", "rename vector consistent",
  ifelse(grepl("DC", rename_label) & observed_programme == "DC", "rename vector consistent",
  ifelse(grepl("Platelet", rename_label) & observed_programme == "Plat", "rename vector consistent",
         "CHECK MANUALLY"))))))))

write.csv(out, args[[3]], row.names = FALSE)
print(out[, c("cluster", "n_cells", "table_label", "rename_label",
              "observed_programme", "verdict")], row.names = FALSE)
cat("\ntop markers per cluster:\n")
print(out[, c("cluster", "top5_markers")], row.names = FALSE)
