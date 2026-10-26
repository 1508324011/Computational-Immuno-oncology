#!/usr/bin/env Rscript
# cluster_correspondence.R -- dump the Seurat per-cell annotation (used against the
# scanpy h5ad obs by cluster_correspondence.py to build the cross-tool contingency table).
# usage: Rscript cluster_correspondence.R <pbmc3k_final.rds> <out_csv>
suppressMessages(library(Seurat))
seu <- readRDS(commandArgs(trailingOnly = TRUE)[[1]])
out <- data.frame(barcode = colnames(seu),
                  seurat_cluster = as.character(seu$seurat_clusters),
                  seurat_label = as.character(Idents(seu)),
                  nFeature_RNA = seu$nFeature_RNA,
                  percent_mt = round(seu$percent.mt, 4),
                  stringsAsFactors = FALSE)
write.csv(out, commandArgs(trailingOnly = TRUE)[[2]], row.names = FALSE)
cat(sprintf("wrote %d cells\n", nrow(out)))
