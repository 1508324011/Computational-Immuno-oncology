#!/usr/bin/env Rscript
# marker_panel_audit.R -- mean log-normalised expression of a canonical PBMC marker panel
# per cluster, computed directly from the "data" layer of the delivered pbmc3k_final.rds
# (Seurat side). Same panel as marker_panel_audit.py so the two tables read side by side.
#
# usage: Rscript marker_panel_audit.R <pbmc3k_final.rds> <out_csv>
args <- commandArgs(trailingOnly = TRUE)
suppressMessages({library(Seurat); library(Matrix)})
seu <- readRDS(args[[1]])

panel <- list(
  "CD4 naive"    = c("IL7R", "CCR7", "LEF1", "SELL", "TCF7", "MAL"),
  "CD4 memory"   = c("IL32", "LTB", "CD2", "S100A4"),
  "CD8 T"        = c("CD8A", "CD8B", "GZMK", "CCL5", "CTSW"),
  "NK"           = c("GNLY", "GZMB", "NKG7", "KLRF1", "PRF1", "FGFBP2", "SPON2"),
  "B"            = c("MS4A1", "CD79A", "CD79B", "CD19", "TCL1A"),
  "CD14+ Mono"   = c("CD14", "LYZ", "S100A8", "S100A9", "FCN1"),
  "FCGR3A+ Mono" = c("FCGR3A", "MS4A7", "LST1"),
  "DC"           = c("FCER1A", "CST3", "CLEC10A", "SERPINF1"),
  "Platelet"     = c("PPBP", "PF4", "GP9"),
  "technical"    = c("RPS12", "RPL32", "RPS27", "PTPRC", "ACTB"))

all_genes <- unique(unlist(panel))
present <- all_genes[all_genes %in% rownames(seu)]
missing <- setdiff(all_genes, present)
cat("markers absent from this object:",
    if (length(missing)) paste(missing, collapse = ", ") else "none", "\n")

data <- LayerData(seu, assay = "RNA", layer = "data")        # log-normalised, genes x cells
idents <- as.character(Idents(seu))
groups <- levels(Idents(seu))

rows <- lapply(present, function(g) {
  v <- as.numeric(data[g, ])
  rec <- data.frame(gene = g,
                    programme = names(panel)[which(sapply(panel, function(p) g %in% p))[1]],
                    stringsAsFactors = FALSE)
  for (grp in groups) rec[[grp]] <- round(mean(v[idents == grp]), 4)
  rec
})
out <- do.call(rbind, rows)
write.csv(out, args[[2]], row.names = FALSE)

cat("cells per cluster:\n"); print(table(Idents(seu)))
cat("mean log-normalised expression (rows = marker, cols = cluster):\n")
print(out[, c("programme", "gene", groups)], row.names = FALSE)
