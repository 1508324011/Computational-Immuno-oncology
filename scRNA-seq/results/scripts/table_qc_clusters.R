#!/usr/bin/env Rscript
# table_qc_clusters.R -- build the two summary tables that the tutorial only ever
# prints to the console: QC accounting (01_qc/qc_summary.csv) and per-cluster sizes
# for both toolchains (03_clustering/cluster_sizes.csv).
#
# usage: Rscript table_qc_clusters.R <run_root> <results_dir>
args <- commandArgs(trailingOnly = TRUE)
run_root <- args[[1]]
res_dir  <- args[[2]]

suppressMessages({library(Seurat); library(dplyr); library(Matrix)})

## ---- 1. QC accounting, computed once from the raw matrix --------------------
sc  <- Read10X(file.path(run_root, "hg19"))
nFeature <- Matrix::colSums(sc > 0)
nCount   <- Matrix::colSums(sc)
mt       <- Matrix::colSums(sc[grepl("^MT-", rownames(sc)), ])
pct_mt   <- 100 * mt / nCount

genes_ge3 <- rownames(sc)[Matrix::rowSums(sc > 0) >= 3]

keep_seurat <- nFeature > 200 & nFeature < 2500 & pct_mt < 5
keep_scanpy <- nFeature >= 200 & nFeature < 2000 & pct_mt < 5

# The *actual* numbers come from the runs themselves (each is printed by the driver and
# recorded in results/09_run_logs/*_cell_log.txt); the recount below is an independent
# cross-check computed on the raw matrix, NOT the source of the reported values.
qc <- data.frame(
  tool = c("Seurat (tutorial Part 2-7, as written)",
           "Seurat (corrected run, MIN_CELLS = 3)",
           "scanpy (tutorial Part 8, n_genes < 2000)",
           "scanpy (corrected run, N_GENES_MAX = 2500)"),
  cell_rule = c("nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5",
                "nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5",
                "n_genes >= 200 & n_genes_by_counts < 2000 & pct_counts_mt < 5",
                "n_genes >= 200 & n_genes_by_counts < 2500 & pct_counts_mt < 5"),
  gene_rule = c("CreateSeuratObject default (min.cells = 0): 32738 features",
                "CreateSeuratObject(min.cells = 3): 13714 features",
                "filter_genes(min_cells = 3): 13714 features",
                "filter_genes(min_cells = 3): 13714 features"),
  cells_in = ncol(sc),
  cells_kept = c(2638L, 2638L, 2635L, 2638L),
  genes_kept = c(32738L, 13714L, 13714L, 13714L),
  evidence = c("09_run_logs/R_tutorial_cell_log.txt (cell 13)",
               "09_run_logs/R_corrected_mincells0_cell_log.txt + min.cells=3 pass",
               "09_run_logs/Py_asis_cell_log.txt (cell 61)",
               "09_run_logs/Py_fixed_2500_cell_log.txt (cell 61)"),
  stringsAsFactors = FALSE)

# Independent recount on the raw matrix (no tool involved) - documented for the
# SC-10 comparison: the same nominal 5% mitochondrial cut keeps 2634 cells here but
# 2635 in the scanpy run, because scanpy recomputes pct_counts_mt on the *gene-filtered*
# matrix (13714 genes) while Seurat uses all 32738 features.
rederive <- data.frame(
  rule = c("nFeature > 200 & nFeature < 2500 & pct_mt < 5  (Seurat rule, all 32738 genes)",
           "nFeature >= 200 & nFeature < 2000 & pct_mt < 5   (scanpy rule, all 32738 genes)",
           "nFeature > 200 & nFeature < 2000 & pct_mt < 5    (Seurat rule, threshold of scanpy)",
           "nFeature >= 200 & nFeature < 2500 & pct_mt < 5   (scanpy rule, threshold of Seurat)"),
  cells_kept = c(sum(keep_seurat), sum(keep_scanpy),
                 sum(nFeature > 200 & nFeature < 2000 & pct_mt < 5),
                 sum(nFeature >= 200 & nFeature < 2500 & pct_mt < 5)),
  cells_dropped_lowFeature  = c(sum(nFeature <= 200), sum(nFeature < 200),
                                sum(nFeature <= 200), sum(nFeature < 200)),
  cells_dropped_highFeature = c(sum(nFeature >= 2500), sum(nFeature >= 2000),
                                sum(nFeature >= 2000), sum(nFeature >= 2500)),
  cells_dropped_by_mt_only  = c(sum(nFeature > 200 & nFeature < 2500 & pct_mt >= 5),
                                sum(nFeature >= 200 & nFeature < 2000 & pct_mt >= 5),
                                sum(nFeature > 200 & nFeature < 2000 & pct_mt >= 5),
                                sum(nFeature >= 200 & nFeature < 2500 & pct_mt >= 5)))

mt_summary <- data.frame(
  metric = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  all_cells_median = c(median(nFeature), median(nCount), median(pct_mt)),
  all_cells_min = c(min(nFeature), min(nCount), min(pct_mt)),
  all_cells_max = c(max(nFeature), max(nCount), max(pct_mt)),
  kept_median = c(median(nFeature[keep_seurat]), median(nCount[keep_seurat]),
                  median(pct_mt[keep_seurat])),
  kept_min = c(min(nFeature[keep_seurat]), min(nCount[keep_seurat]), min(pct_mt[keep_seurat])),
  kept_max = c(max(nFeature[keep_seurat]), max(nCount[keep_seurat]), max(pct_mt[keep_seurat])))

write.csv(qc, file.path(res_dir, "01_qc", "qc_summary.csv"), row.names = FALSE)
write.csv(rederive, file.path(res_dir, "01_qc", "qc_rederivation_rawmatrix.csv"), row.names = FALSE)
write.csv(mt_summary, file.path(res_dir, "01_qc", "qc_metric_distribution.csv"), row.names = FALSE)

## ---- 2. cluster sizes -------------------------------------------------------
seu <- readRDS(file.path(run_root, "rr_fixed", "pbmc3k_final.rds"))
seu_tab <- as.data.frame(table(cluster = Idents(seu)), responseName = "n_cells")
seu_tab$tool <- "Seurat (resolution 0.5, 10 PCs)"
seu_tab$share <- round(seu_tab$n_cells / sum(seu_tab$n_cells), 4)

py_obs <- read.csv(file.path(run_root, "py_fixed", "pbmc3k_obs.csv"), row.names = 1)
py_tab <- as.data.frame(table(cluster = py_obs$leiden), responseName = "n_cells")
py_tab$tool <- "scanpy (leiden 0.9, 40 PCs, n_genes<2000)"
py_tab$share <- round(py_tab$n_cells / sum(py_tab$n_cells), 4)

py_obs2 <- read.csv(file.path(run_root, "py_align2500", "pbmc3k_obs.csv"), row.names = 1)
py_tab2 <- as.data.frame(table(cluster = py_obs2$leiden), responseName = "n_cells")
py_tab2$tool <- "scanpy (leiden 0.9, 40 PCs, n_genes<2500)"
py_tab2$share <- round(py_tab2$n_cells / sum(py_tab2$n_cells), 4)

out <- rbind(seu_tab[, c("tool", "cluster", "n_cells", "share")],
             py_tab[, c("tool", "cluster", "n_cells", "share")],
             py_tab2[, c("tool", "cluster", "n_cells", "share")])
write.csv(out, file.path(res_dir, "03_clustering", "cluster_sizes.csv"), row.names = FALSE)

cat("--- qc_summary.csv (from the runs) ---\n"); print(qc[, c("tool", "cells_in", "cells_kept", "genes_kept")])
cat("--- qc_rederivation_rawmatrix.csv ---\n"); print(rederive[, c("rule", "cells_kept")])
cat("--- Seurat clusters ---\n"); print(seu_tab)
cat("--- scanpy clusters (n_genes<2000) ---\n"); print(py_tab)
cat("--- scanpy clusters (n_genes<2500) ---\n"); print(py_tab2)
