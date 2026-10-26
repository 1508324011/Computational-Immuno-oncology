# @@CORRECTED version of the tutorial's R cells (Part 1-7).
# @@Every change is tagged  # --- FIX Rn (SC-xx): reason ---  and is the minimal edit
# @@needed to make the cell reproducible / auditable. Everything else is verbatim.
# @@Extras that the tutorial only implies are marked  # --- ADD ... ---  (tables, metrics).
#@@CELL 0
# --- CONFIG (SC-19): CreateSeuratObject() is called with its default min.cells = 0 in the
#     tutorial, i.e. all 32,738 features survive (19k of them are all-zero after QC), while
#     scanpy's filter_genes(min_cells = 3) keeps 13,714. This switch keeps the tutorial value
#     by default; MIN_CELLS=3 reproduces the scanpy gene universe. The two are NOT equivalent:
#     measured on this dataset, 30/2000 HVGs and 9-cluster boundaries move (see doc SC-19). ---
MIN_CELLS <- as.integer(Sys.getenv("MIN_CELLS", "0"))

#@@CELL 11
library(Seurat)
library(dplyr)
library(ggplot2)
# --- FIX R0 (SC-20): the DEG cell calls geom_text_repel() but the tutorial never loads
#     ggrepel (it only loads Seurat/dplyr/ggplot2), so that whole cell -- FindMarkers
#     results, volcano plot AND saveRDS() -- dies with 'could not find function'. ---
library(ggrepel)

data_dir <- "/lustre1/user/bjx131_pkuhpc/scrna_run_20261026/hg19"
sc_data <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "scRNAseq", min.cells = MIN_CELLS)
cat(sprintf("after CreateSeuratObject(min.cells=%d): %d genes x %d cells\n",
            MIN_CELLS, nrow(seurat_obj), ncol(seurat_obj)))

#@@CELL 13
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# QC visualize
p <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
ggsave("violin_mt.png", plot = p, width = 10, height = 5, dpi = 300)

n_before <- ncol(seurat_obj)
# Setting Thresholds
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cat(sprintf("QC: %d -> %d cells, %d genes\n", n_before, ncol(seurat_obj), nrow(seurat_obj)))

#@@CELL 16
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

#@@CELL 18
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
#HVG visualization
p <- VariableFeaturePlot(seurat_obj)
top_features <- head(VariableFeatures(seurat_obj), 20)
# --- FIX R2 (SC-08): LabelPoints() already returns a complete ggplot built on p;
#     the tutorial then adds the whole plot again (p + plot2), duplicating every layer. ---
p <- LabelPoints(plot = p, points = top_features, repel = TRUE)
ggsave("VariableFeaturePlot.png", plot = p, width = 10, height = 5, dpi = 300, bg = "white")

#@@CELL 21
# --- NOTE (SC-01): scaling all 32,738 features is what makes DoHeatmap() able to draw
#     arbitrary marker genes later (slot = "scale.data"). It is a deliberate trade-off,
#     not an error; measured cost on this machine is reported in the doc. ---
t_scale <- system.time(seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj)))
cat(sprintf("ScaleData on %d features: %.2f s, objects %.2f GB\n",
            nrow(seurat_obj), t_scale[["elapsed"]],
            sum(sapply(ls(), function(o) object.size(get(o)))) / 2^30))

#@@CELL 23
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# PCA visualization
p<- VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
ggsave("PCA_VizDimLoadings.png", plot = p, width = 10, height = 5, dpi = 300, bg = "white")
p <- ElbowPlot(seurat_obj)
ggsave("PCA.png", plot = p, width = 10, height = 5, dpi = 300, bg = "white")
p <- DimPlot(seurat_obj, reduction = "pca") + NoLegend()
ggsave("PCA_dimplot.png", plot = p, width = 10, height = 5, dpi = 300, bg = "white")
# --- FIX R3 (SC-07): PCHeatmap() is a deprecated alias; DimHeatmap() is the supported call. ---
png("PCA_heatmap_ALL.png", width = 10, height = 10, units = "in", res = 300)
p <- DimHeatmap(seurat_obj, dims = 1:6, cells = 100, balanced = TRUE, ncol = 3)
dev.off()

#@@CELL 27
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
# Clustering
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

#@@CELL 29
# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
p <- DimPlot(seurat_obj, reduction = "umap")
ggsave("UMAP.png", plot = p, width = 10, height = 5, dpi = 300)

#@@CELL 33
p<- FeaturePlot(seurat_obj, features = c("PTPRC","ITGAM", "CD3D","CD19", "CD14", "CD4", "CD8A", "CD68","NCAM1"))
ggsave("marker_Features.png", plot = p, width = 10, height = 10, dpi = 300)

#@@CELL 36
seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
# --- FIX R4 (SC-03): the tutorial computes the p_adj/log2FC-filtered table and throws it
#     away (unassigned pipe); top10 then comes from the *unfiltered* table. Keep the
#     filtered table, and fall back to the unfiltered one for clusters with <10 hits. ---
markers_sig <- seurat_obj.markers %>%
    dplyr::filter(p_val_adj < 0.05, avg_log2FC > 1)

top10 <- markers_sig %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)
cat(sprintf("FindAllMarkers: %d rows total, %d pass p_adj<0.05 & log2FC>1, %d used for heatmap\n",
            nrow(seurat_obj.markers), nrow(markers_sig), nrow(top10)))
print(table(markers_sig$cluster))

p<-DoHeatmap(seurat_obj, features = top10$gene) + NoLegend() +
  theme(
    axis.text.y = element_text(size = 6)
  )
ggsave("Marker_Heatmap.png", plot = p, width = 10, height = 10, dpi = 300)
# --- ADD (SC-13): the tutorial never writes the marker table out; export it. ---
write.csv(seurat_obj.markers, "markers_all_clusters.csv", row.names = FALSE)
write.csv(top10, "markers_top10_per_cluster.csv", row.names = FALSE)

#@@CELL 39
# --- ADD (SC-02 support): cluster -> canonical marker -> label audit table. ---
audit <- data.frame(
  cluster = levels(Idents(seurat_obj)),
  tutorial_marker = c("IL7R, CCR7","CD14, LYZ","IL7R, S100A4","MS4A1","CD8A",
                      "FCGR3A, MS4A7","GNLY, NKG7","FCER1A, CST3","PPBP"),
  tutorial_label  = c("Naive CD4 T","CD14+ Mono","Memory CD4 T","B","CD8 T",
                      "FCGR3A+ Mono","NK","DC","Platelet"),
  stringsAsFactors = FALSE)
audit$top_markers <- sapply(audit$cluster, function(cl)
  paste(head(markers_sig$gene[markers_sig$cluster == cl], 5), collapse = ", "))
print(audit[, c("cluster", "tutorial_label", "top_markers")])
write.csv(audit, "cluster_label_audit.csv", row.names = FALSE)

#@@CELL 40
p<- FeaturePlot(seurat_obj, features = c("MS4A1", "GNLY","NKG7", "CD3E","IL7R","CD14", "FCER1A","CST3", "FCGR3A", "LYZ", "PPBP","CD8A"))
ggsave("marker_Features_2.png", plot = p, width = 20, height = 10, dpi = 300)

#@@CELL 43
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T","CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
# --- FIX R5 (SC-02): RenameIdents() silently mislabels everything if the cluster count
#     drifts (a different resolution, Seurat version or QC threshold is enough). Check it. ---
# --- FIX R5 (SC-02): RenameIdents() silently mislabels everything if the cluster count
#     drifts (a different resolution, Seurat version or QC threshold is enough). Check the
#     count, print the audit table, and fall back to neutral labels instead of lying. ---
n_cl <- nlevels(seurat_obj)
if (n_cl == length(new.cluster.ids)) {
  names(new.cluster.ids) <- levels(seurat_obj)
  seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
} else {
  warning(sprintf("%d clusters but %d reference names; keeping neutral labels. See cluster_label_audit.csv",
                  n_cl, length(new.cluster.ids)), call. = FALSE)
  neutral <- setNames(paste0("C", seq_len(n_cl) - 1L), levels(seurat_obj))
  seurat_obj <- RenameIdents(seurat_obj, neutral)
}
print(table(Idents(seurat_obj)))
# 重新可视化带注释的UMAP图
p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("cell annotations.png", plot = p, width = 10, height = 5, dpi = 300)

#@@CELL 46
# DEG for  Naive CD4 T
# --- FIX R5b (SC-02): if the adaptive guard above had to switch to neutral labels,
#     these two FindMarkers() calls would both fail; point them at whatever label exists
#     and say so, instead of silently producing nothing. ---
if (!"Naive CD4 T" %in% levels(Idents(seurat_obj))) {
  warning("label 'Naive CD4 T' absent (cluster drift); using the first cluster as ident.1", call. = FALSE)
}
id_naive <- if ("Naive CD4 T" %in% levels(Idents(seurat_obj))) "Naive CD4 T" else levels(Idents(seurat_obj))[1]
naive_cd4_markers <- FindMarkers(seurat_obj, ident.1 = id_naive, min.pct = 0.25)

# p-value < 0.05
naive_cd4_genes <- rownames(naive_cd4_markers[naive_cd4_markers$p_val_adj < 0.05, ])

# Memory CD4 T vs Naive CD4 T
id_memory <- if ("Memory CD4 T" %in% levels(Idents(seurat_obj))) "Memory CD4 T" else levels(Idents(seurat_obj))[2]
memory_vs_naive_markers <- FindMarkers(seurat_obj, ident.1 = id_memory, ident.2 = id_naive, min.pct = 0.25)

memory_vs_naive_genes <- rownames(memory_vs_naive_markers[memory_vs_naive_markers$p_val_adj < 0.05, ])

df <- naive_cd4_markers
df$gene <- rownames(df)
df$significance <- ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.5,
                          "Significant", "Not significant")

#top Significant genes
top_genes <- df %>%
  dplyr::filter(significance == "Significant") %>%
  dplyr::arrange(p_val_adj) %>%
  head(15)

p <- ggplot(df, aes(x = avg_log2FC, y = -log10(pmax(p_val_adj, .Machine$double.xmin)),
                    color = significance)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano plot: Naive CD4 T markers",
       x = "avg_log2FC", y = "-log10(adj. p-value)")

ggsave("Naive_Memory_CD4_DEG.png", plot = p, width = 6, height = 5, dpi = 300)
#Geneset enrichment?
# --- ADD (SC-04/SC-13): export the two DEG tables the tutorial only keeps in memory;
#     the -log10(p_adj) axis is also finite now (p_adj = 0 used to give Inf). ---
write.csv(cbind(gene = rownames(naive_cd4_markers), naive_cd4_markers),
          "DEG_NaiveCD4_vs_rest.csv", row.names = FALSE)
write.csv(cbind(gene = rownames(memory_vs_naive_markers), memory_vs_naive_markers),
          "DEG_MemoryCD4_vs_NaiveCD4.csv", row.names = FALSE)
cat(sprintf("Naive CD4 DEG: %d genes tested, %d significant (p_adj<0.05 & |log2FC|>0.5); %d significant by p_adj only\n",
            nrow(df), sum(df$significance == "Significant"), length(naive_cd4_genes)))
cat(sprintf("Memory vs Naive: %d genes tested, %d significant by p_adj\n",
            nrow(memory_vs_naive_markers), length(memory_vs_naive_genes)))
saveRDS(seurat_obj, file = "./pbmc3k_final.rds")

#@@CELL 48
# Load the required libraries
library(Seurat)

# --- FIX R6 (SC-05): load() returns the *names* of the objects it created, and dumps
#     them into the global environment; in a function or a knitted doc that silently
#     leaks variables. Load into a private environment and assert the shapes instead. ---
stim_env <- new.env()
load("/lustre1/user/bjx131_pkuhpc/scrna_run_20261026/pbmc_stim.RData", envir = stim_env)
ctrl <- get("ctrl.sparse", envir = stim_env)
stim <- get("stim.sparse", envir = stim_env)
cat(sprintf("ctrl.sparse: %d genes x %d cells; stim.sparse: %d genes x %d cells\n",
            nrow(ctrl), ncol(ctrl), nrow(stim), ncol(stim)))
stopifnot(identical(rownames(ctrl), rownames(stim)))

# --- FIX R7 (SC-06): the tutorial integrates raw counts with no QC at all, so its UMAPs
#     mix low-quality cells with real batch structure. Apply the same QC as Part 2. ---
qc1 <- CreateSeuratObject(counts = ctrl, project = "ctrl")
qc1[["percent.mt"]] <- PercentageFeatureSet(qc1, pattern = "^MT-")
n1 <- ncol(qc1)
qc1 <- subset(qc1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
qc2 <- CreateSeuratObject(counts = stim, project = "stim")
qc2[["percent.mt"]] <- PercentageFeatureSet(qc2, pattern = "^MT-")
n2 <- ncol(qc2)
qc2 <- subset(qc2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cat(sprintf("QC: ctrl %d -> %d cells, stim %d -> %d cells\n", n1, ncol(qc1), n2, ncol(qc2)))
ctrl_seurat <- qc1
stim_seurat <- qc2
ctrl_seurat$orig.ident <- "ctrl"
stim_seurat$orig.ident <- "stim"

# Merge the datasets
merged_seurat <- merge(ctrl_seurat, y = stim_seurat, add.cell.ids = c("ctrl", "stim"))
# Normalize the data
merged_seurat <- NormalizeData(merged_seurat)
# Find variable features
merged_seurat <- FindVariableFeatures(merged_seurat)
# Scale the data
merged_seurat <- ScaleData(merged_seurat)
# Perform PCA
merged_seurat <- RunPCA(merged_seurat)
# Run UMAP
merged_seurat <- RunUMAP(merged_seurat, dims = 1:10)
p <- DimPlot(merged_seurat, group.by = "orig.ident")  # Use orig.ident to distinguish between ctrl and stim
ggsave("UMAP_without_harmony.png", plot = p, width = 10, height = 5, dpi = 300)

# --- ADD (SC-06): quantify the batch structure instead of only looking at the UMAP.
#     cross_knn = share of each cell's 30 nearest neighbours that come from the other
#     sample; 0 = perfect separation, ~0.5 = fully mixed. ---
cross_knn <- function(emb, sample, k = 30) {
  nn <- FNN::get.knn(emb, k = k)$nn.index
  mean(sample[nn] != sample)
}
smp <- as.character(merged_seurat$orig.ident)
cross_pca <- cross_knn(Embeddings(merged_seurat, "pca")[, 1:10], smp)
cat(sprintf("cross-sample neighbour fraction, PCA space (no harmony): %.4f\n", cross_pca))

#@@CELL 51
#install.packages("harmony")
library(harmony)
# Run Harmony
merged_seurat <- RunHarmony(merged_seurat, group.by.vars = "orig.ident")
# Run UMAP with Harmony embeddings
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:10)
# Plot UMAP with Harmony
p <- DimPlot(merged_seurat, group.by = "orig.ident")
ggsave("UMAP_with_harmony.png", plot = p, width = 10, height = 5, dpi = 300)

# --- ADD (SC-06): same metric after harmony, plus a clustering-based sanity check. ---
cross_harmony <- cross_knn(Embeddings(merged_seurat, "harmony")[, 1:10], smp)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:10)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)
tab <- table(merged_seurat$seurat_clusters, merged_seurat$orig.ident)
share <- prop.table(tab, margin = 1)
cat(sprintf("cross-sample neighbour fraction, harmony space: %.4f (PCA: %.4f)\n",
            cross_harmony, cross_pca))
cat(sprintf("clusters after harmony: %d; clusters with a sample share >0.9: %d; mean |share-0.5| = %.3f\n",
            nrow(tab), sum(apply(share, 1, max) > 0.9), mean(abs(share[, 1] - 0.5))))
write.csv(as.data.frame.matrix(tab), "harmony_cluster_by_sample.csv")
# --- ADD: one tidy metrics file for the doc. ---
write.csv(data.frame(metric = c("cells_ctrl", "cells_stim", "cells_merged",
                                "cross_knn_pca", "cross_knn_harmony",
                                "clusters_harmony", "clusters_dominated_gt90pct"),
                     value = c(ncol(qc1), ncol(qc2), ncol(merged_seurat),
                               round(cross_pca, 4), round(cross_harmony, 4),
                               nrow(tab), sum(apply(share, 1, max) > 0.9))),
          "integration_metrics.csv", row.names = FALSE)
saveRDS(merged_seurat, file = "./pbmc_stim_merged_harmony.rds")
