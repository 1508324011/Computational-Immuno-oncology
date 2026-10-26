# @@SOURCE class/scRNA-seq/scRNA_seq_analyze_totorial_2026.html (R cells, verbatim)
# @@NOTE only change: tutorial placeholder "yourfolder/hg19" -> absolute run path
# @@NOTE cell ids = html cell numbers (1-based, matching the nbconvert dump)
#@@CELL 11
library(Seurat)
library(dplyr)
library(ggplot2)

data_dir <- "/lustre1/user/bjx131_pkuhpc/scrna_run_20261026/hg19"
sc_data <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "scRNAseq")

#@@CELL 13
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# QC visualize
p <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
ggsave("violin_mt.png", plot = p, width = 10, height = 5, dpi = 300)

# Setting Thresholds
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#@@CELL 16
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

#@@CELL 18
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
#HVG visualization
p <- VariableFeaturePlot(seurat_obj)
top_features <- head(VariableFeatures(seurat_obj), 20)
plot2 <- LabelPoints(plot = p, points = top_features, repel = TRUE)
p <- p + plot2
ggsave("VariableFeaturePlot.png", plot = p, width = 10, height = 5, dpi = 300, bg = "white")

#@@CELL 21
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

#@@CELL 23
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# PCA visualization
p<- VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
ggsave("PCA_VizDimLoadings.png", plot = p, width = 10, height = 5, dpi = 300, bg = "white")
p <- ElbowPlot(seurat_obj)
ggsave("PCA.png", plot = p, width = 10, height = 5, dpi = 300, bg = "white")
p <- DimPlot(seurat_obj, reduction = "pca") + NoLegend()
ggsave("PCA_dimplot.png", plot = p, width = 10, height = 5, dpi = 300, bg = "white")
#p<- DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
png("PCA_heatmap_ALL.png", width = 10, height = 10, units = "in", res = 300)
p<- PCHeatmap(seurat_obj, dims = 1:6, cells = 100, balanced = TRUE, ncol = 3)
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
seurat_obj.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

top10 <- seurat_obj.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

p<-DoHeatmap(seurat_obj, features = top10$gene) + NoLegend() +
  theme(
    axis.text.y = element_text(size = 6)
  )
ggsave("Marker_Heatmap.png", plot = p, width = 10, height = 10, dpi = 300)

#@@CELL 39
#Cluster ID 	Markers 	Cell Type
#0 	IL7R, CCR7 	Naive CD4+ T
#1 	CD14, LYZ 	CD14+ Mono
#2 	IL7R, S100A4 	Memory CD4+
#3 	MS4A1 	B
#4 	CD8A 	CD8+ T
#5 	FCGR3A, MS4A7 	FCGR3A+ Mono
#6 	GNLY, NKG7 	NK
#7 	FCER1A, CST3 	DC
#8 	PPBP 	Platelet

#@@CELL 40
p<- FeaturePlot(seurat_obj, features = c("MS4A1", "GNLY","NKG7", "CD3E","IL7R","CD14", "FCER1A","CST3", "FCGR3A", "LYZ", "PPBP","CD8A"))
ggsave("marker_Features_2.png", plot = p, width = 20, height = 10, dpi = 300)

#@@CELL 43
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T","CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
# 重新可视化带注释的UMAP图dpi = 300)
p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("cell annotations.png", plot = p, width = 10, height = 5, dpi = 300)

#@@CELL 46
# DEG for  Naive CD4 T
naive_cd4_markers <- FindMarkers(seurat_obj, ident.1 = "Naive CD4 T", min.pct = 0.25)

# p-value < 0.05）
naive_cd4_genes <- rownames(naive_cd4_markers[naive_cd4_markers$p_val < 0.05, ])

# Memory CD4 T vs Naive CD4 T
memory_vs_naive_markers <- FindMarkers(seurat_obj, ident.1 = "Memory CD4 T", ident.2 = "Naive CD4 T", min.pct = 0.25)

memory_vs_naive_genes <- rownames(memory_vs_naive_markers[memory_vs_naive_markers$p_val < 0.05, ])

df <- naive_cd4_markers
df$gene <- rownames(df)
df$significance <- ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.5,
                          "Significant", "Not significant")

#top Significant genes
top_genes <- df %>%
  dplyr::filter(significance == "Significant") %>%
  dplyr::arrange(p_val_adj) %>%
  head(15)  #

p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano plot: Naive CD4 T markers",
       x = "avg_log2FC", y = "-log10(adj. p-value)")

ggsave("Naive_Memory_CD4_DEG.png", plot = p, width = 6, height = 5, dpi = 300)
#Geneset enrichment?

saveRDS(seurat_obj, file = "./pbmc3k_final.rds")

#@@CELL 48
# Load the required libraries
library(Seurat)

# Load your data
data_pbmc_multi <- load("/lustre1/user/bjx131_pkuhpc/scrna_run_20261026/pbmc_stim.RData")

# Assuming ctrl.sparse and stim.sparse are the loaded objects
ctrl <- get("ctrl.sparse")
stim <- get("stim.sparse")

# Create Seurat objects
ctrl_seurat <- CreateSeuratObject(counts = ctrl)
stim_seurat <- CreateSeuratObject(counts = stim)
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
