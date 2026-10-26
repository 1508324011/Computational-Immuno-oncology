#!/usr/bin/env python3
"""Insert the figures + analysis into the scRNA-seq corrected editions (EN and ZH).

Each block is anchored to a unique sentence that already exists in the document; the
figure(s) land right after the paragraph that ends with it, so the prose keeps its flow
and every figure sits next to the result it explains. Figure numbers are identical in
the two languages (EN is authoritative; ZH is generated from the same block list).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def insert_after_para(text, anchor, block, label):
    n = text.count(anchor)
    if n != 1:
        raise SystemExit(f"[{label}] anchor found {n}x: {anchor[:70]!r}")
    i = text.index(anchor)
    j = text.index("\n\n", i)
    return text[: j + 2] + block.strip("\n") + "\n\n" + text[j + 2 :]


# ---------------------------------------------------------------- EN blocks ----
EN = [
    (
        "It is recorded in `results/01_qc/qc_rederivation_rawmatrix.csv`.",
        r"""
![Figure 1 — Seurat QC distributions for all 2,700 cells: genes per cell, UMIs per cell, and mitochondrial percentage, with the tutorial's three thresholds.](./figures/01_violin_mt.png){width=100%}

Read the violins against the rule `nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5`. The nFeature distribution has the healthy shape of a shallow 2016 run — median 817 genes, a long left tail that never reaches 200 (hence **0** cells lost on that side) and a thin right shoulder above 2,500 (5 cells, the only candidates for doublets this pipeline ever removes). The mitochondrial violin is nearly flat below 5% with a small right shoulder; that shoulder is the 57 cells removed. There is no second mode at high UMI, which is why nothing more elaborate than these three numbers is warranted here.

![Figure 2 — The same QC metrics computed by scanpy (`calculate_qc_metrics`), after `filter_cells(min_genes=200)` and `filter_genes(min_cells=3)`.](./figures/02_scanpy_violin_qc.png){width=100%}

The Python panel is not a copy of the R panel and the difference matters: these distributions are computed on the **13,714-gene** filtered matrix, while Seurat's are computed on all 32,738 features. The mitochondrial percentage is a ratio, so removing 19,024 all-zero genes moves its denominator for the borderline cells — this is exactly the one-cell discrepancy that makes Seurat keep 2,638 and scanpy 2,635 (SC-10).

![Figure 3 — `total_counts` versus `n_genes_by_counts` for every cell.](./figures/03_scanpy_scatter_genes.png){width=100%}

This pair is the doublet test the tutorial never names. The main cloud is almost perfectly linear — every extra UMI buys genes at a constant rate — and doublets sit *above* that line (two cells' worth of genes for the same UMI count). Here the cloud is tight and there is no separated upper branch, so 2,500 genes as the upper cut is doing the doublet job by proxy.

![Figure 4 — `total_counts` versus `pct_counts_mt`.](./figures/04_scanpy_scatter_mt.png){width=100%}

The opposite reading applies: dying cells lose cytoplasmic mRNA while mitochondrial transcripts stay inside the still-intact organelle, so they drift up and to the left. The vertical cut at 5% is drawn through the sparse upper cloud rather than through a gap, which is the honest description of any mitochondrial threshold — it is a convention, not a physical boundary.
""",
    ),
    (
        "SC-08 is here: `p <- p + plot2` renders **two panels**",
        r"""
![Figure 5 — `VariableFeaturePlot.png` as the tutorial writes it: two panels, because `p <- p + plot2` prints the old plot beside the labelled one.](./figures/05_variable_features_aswritten.png){width=100%}

Both halves show the same 32,738 genes with the 2,000 HVGs in red; the left half has no gene labels, the right half labels twenty. This is the cleanest example in the module of code that runs, produces a usable-looking file, and still does something other than what its author meant — and the next figure shows that the teacher's own reference figure carries the same artefact.

![Figure 6 — The reference `VariableFeaturePlot.png` shipped with the course: the same two-panel layout, i.e. the artefact was shipped as the expected output.](./figures/ref_variable_featureplot_2panel.png){width=100%}

Comparing this with Figure 5 is how "reproduction" should be judged: our as-written run agrees with the reference at a structural correlation of 0.9751 (identical size, identical point cloud, identical labels), because both were produced by the same two-panel bug. A student who "fixes" the plot and then compares to the reference would conclude, wrongly, that their analysis had diverged.

![Figure 7 — The corrected single panel (`VariableFeaturePlot_corrected.png`).](./figures/06_variable_features_corrected.png){width=100%}

The corrected version draws one scatter with the twenty labels attached by `LabelPoints(repel = TRUE)`. Its correlation with the reference drops to 0.0713 *by construction*: the panels are different pictures of the same data. The labels themselves are the sanity check — PPBP, LYZ, S100A9, IGLL5, FTL, GNLY, PF4, FTH1, GNG11, S100A8, FCER1A, HLA-DRA, CD74, CLU, NKG7, CST3, GZMB, CCL4, C1QA, HLA-DPB1 are all leukocyte genes, which is what `vst` should select on PBMCs.

![Figure 8 — scanpy's HVG selection on the other half: normalised dispersion against mean expression, 1,833 genes selected.](./figures/07_scanpy_hvg.png){width=100%}

The two methods are not interchangeable and this figure is where you can see why: `vst` standardises a loess fit of variance on log-mean and can promote a gene with modest mean (PF4, PPBP), while scanpy's `min_mean 0.0125 / max_mean 3 / min_disp 0.5` cuts on the dispersion of the normalised, binned data and returns 1,833 genes from a 13,714-gene universe. Neither number is wrong; they answer "which genes vary more than their expression level predicts" with different variance models.
""",
    ),
    (
        "but the elbow is drawn for 20 PCs and the choice of 10 is never justified in the text.",
        r"""
![Figure 9 — `PCA_VizDimLoadings` for PC1 and PC2: the genes with the largest positive and negative loadings in each component.](./figures/08_pca_loadings.png){width=100%}

Loadings are the interpretable half of a PCA. PC1 is dominated by the monocyte-versus-lymphocyte axis (LYZ, S100A8, S100A9, CD14 on one side; IL7R, CD3D, LTB on the other), PC2 separates FCGR3A+ monocytes and NK cells. If a component had no coherent gene list, that would be the signal to stop using it — here the first components are biologically legible, which is a weak but real argument for `dims = 1:10`.

![Figure 10 — `PCA.png`: the standard-deviation elbow for 20 components.](./figures/09_pca_elbow.png){width=100%}

The elbow bends around PC 7-10 and then flattens, which is the (informal) justification for the `dims = 1:10` used by `FindNeighbors`, `FindClusters` and `RunUMAP`. Note what the plot does *not* do: it shows 20 components although `RunPCA` computed 50, and it offers no criterion — the choice of 10 is a visual convention, and SC-07 records that the tutorial never reconciles the 20 drawn with the 10 used.

![Figure 11 — `PCA_dimplot`: every cell placed on PC1 and PC2, coloured by Seurat cluster.](./figures/10_pca_dimplot.png){width=100%}

Before any UMAP exists, the clusters already separate on two components: the monocyte arm, the B/NK arm and the T-cell arm occupy different quadrants. What PCA cannot do is separate naive from memory CD4 T — they overlap in the middle — and that is precisely the split that needs ten components and a graph.

![Figure 12 — `PCA_heatmap_ALL`: the top 6 components, each drawn for 100 cells taken evenly from the extremes of that component.](./figures/11_pca_heatmap.png){width=100%}

Each row is a gene, each column a cell; the block structure is what a "good" component looks like — a set of co-regulated genes switching on together across the sampled cells. This figure is also our strongest reproduction evidence: its structural correlation with the teacher's reference PNG is **1.0000**, i.e. the same PCA, the same balanced sampling and the same plotting path reproduce the file essentially bit for bit.

![Figure 13 — scanpy's PCA (PC1/PC2, coloured by CST3) computed on 1,833 HVGs after `sc.pp.scale(max_value = 10)`.](./figures/12_scanpy_pca.png){width=100%}

Comparing this with Figure 11 is instructive: the same biology appears, but the axes are oriented and scaled differently, because scanpy scaled each gene to unit variance (clipping at 10) over 1,833 genes while Seurat scaled 32,738 and then projected the 2,000 HVGs. Neither orientation is canonical — components are only defined up to sign and rotation — so a "different-looking PCA" between the two halves is expected and not a bug.

![Figure 14 — `pca_variance_ratio` — the figure the tutorial's own code never produced.](./figures/13_scanpy_pca_ratio.png){width=100%}

This panel is missing from the 28 reference figures for a mundane reason: the cell that draws it contains the doubled comma (SC-09), and Python compiles a cell as a unit, so the `SyntaxError` also deleted the PCA, the neighbour graph and the UMAP that were written in the same cell. Ours comes from the corrected pass; the cumulative curve reaching ~95% only around PC 40 is why the tutorial's `n_pcs = 40` for the neighbour graph is defensible on this dataset.
""",
    ),
    (
        "at `resolution = 0.5` it survives, at 0.3 it would not.",
        r"""
![Figure 15 — The Seurat 9-cluster UMAP of 2,638 cells (`UMAP.png`).](./figures/14_umap_seurat.png){width=100%}

The topology is the PBMC 3k textbook picture: one large T-cell lobe (clusters 0, 1 and 4), a monocyte arm (2 and 5) that is separated by a long bridge of empty space, a compact B island (3), and the two small satellites — NK (6) and the immune-poor corner holding DC (7) and platelets (8). What to check first is not the colours but the two tiny clusters: 28 and 14 cells out of 2,638 are 1.1% and 0.5%, they are visually separate from everything, and their marker panels (FCER1A/CLEC10A vs PPBP/PF4) are unambiguous. The visual caveat: UMAP distances between clusters are not quantitative, so "cluster 2 is far from cluster 0" must never be read as "monocytes are more different from CD4 T than B cells are".

![Figure 16 — The same data clustered by scanpy: leiden at resolution 0.9 plus twelve marker panels (`umapscanpy_markers.png`).](./figures/15_scanpy_markers_grid.png){width=100%}

One figure that replaces nine: the first panel is the leiden solution (8 groups at the tutorial's threshold), and the rest are the marker genes overlaid on the same embedding. Read it as a consistency check between two independent tools — B cells light up MS4A1/CD79A in exactly one group, monocytes light up CD14/LYZ/FCGR3A, platelets light up PPBP — and note where it disagrees with Seurat: NK markers (GNLY, NKG7, KLRF1) are spread over two groups here, because the 40-PC neighbour graph merges CD8 T and NK into one neighbourhood that this resolution then splits on a different boundary.
""",
    ),
    (
        "They are correct *by luck*:",
        r"""
![Figure 17 — `marker_Features` with the nine lineage markers the tutorial names in its annotation table.](./figures/16_featureplot_markers.png){width=100%}

This is the figure a human actually annotates from, so it is worth reading methodically. Go marker by marker and note which cluster lights up *and which ones stay dark*: CD14 and LYZ on cluster 2 only, MS4A1 on cluster 3 only, GNLY/NKG7 on cluster 6, CD3E/IL7R on 0 and 1, FCER1A/CST3 on 7, PPBP on 8. Two things follow. First, seven of nine clusters are unambiguous on this panel alone. Second, CD3E/IL7R cannot separate clusters 0 and 1 — both are T cells — which is precisely why the tutorial needs the second panel (Figure 18) and the DEG analysis (Figure 23) to tell naive from memory.

![Figure 18 — `marker_Features_2`: twelve additional markers on the same embedding.](./figures/17_featureplot_markers2.png){width=100%}

The extra markers resolve what the first panel could not. IGFBP7/TNFRSF4 style markers aside, this panel is where FCGR3A+ monocytes (cluster 5: FCGR3A, MS4A7, LST1) get separated from classical monocytes (cluster 2), where CCR7 and LEF1 mark cluster 1 as the naive CD4 compartment and S100A4 marks cluster 0 as the memory one, and where KLRB1/KLRD1 confirm the NK identity of cluster 6. The practical lesson: annotation needs *pairs* of markers per lineage, because any single marker is a proxy and several are shared (LYZ is expressed at 0.96 in naive CD4 T cells as well as 5.09 in monocytes).

![Figure 19 — The final annotated UMAP of the Seurat object (`cell annotations.png`).](./figures/20_umap_annotated_seurat.png){width=100%}

Nine labels drawn on the data, one per cluster, matching the audit table in §9.1. Use this figure as the reference for everything downstream: the monocyte arm is `CD14+ Mono` (cluster 2, 471 cells) and `FCGR3A+ Mono` (cluster 5, 175), the upper T lobe is `Naive CD4 T` (0, 657) plus `Memory CD4 T` (1, 503), and the far cluster 4 is `CD8 T` (307). Any claim in a later section that touches a compartment can be checked against this colour, which is why the mislabelling risk in §13.2 (cluster 1/2 swapped after a 30-HVG perturbation) would have looked completely normal in this picture.
""",
    ),
    (
        "the resulting heatmap differs from the reference at 0.9865.",
        r"""
![Figure 20 — `Marker_Heatmap.png` as written: the top 10 genes per cluster by `avg_log2FC`, taken from the **unfiltered** marker table.](./figures/18_heatmap_aswritten.png){width=100%}

The block structure is reassuring — each cluster has a bright gene set that is dark everywhere else — but this figure is the visual form of SC-03. Because the filter was computed and discarded (`dplyr::filter()` inside a pipe whose value nobody assigned), `slice_max(avg_log2FC)` ranks by effect size alone and can, and does, put a gene with `p_val_adj = 0.2` next to a gene with `p_val_adj = 1e-40`. A reader cannot see that from the colours: the values plotted are normalised counts, so a non-significant gene still looks convincing if it is high in a few cells.

![Figure 21 — The corrected heatmap: 90 genes chosen from the 3,483 markers that pass `p_val_adj < 0.05 & avg_log2FC > 1`.](./figures/19_heatmap_corrected.png){width=100%}

Same nine blocks, same colours, different gene list. The interesting comparison is which blocks *stay* sharp: B cells (CD79A/MS4A1/TCL1A), monocytes (LYZ/S100A8/S100A9) and platelets (PPBP/PF4) look the same in both versions because their effects are overwhelming; the CD4 blocks change most, because "naive vs memory" is a subtle contrast (log2FC of 2-4 rather than 20) and that is exactly where a significance filter earns its keep. Correlation with the reference: 0.7586 — the visible gap is the removed genes, not a different clustering.

![Figure 22 — The reference `Marker_Heatmap.png` shipped with the course: 3,000 × 6,000 px, i.e. 10 × 20 inches at 300 dpi.](./figures/ref_marker_heatmap_3000x6000.png){width=50%}

The tutorial's `ggsave("Marker_Heatmap.png", width = 10, height = 10)` can only produce 3,000 × 3,000. The reference is twice as tall (SC-19), so it was **not** produced by the code in the tutorial — presumably a taller canvas was used to give the nine blocks room. This is why its comparison score with our 3,000 × 3,000 output is 0.7250 while everything else in the module sits above 0.96: the difference is aspect ratio, not data.
""",
    ),
    (
        "SC-04: the `-log10(p_val_adj)` axis",
        r"""
![Figure 23 — `Naive_Memory_CD4_DEG.png`: the volcano plot the tutorial asks for and never gets.](./figures/21_volcano_naive_cd4.png){width=100%}

X is `avg_log2FC` (naive CD4 T versus everything else), Y is `-log10(p_val_adj)`, grey means "not significant by the code's own criteria", red means `p_val_adj < 0.05 & |log2FC| > 0.5`. Of the 724 genes tested, 353 pass the adjusted-p criterion and 199 pass both — the vertical gap between the red envelope and the grey mass is the effect-size filter doing its work. The top-left corner is empty, which is the biological statement: naive CD4 T cells are defined by *having* IL7R, LEF1, CCR7 and *lacking* monocyte/B/NK genes, so their positive markers are the ones at the right. Two code details worth noticing: this cell only exists because `library(ggrepel)` was added (SC-16), and the Y axis is finite only because `p_val_adj` was floored at `.Machine$double.xmin` (SC-04) — with 353 tests across ~13k genes it is one rounding away from `Inf`.
""",
    ),
    (
        "This is what \"harmony worked\" looks like as a number instead of a picture",
        r"""
![Figure 24 — `UMAP_without_harmony.png`: the 25,894 merged ctrl+stim cells, coloured by `orig.ident`.](./figures/26_umap_no_harmony.png){width=100%}

Two colours, two lobes: the control and stimulated samples each occupy their own region and touch only at a thin seam. This is textbook batch confounding — the stimulation is a biological effect *and* a technical shift (different 10x run, different capture efficiency), and the uncorrected embedding lets the technical part dominate. It is also the figure the tutorial draws with **no QC at all** (SC-06): no mitochondrial cut, no feature cut, so the two lobes each carry their own burden of low-quality cells.

![Figure 25 — `UMAP_with_harmony.png`: the same cells after `RunHarmony(group.by.vars = "orig.ident")`.](./figures/27_umap_harmony.png){width=100%}

The colours interleave across every island, which is the visual claim "batch removed". Treat the picture as a hypothesis and the numbers in the table above as the test: the cross-sample neighbour fraction rises from 0.0287 (PCA space) to 0.4051 (harmony space), and after clustering the corrected pass finds 12 clusters of which **none** has more than 90% of its cells from one sample (mean |ctrl share − 0.5| = 0.030). One caveat that this figure cannot show: harmony uses all 50 PCs here while the UMAP uses the first 10 (SC-20), and two runs of the same code produce structurally different images (0.5631 between our own runs, SC-18) — so no conclusion about mixing should ever rest on pixel similarity.
""",
    ),
    (
        "the corrected script refuses to guess, keeping numeric labels and writing `cluster_label_audit_py.csv` for a human decision.",
        r"""
![Figure 26 — `rank_genes_groups` (Wilcoxon) for the eight leiden groups, 25 genes each.](./figures/22_scanpy_rank.png){width=100%}

The Python counterpart of `FindAllMarkers`, and the table the tutorial never writes to disk (SC-13). Each column is one group with its top genes ranked by score; read the columns top-down and ask whether the list is coherent — RPS12/RPS27/RPL32 at the top of group 0 means that group is defined largely by ribosomal content, which is a technical axis, not a cell type. This is the figure that makes the positional-naming bug (SC-12) visible: group 1's top genes (LDHB, LTB, CD3D, IL32, TPT1) are CD4 T-cell genes, yet the tutorial names group 1 `CD8 T`.

![Figure 27 — `rank_genes_groups_violin` for group 0: expression of its top ten genes, cluster by cluster.](./figures/23_scanpy_rank_violin.png){width=100%}

The violin view answers the question a ranked list cannot: is the "marker" specific, or merely higher? If a gene's distribution is shifted in every group with a long tail, it is a housekeeping-like difference; if it has a single tight peak in one group and a floor elsewhere, it is a real marker. Use this panel to sanity-check the first column of Figure 26 before trusting any label derived from it.

![Figure 28 — `dotplot`: 17 canonical markers across the eight leiden groups. Colour = mean expression, size = fraction of cells expressing.](./figures/24_scanpy_dotplot.png){width=100%}

This is the single most informative annotation figure in the Python half, and it directly contradicts two of the tutorial's eight names. Column `CD8 T` shows CD8A/CD8B at a fraction of cells and low mean, while IL7R is bright — a CD4 profile; column `NK` shows LYZ/CST3 and CD3E rather than GNLY/NKG7 dominance, i.e. a mixed CD8/NK group. The visual signature of a *correct* label is a dark small dot (high mean, few cells) or a large pale dot (moderate mean, many cells) confined to one column; two columns responding to the same marker set is the signal that the clustering merged two populations.

![Figure 29 — The final scanpy UMAP with the tutorial's eight labels drawn on the data.](./figures/25_scanpy_final_umap.png){width=100%}

Compare this with Figure 19 and the two methods' disagreements become concrete rather than abstract: here the T lobe is a single `CD4 T` group plus a `CD8 T` group that is actually CD4, NK is folded into CD8, and the platelet satellite is a 12-cell island. The picture is internally consistent and biologically wrong in two places — the exact failure mode that motivates the audit tables in `results/04_annotation/`.
""",
    ),
]

# ---------------------------------------------------------------- ZH blocks ----
ZH = [
    (
        "记录在 `results/01_qc/qc_rederivation_rawmatrix.csv`。",
        r"""
![图 1 — Seurat 对全部 2,700 个细胞的质控分布：每细胞基因数、UMI 数、线粒体占比，以及教材的三条阈值。](./figures/01_violin_mt.png){width=100%}

对着规则 `nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5` 来读。基因数分布是 2016 年浅测序数据的典型形状——中位 817，左侧长尾从未触到 200（所以那一侧**一个**细胞都没丢），右侧有一层薄薄的肩（5 个细胞，也是这套流程唯一会当成双细胞删掉的东西）。线粒体占比整体压在 5% 以下，右侧小肩就是被删的 57 个细胞。高 UMI 处**没有**第二个峰，这正是"三个数字就够了、不必上双细胞检测器"的依据。

![图 2 — scanpy 用 `calculate_qc_metrics` 算出的同样三个指标（在 `filter_cells(min_genes=200)` 与 `filter_genes(min_cells=3)` 之后）。](./figures/02_scanpy_violin_qc.png){width=100%}

这张面板不是上面那张的复制品，而差异恰好是要点：它的分布是在**过滤后的 13,714 基因**矩阵上算的，Seurat 那张用的是全部 32,738 个特征。线粒体占比是个比值，去掉 19,024 个全零基因就改变了边界细胞的分母——这正是"同样写 `< 5%`、Seurat 留 2,638 而 scanpy 留 2,635"的那一个细胞的来源（SC-10）。

![图 3 — 每个细胞的 `total_counts` 对 `n_genes_by_counts`。](./figures/03_scanpy_scatter_genes.png){width=100%}

这张图就是教材没点名的双细胞检验。主云几乎完全线性——每多一个 UMI 就按固定速率多拿基因——双细胞会落在线上方（同样的 UMI 数却拿到两个细胞的基因量）。这里主云紧致、上方没有分离出来的分支，所以"基因数上限 2,500"实际上是在代理双细胞过滤。

![图 4 — `total_counts` 对 `pct_counts_mt`。](./figures/04_scanpy_scatter_mt.png){width=100%}

读法相反：濒死细胞丢掉胞质 mRNA，线粒体转录本仍留在尚完整的细胞器里，于是往左上漂。5% 那条竖线是穿过稀疏上云画的，而不是穿过某个空隙——这是对任何线粒体阈值诚实的描述：它是一条约定，不是物理边界。
""",
    ),
    (
        "SC-08 在此：`p <- p + plot2` 画出的是**两联图**",
        r"""
![图 5 — 教材原文的 `VariableFeaturePlot.png`：两联图，因为 `p <- p + plot2` 把旧图和有标注的图并排打印了。](./figures/05_variable_features_aswritten.png){width=100%}

左右两半都是同一批 32,738 个基因、2,000 个 HVG 标红；左半没有基因名，右半标了二十个。这是全模块最干净的"代码能跑、产出看起来能用、但做的不是作者本意"的例子——而且下一张图会证明教师自己的参考图也带着同一个痕迹。

![图 6 — 课程随附的参考 `VariableFeaturePlot.png`：同样的两联版式，也就是说这个痕迹被当作"预期产出"发了出去。](./figures/ref_variable_featureplot_2panel.png){width=100%}

把它和图 5 并排看，才是判断"复现"的正确方式：我们的原样跑与之结构相关 0.9751（尺寸相同、点云相同、标注相同），因为两者是同一个两联 bug 的产物。反过来，一个学生"修好"这张图再与参考图对比，会得出"我的分析与老师不同"的错误结论。

![图 7 — 修正后的单联版（`VariableFeaturePlot_corrected.png`）。](./figures/06_variable_features_corrected.png){width=100%}

修正版用 `LabelPoints(repel = TRUE)` 画一张散点并把二十个标注挂上去。它与参考图的相关度按构造掉到 0.0713：两张图是同一份数据的不同画面。标注本身是体检项——PPBP、LYZ、S100A9、IGLL5、FTL、GNLY、PF4、FTH1、GNG11、S100A8、FCER1A、HLA-DRA、CD74、CLU、NKG7、CST3、GZMB、CCL4、C1QA、HLA-DPB1 全是白细胞基因，说明 `vst` 在 PBMC 上选对了。

![图 8 — 另一侧的 scanpy 选 HVG：标准化离散度对平均表达，选出 1,833 个基因。](./figures/07_scanpy_hvg.png){width=100%}

两种方法不可互换，这张图就是原因所在：`vst` 对"方差 ~ log 均值"的 loess 拟合做标准化，能提拔均值不高但极端的基因（PF4、PPBP）；scanpy 的 `min_mean 0.0125 / max_mean 3 / min_disp 0.5` 卡的是分箱后标准化数据的离散度，在 13,714 个基因的宇宙里返回 1,833 个。两个数都不错——它们用不同的方差模型回答同一个问题："哪些基因的变化超出其表达水平所预测的程度"。
""",
    ),
    (
        "但肘图画的是 20 维，而正文从未解释为什么取 10。",
        r"""
![图 9 — `PCA_VizDimLoadings`：PC1 与 PC2 上正负载荷最大的基因。](./figures/08_pca_loadings.png){width=100%}

载荷是 PCA 可解释的那一半。PC1 由单核-淋巴的对比主导（一侧 LYZ、S100A8、S100A9、CD14，另一侧 IL7R、CD3D、LTB），PC2 把 FCGR3A+ 单核与 NK 分开。如果某个主成分拿不出一份连贯的基因清单，那就是"该停止使用它"的信号——这里前几个成分在生物学上可读，这是选 `dims = 1:10` 的一个虽弱但真实的理由。

![图 10 — `PCA.png`：20 个主成分的标准差肘图。](./figures/09_pca_elbow.png){width=100%}

肘部弯在 PC 7-10 附近然后变平，这就是 `dims = 1:10`（用于 `FindNeighbors`、`FindClusters`、`RunUMAP`）的非正式依据。注意这张图**没有**做什么：它只画 20 维，而 `RunPCA` 算了 50 维；它也不给判据——取 10 是视觉惯例，SC-07 记录的就是教材从未把"画 20"和"用 10"对上。

![图 11 — `PCA_dimplot`：所有细胞落在 PC1/PC2 平面上，按 Seurat 群着色。](./figures/10_pca_dimplot.png){width=100%}

在 UMAP 出现之前，群就已经在两个主成分上分开了：单核那一臂、B/NK 那一臂、T 细胞那一区各占一个象限。PCA 做不到的是把 naive 与 memory CD4 分开——它们在中间重叠——而这恰恰需要十个主成分加一张图。

![图 12 — `PCA_heatmap_ALL`：前 6 个主成分，每个成分取该成分两端均匀分布的 100 个细胞。](./figures/11_pca_heatmap.png){width=100%}

每行一个基因、每列一个细胞；块状结构就是"好成分"的样子——一组共调控基因在被抽样的细胞里一起亮起来。这张图同时是我们最强的复现证据：与教师参考 PNG 的结构相关为 **1.0000**，即同一份 PCA、同样的平衡抽样、同一条绘图路径，文件几乎逐比特重现。

![图 13 — scanpy 的 PCA（PC1/PC2，按 CST3 上色），在 1,833 个 HVG 上、`sc.pp.scale(max_value = 10)` 之后计算。](./figures/12_scanpy_pca.png){width=100%}

与图 11 对照很有教益：同样的生物学，但轴的朝向与尺度不同——scanpy 在 1,833 个基因上把每个基因标准化到单位方差（截断 10），Seurat 先 scale 32,738 个再投影 2,000 个 HVG。哪个朝向都不是"标准"，主成分本身只定义到符号和旋转——所以两侧 PCA 长得不一样是预期结果，不是 bug。

![图 14 — `pca_variance_ratio`：教材自己的代码从未画出的那张图。](./figures/13_scanpy_pca_ratio.png){width=100%}

它在 28 张参考图里缺席的原因很朴素：画它的那个单元带着双逗号（SC-09），而 Python 以整单元格编译，于是 `SyntaxError` 把同一单元格里写着的 PCA、近邻图和 UMAP 一起抹掉了。我们这张来自修正版；累计曲线要到 PC 40 附近才到约 95%，这正是"近邻图用 `n_pcs = 40`"在这份数据上说得通的理由。
""",
    ),
    (
        "血小板群只占 0.5%，在 `resolution = 0.5` 下能保留，0.3 就没了。",
        r"""
![图 15 — 2,638 个细胞的 Seurat 9 群 UMAP（`UMAP.png`）。](./figures/14_umap_seurat.png){width=100%}

拓扑就是 PBMC 3k 的教科书图景：一大片 T 细胞叶（0、1、4），一条被大片空白隔开的单核臂（2、5），一座紧凑的 B 岛（3），加两颗小卫星——NK（6）和最外角上的 DC（7）与血小板（8）。第一件要检查的不是颜色，而是两个小群：2,638 个细胞里的 28 个与 14 个（1.1% 与 0.5%），它们在图上与一切清晰分离，marker 面板（FCER1A/CLEC10A 对 PPBP/PF4）毫不含糊。视觉上的告诫是：UMAP 的群间距离不是定量的，所以"cluster 2 离 cluster 0 很远"绝不能读成"单核与 CD4 T 的差异比 B 细胞更大"。

![图 16 — 同一份数据的 scanpy 版本：分辨率 0.9 的 leiden 加十二个 marker 面板（`umapscanpy_markers.png`）。](./figures/15_scanpy_markers_grid.png){width=100%}

一张顶九张：第一格是 leiden 解（教材阈值下 8 群），其余格是同一嵌入上叠加的 marker 基因。把它当作两个独立工具的一致性检查来读——B 细胞只在某一群里点亮 MS4A1/CD79A，单核点亮 CD14/LYZ/FCGR3A，血小板点亮 PPBP——同时注意它与 Seurat 不一致的地方：NK 的 marker（GNLY、NKG7、KLRF1）在这里跨了两个群，因为 40 维 PC 的近邻图把 CD8 T 与 NK 并成一个邻域，这个分辨率再按另一条边界把它切开。
""",
    ),
    (
        "对得**靠运气**：",
        r"""
![图 17 — `marker_Features`：教材注释表里点名的那九个谱系 marker。](./figures/16_featureplot_markers.png){width=100%}

人才是真正用这张图做注释的，所以值得逐格读。一个 marker 一个 marker 地看，注意它点亮了哪一群、**以及哪些群仍然是暗的**：CD14 与 LYZ 只在 cluster 2，MS4A1 只在 3，GNLY/NKG7 在 6，CD3E/IL7R 在 0 和 1，FCER1A/CST3 在 7，PPBP 在 8。由此可以得出两件事：第一，九个群里有七个光凭这张面板就已经没有歧义；第二，CD3E/IL7R **无法**把 0 与 1 分开——两者都是 T 细胞——这正是教材需要第二张面板（图 18）和差异分析（图 23）来区分 naive 与 memory 的原因。

![图 18 — `marker_Features_2`：同一嵌入上的另外十二个 marker。](./figures/17_featureplot_markers2.png){width=100%}

补充的 marker 解决了第一张面板解决不了的问题：FCGR3A+ 单核（cluster 5：FCGR3A、MS4A7、LST1）靠这张面板才与经典单核（cluster 2）分开；CCR7 与 LEF1 把 cluster 1 定为 naive CD4、S100A4 把 cluster 0 定为 memory；KLRB1/KLRD1 坐实 cluster 6 的 NK 身份。实用教训是：每个谱系至少要看**一对** marker，因为单个 marker 只是代理，而且不少是共享的（LYZ 在 naive CD4 T 里也表达 0.96，在单核里是 5.09）。

![图 19 — Seurat 对象的最终注释 UMAP（`cell annotations.png`）。](./figures/20_umap_annotated_seurat.png){width=100%}

九个标签直接画在数据上，与 §9.1 的审计表一一对应。把它当作下游所有论述的对照：单核臂是 `CD14+ Mono`（cluster 2，471 个细胞）与 `FCGR3A+ Mono`（cluster 5，175），上部 T 叶是 `Naive CD4 T`（0，657）加 `Memory CD4 T`（1，503），最远的 cluster 4 是 `CD8 T`（307）。后文任何涉及某个腔室的说法都可以用这里的颜色核对——这正是为什么 §13.2 里"30 个 HVG 扰动后 cluster 1/2 互换"这种错标在这张图里看起来完全正常。
""",
    ),
    (
        "由此得到的热图与参考图相关 0.9865。",
        r"""
![图 20 — 原样的 `Marker_Heatmap.png`：每群按 `avg_log2FC` 取前 10 个基因，取自带**未过滤**的 marker 表。](./figures/18_heatmap_aswritten.png){width=100%}

块状结构看起来令人放心——每群都有一组亮的基因、在别处都暗——但这张图正是 SC-03 的视觉形态。因为过滤算完就丢（`dplyr::filter()` 写在一条没人赋值的管道里），`slice_max(avg_log2FC)` 只按效应量排序，完全可以把 `p_val_adj = 0.2` 的基因与 `p_val_adj = 1e-40` 的基因并排放。读者从颜色上看不出这一点：图里画的是标准化计数，一个不显著的基因只要在少数细胞里高，看起来照样很有说服力。

![图 21 — 修正版热图：从通过 `p_val_adj < 0.05 & avg_log2FC > 1` 的 3,483 个 marker 里选 90 个。](./figures/19_heatmap_corrected.png){width=100%}

同样的九个块、同样的颜色、不同的基因清单。有意思的对比是哪些块**保持**清晰：B 细胞（CD79A/MS4A1/TCL1A）、单核（LYZ/S100A8/S100A9）与血小板（PPBP/PF4）在两版里几乎一样，因为它们的效应是压倒性的；变化最大的是 CD4 那两块，因为"naive vs memory"是细微对比（log2FC 2-4 而非 20），这正是显著性过滤的价值所在。对参考图相关 0.7586——可见差距来自被剔掉的基因，而不是聚类变了。

![图 22 — 课程随附的参考 `Marker_Heatmap.png`：3,000 × 6,000 px，即 300 dpi 下 10 × 20 英寸。](./figures/ref_marker_heatmap_3000x6000.png){width=50%}

教材的 `ggsave("Marker_Heatmap.png", width = 10, height = 10)` 只能产出 3,000 × 3,000。参考图高了一倍（SC-19），所以它**不是**教材代码的产物——想必是用了更高的画布给九个块腾地方。这就是为什么它与我们 3,000 × 3,000 输出的比对分数是 0.7250，而本模块其余图都在 0.96 以上：差别在长宽比，不在数据。
""",
    ),
    (
        "SC-04：`-log10(p_val_adj)` 在本数据上是有界值",
        r"""
![图 23 — `Naive_Memory_CD4_DEG.png`：教材要求画、却从未得到的那张火山图。](./figures/21_volcano_naive_cd4.png){width=100%}

横轴是 `avg_log2FC`（naive CD4 T 对其余全部细胞），纵轴是 `-log10(p_val_adj)`，灰色表示"按代码自己的标准不显著"，红色是 `p_val_adj < 0.05 & |log2FC| > 0.5`。724 个受检基因里 353 个通过校正 p 值、199 个同时通过两条——红色包络与灰色团之间那段竖向空隙，就是效应量过滤在起作用。左上角是空的，这是生物学陈述：naive CD4 T 由"**有** IL7R、LEF1、CCR7"和"**没有**单核/B/NK 基因"共同定义，所以正向 marker 都落在右侧。两个代码细节值得注意：这个单元之所以存在，是因为补了 `library(ggrepel)`（SC-16）；纵轴之所以有界，是因为把 `p_val_adj` 用 `.Machine$double.xmin` 托了底（SC-04）——在约 13k 个基因上做 353 次检验，它离 `Inf` 只差一次舍入。
""",
    ),
    (
        "这就是\"harmony 起作用了\"的数字版",
        r"""
![图 24 — `UMAP_without_harmony.png`：合并后的 25,894 个 ctrl+stim 细胞，按 `orig.ident` 上色。](./figures/26_umap_no_harmony.png){width=100%}

两种颜色、两个瓣：control 与 stimulated 各占一块，只在一条细缝处接触。这是教科书式的批次混淆——刺激既是生物学效应，也是技术漂移（两次 10x 运行、捕获效率不同），未校正的嵌入让技术那部分占了主导。这张图同时也是教材**完全不质控**画出来的（SC-06）：没有线粒体阈值、没有基因数阈值，于是两个瓣各自背着自己的低质量细胞。

![图 25 — `UMAP_with_harmony.png`：同样这批细胞 `RunHarmony(group.by.vars = "orig.ident")` 之后。](./figures/27_umap_harmony.png){width=100%}

颜色在每一座岛上都交错开，这就是"批次被去掉了"的视觉主张。把图当作假设、把上面那张表当作检验：跨样本近邻比例从 PCA 空间的 0.0287 升到 harmony 空间的 0.4051；修正版聚类后得到 12 个群，**没有一个**群的单样本占比超过 90%（ctrl 占比偏离 0.5 的平均绝对值 0.030）。这张图展示不了的两点告诫：harmony 用了全部 50 个 PC 而 UMAP 只用前 10（SC-20）；同一份代码跑两次会给出结构不同的图（我们自己两次之间 0.5631，SC-18）——所以任何关于混合程度的结论都不该建立在像素相似度上。
""",
    ),
    (
        "修正版拒绝猜测，保留数字标签并写出 `cluster_label_audit_py.csv` 交给人判断。",
        r"""
![图 26 — `rank_genes_groups`（Wilcoxon）对 8 个 leiden 群、每群 25 个基因。](./figures/22_scanpy_rank.png){width=100%}

`FindAllMarkers` 的 Python 对应物，也是教材从不落盘的那张表（SC-13）。每一列是一个群、按得分排的头部基因；自上而下读每一列并问"这份清单连贯吗"——group 0 顶部是 RPS12/RPS27/RPL32，说明这一群很大程度由核糖体含量定义，那是技术轴而非细胞类型。这张图让位置式命名 bug（SC-12）变得可见：group 1 的头部基因（LDHB、LTB、CD3D、IL32、TPT1）是 CD4 T 的基因，而教材把 group 1 命名为 `CD8 T`。

![图 27 — `rank_genes_groups_violin`（group 0）：它的前十个基因在各群中的表达分布。](./figures/23_scanpy_rank_violin.png){width=100%}

小提琴图回答了排名清单回答不了的问题：这个"marker"是特异的，还是只是更高一点？如果某基因在所有群里都整体右移且长尾，那是类管家式的差异；如果它在某一群里是单个紧峰、在别处贴着地板，那才是真 marker。用这一格去核对图 26 的第一列，再决定是否相信据此推出的标签。

![图 28 — `dotplot`：17 个经典 marker × 8 个 leiden 群。颜色 = 平均表达，大小 = 表达细胞比例。](./figures/24_scanpy_dotplot.png){width=100%}

这是 Python 侧信息量最大的注释图，而且它直接反驳了教材八个名字里的两个。`CD8 T` 这一列里 CD8A/CD8B 只有少量细胞表达且均值低，而 IL7R 很亮——这是 CD4 的谱型；`NK` 这一列里主导的是 LYZ/CST3 与 CD3E，而不是 GNLY/NKG7，即一个 CD8/NK 混合群。**正确**标签的视觉签名是"小而深的点（均值高、细胞少）"或"大而浅的点（均值中等、细胞多）"且只出现在一列；如果两列对同一批 marker 都有响应，那就是聚类把两个群体合并了的信号。

![图 29 — scanpy 的最终 UMAP，把教材那八个标签画在数据上。](./figures/25_scanpy_final_umap.png){width=100%}

把它与图 19 并排，两个方法的分歧就从抽象变具体了：这里 T 叶是一个 `CD4 T` 群加一个其实属于 CD4 的 `CD8 T` 群，NK 被折进 CD8，血小板是 12 个细胞的孤岛。这张图内部自洽、而在两个位置生物学上错误——这正是 `results/04_annotation/` 里那些审计表存在的理由。
""",
    ),
]


def apply(md, blocks, tag):
    p = HERE / md
    text = p.read_text()
    for anchor, block in blocks:
        text = insert_after_para(text, anchor, block, tag)
    p.write_text(text)
    n = text.count("./figures/")
    print(f"{md}: inserted {len(blocks)} blocks, {n} figure links")


if __name__ == "__main__":
    apply("scRNA_seq_analyze_corrected.md", EN, "EN")
    apply("scRNA_seq_analyze_corrected_zh.md", ZH, "ZH")
