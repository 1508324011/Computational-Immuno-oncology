# Figure comparison: our re-run versus the 28 reference PNGs

The course directory ships 28 reference figures (the teacher's 2025 output). 25 of them are
products of the tutorial's own code; three are lecture/demo material that no tutorial cell
produces (see `ref/scRNA_seq_analyze_corrected.md` §14).

Three machine-readable tables sit next to this file:

| file | what it compares |
|---|---|
| `compare_aswritten_vs_reference.tsv` | the tutorial run **as written** against the reference |
| `compare_corrected_vs_reference.tsv` | the **corrected** run against the reference |
| `compare_appendix_mincells3_vs_reference.tsv` | the **appendix** runs (gene universe aligned to scanpy, threshold aligned to Seurat) |
| `compare_aswritten_vs_corrected.tsv` | our own two runs against each other (byte identity + structural correlation) |

Method: pixel dimensions from the PNG header, then a 64 × 64 grayscale downsample,
mean-subtracted and L2-normalised, dotted with the reference (so 1.0000 = same structure
after downsampling; PNG re-encoding, anti-aliasing and font substitution all survive it).
For point-cloud panels the script also counts distinct saturated colours and their pixel
shares, which is what a reader actually sees in a UMAP.

---

## 1. As-written run: 18 of 25 figures, and how close they are

| figure | px (ours / ref) | structural correlation |
|---|---|---|
| `PCA_heatmap_ALL.png` | 3000×3000 / 3000×3000 | **1.0000** |
| `violin_mt.png` | 3000×1500 / 3000×1500 | 0.9982 |
| `marker_Features_2.png` | 6000×3000 / 6000×3000 | 0.9965 |
| `marker_Features.png` | 3000×3000 / 3000×3000 | 0.9957 |
| `UMAP_with_harmony.png` | 3000×1500 / 3000×1500 | 0.9953 |
| `UMAP_without_harmony.png` | 3000×1500 / 3000×1500 | 0.9941 |
| `UMAP.png` | 3000×1500 / 3000×1500 | 0.9914 |
| `cell annotations.png` | 3000×1500 / 3000×1500 | 0.9898 |
| `PCA_dimplot.png` | 3000×1500 / 3000×1500 | 0.9875 |
| `PCA_VizDimLoadings.png` | 3000×1500 / 3000×1500 | 0.9775 |
| `VariableFeaturePlot.png` | 3000×1500 / 3000×1500 | 0.9751 |
| `PCA.png` | 3000×1500 / 3000×1500 | 0.9603 |
| `Marker_Heatmap.png` | 3000×**3000** / 3000×**6000** | 0.7250 |
| `Naive_Memory_CD4_DEG.png` | — / 1800×1500 | not produced (cell 46 died) |
| scanpy `highest_expr_genes` | 474×486 / 474×486 | **1.0000** |
| scanpy `scatter` × 2 | 612×413, 630×413 | **1.0000** |
| scanpy `HVG` | 1061×432 / 1061×432 | 0.9998 |
| scanpy `violin_plot_QC` | 1508×490 / 1508×490 | 0.9992 |
| scanpy `PCA`, `markers`, `rank`, `rank_violin`, `dotplot`, `final_umap` | — | not produced (cell 67 SyntaxError) |

Interpretation:

* **The deterministic half of the pipeline is bit-reproducible.** `PCA_heatmap_ALL.png` is
  1.0000 because it is the same PCA on the same 100 balanced cells per component: exact
  reproduction is possible and was achieved. The two scanpy scatter plots and the
  highest-expression bar chart are likewise 1.0000.
* **Everything derived from a graph or an embedding sits at 0.96-0.998.** Louvain and UMAP
  add a small amount of layout noise between any two runs, including two of ours.
* **`Marker_Heatmap.png` at 0.7250 is an aspect-ratio artefact, not a data difference.**
  The reference is 3,000 × 6,000 px (10 × 20 in at 300 dpi) while the tutorial's
  `ggsave(..., width = 10, height = 10)` gives 3,000 × 3,000. Squeezing a 2:1 image into a
  1:1 comparison costs most of the correlation.
* **The two Part-7 UMAPs (0.9941, 0.9953) are better than our own run-to-run agreement**
  (§3), i.e. the batch-corrected embedding is reproducible here even though uwot is not
  bitwise-reproducible in general.

---

## 2. Corrected run: 25 of 25, plus the figure the tutorial never rendered

Changes with respect to the as-written column:

| figure | correlation | why it moved |
|---|---|---|
| `VariableFeaturePlot.png` | 0.9751 → **0.0713** | the correction removes `p + plot2`, which draws two panels. The reference ships with the two-panel artefact (SC-08), so the "correct" figure necessarily looks different. The single-panel version is kept as `VariableFeaturePlot_corrected.png`. |
| `Marker_Heatmap.png` | 0.7250 → 0.7586 | the heatmap now uses only genes with `p_val_adj < 0.05 & avg_log2FC > 1` (3,483 of 11,986) before taking the top 10 per cluster (SC-03). |
| `Naive_Memory_CD4_DEG.png` | — → 0.3255 | produced for the first time (`library(ggrepel)`); the label layout of `geom_text_repel` is randomised, so a low correlation is expected even for an identical volcano. |
| `pca_variance_ratioscanpy_PCA_ratio.png` | — → *no reference* | the tutorial's own code never rendered it because of the `SyntaxError`; the reference set consequently lacks the file. |
| scanpy `final_umap.png` | — → 0.4915 | depends on 8 leiden communities; the width difference (516 vs 533 px) reflects the legend labels being drawn on the data. |
| scanpy `markers.png` | — → 0.6824 | 13-panel feature grid on the leiden embedding; panel geometry identical, point cloud already differs at this zoom. |
| scanpy `rank_genes_groups.png` | — → 0.7669 | same 8 groups, same height (822 px); the width differs by 13 px because the gene labels differ at the 25th position. |

---

## 3. As-written versus corrected: our own two runs

| figure | byte-identical | structural correlation |
|---|---|---|
| `violin_mt.png` | **yes** | 1.0000 |
| `PCA_VizDimLoadings.png` | **yes** | 1.0000 |
| `PCA.png` | **yes** | 1.0000 |
| `PCA_dimplot.png` | **yes** | 1.0000 |
| `PCA_heatmap_ALL.png` | **yes** | 1.0000 |
| `UMAP.png` | **yes** | 1.0000 |
| `marker_Features.png` | **yes** | 1.0000 |
| `marker_Features_2.png` | **yes** | 1.0000 |
| `cell annotations.png` | **yes** | 1.0000 |
| `Marker_Heatmap.png` | no | 0.9865 |
| `UMAP_without_harmony.png` | no | 0.5631 |
| `UMAP_with_harmony.png` | no | 0.5676 |
| `VariableFeaturePlot.png` | no | 0.0825 |

Nine figures are byte-identical, which is the strongest available statement that the
corrected script does not perturb the analysis. The three that differ are the two
deliberate changes and the part-7 UMAPs.

**The Part-7 UMAPs are not reproducible, and this is measurable.** Our two runs use
identical input, identical code and identical arguments, yet agree only at 0.5631 and
0.5676 — statistically indistinguishable from their agreement with the teacher's reference
figures (0.5664 and 0.5689). `RunUMAP()` calls uwot, which is multithreaded and has no
thread-invariant seeding, so a 25,711-cell embedding is redrawn on every run. Any claim of
"the UMAP matches" for this dataset is therefore unfounded; what can be compared is the
cluster composition (12 clusters, none dominated by one sample, mean |ctrl share − 0.5| =
0.030; see `results/06_integration/integration_metrics.csv`).

---

## 4. Appendix runs: what happens if the two halves are aligned

`MIN_CELLS=3` on the R side (13,714 features, the scanpy gene universe):

| figure | correlation with reference |
|---|---|
| `PCA.png` | 0.9603 → 0.8974 |
| `PCA_dimplot.png` | 0.9875 → 0.3929 |
| `PCA_heatmap_ALL.png` | 1.0000 → 0.9371 |
| `UMAP.png` | 0.9914 → 0.1176 |
| `cell annotations.png` | 0.9898 → 0.1173 |

Removing 19,024 all-zero features changes 30 of the 2,000 HVGs (98.5% overlap) and that is
enough to renumber Louvain's communities: the clustering keeps 9 clusters but with sizes
684/481/476/344/291/162/155/32/13 instead of 657/503/471/339/307/175/144/28/14, and the
tutorial's fixed name vector then swaps the monocyte and memory-CD4 clusters (957 cells,
36% of the dataset). This is why the delivered corrected run keeps `MIN_CELLS = 0`: the
gene-universe alignment is a legitimate choice, but it is not free and it must not be made
silently (`ref/scRNA_seq_analyze_corrected.md` §13.2).

`N_GENES_MAX=2500` on the Python side (the Seurat threshold) gives 2,638 cells — exactly
the Seurat count — and 9 leiden communities instead of 8, which no longer fits the eight
hard-coded names. That run is kept as the appendix `py_align2500`; its figures are in the
third table.
