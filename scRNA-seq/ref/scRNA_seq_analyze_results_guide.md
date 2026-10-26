# scRNA-seq 分析结果讲解

本文件导读 `scRNA-seq/results/` 下的全部产物，并给出可以亲手复现的练习。配套的矫正增强版文档在 `ref/scRNA_seq_analyze_corrected.md`（英文）与 `ref/scRNA_seq_analyze_corrected_zh.md`（中文）。

## 1. 怎么读这个目录

```
results/
├── 01_qc/                  Part 2：质控图 + 两份账目表（运行值与原始矩阵重算值）
├── 02_hvg_pca/             Part 3-4：高变基因图、PCA 四张图；两侧、两版并存
├── 03_clustering/          Part 5：UMAP 图 + 三种聚类方案的群大小表
├── 04_annotation/          Part 6：注释图、marker 表、审计表、跨工具列联表
├── 05_deg/                 Part 6：火山图 + 两张 DEG 表 + scanpy 排名/小提琴/点图
├── 06_integration/         Part 7：去批次前后 UMAP + 混合度指标 + 分样本群占比
├── 07_objects/             Part 8：pbmc3k.h5ad + 逐细胞元数据表
├── 08_figure_comparison/   与课程随附 28 张参考图的机器可读/书面比对
├── 09_run_logs/            六次运行的逐单元日志与状态表
└── scripts/                重新生成以上全部内容的 13 个脚本
```

约定：

- 文件名里带 `_corrected` 的是修正版；不带的是**原样跑**的产出。
- Seurat 侧的图直接在 `scripts/run_cells.R` 的工作目录里生成；scanpy 侧的图在 `figures/` 子目录里，文件名保留了 scanpy 的 `save=` 前缀（`highest_expr_genesscanpy_...png`），因为这就是教材代码的行为。
- 每个 `09_run_logs/*_cell_log.txt` 里，每个单元都有 `--status--`、`--warnings`、`--state--` 三段；`--state--` 记录该单元结束后对象是多少细胞、多少群、什么标签。

## 2. 逐类导读

### 2.1 `01_qc/`

`violin_mt.png` 是 Seurat 的三个小提琴图（基因数、UMI、线粒体占比），`violinscanpy_violin_plot_QC.png` 与两张 `scatterscanpy_scatter_*.png` 是 scanpy 侧对应的质控视图。两张散点的用途是**在定阈值之前**看形状：`total_counts` vs `n_genes_by_counts` 应该是线性关系，偏离线性且在右上的点是双细胞；`total_counts` vs `pct_counts_mt` 右上角是濒死细胞。

表比图重要：

- `qc_summary.csv`：四次运行（Seurat 原样、Seurat 对齐基因全集、scanpy 教材阈值、scanpy 对齐阈值）各自的规则、输入/保留细胞数、基因数，最后一列给出该数字出自哪份日志。
- `qc_rederivation_rawmatrix.csv`：在原始矩阵上按同样数字独立重算，用来解释"为什么同样写 `< 5%`，Seurat 留 2,638 个细胞而 scanpy 留 2,635 个"。结论：线粒体占比在 scanpy 里是按**过滤后**的 13,714 基因矩阵算的，边界细胞的占比会落在 5% 两侧。
- `qc_metric_distribution.csv`：三个指标在全部细胞与保留细胞上的中位/最小/最大，用来判断阈值卡得松还是紧。

### 2.2 `02_hvg_pca/`

Seurat 侧：`VariableFeaturePlot.png`（原样，**两联**——左边是 `p`，右边是 `p + 标注`，教材 `p <- p + plot2` 的真实效果）、`VariableFeaturePlot_corrected.png`（单联，修正后）、`PCA_VizDimLoadings.png`（PC1/PC2 的载荷，可直接读出哪些基因在驱动主成分）、`PCA.png`（肘图，20 维）、`PCA_dimplot.png`（前两个 PC 上的细胞散点，9 个群已经能看出轮廓）、`PCA_heatmap_ALL.png`（6 个 PC × 每个 100 个平衡细胞的热图，Seurat 社区里常说的"看 PC 上有没有结构"就看它）。

scanpy 侧：`filter_genes_dispersionscanpy_HVG.png`（均值-离散度散点，HVG 被标红）、`pcascanpy_PCA.png`（PC1/PC2，按 CST3 上色）、`pca_variance_ratioscanpy_PCA_ratio.png`（碎石图，**参考图里没有这一张**——教材那个单元格因双逗号整段没跑）。

### 2.3 `03_clustering/`

`UMAP.png` 是 Seurat 9 群的 UMAP；`umapscanpy_markers.png` 是 scanpy 的 13 联图（leiden + 12 个 marker）。

`cluster_sizes.csv` 把三种聚类方案并排放：

- Seurat：9 群，657 / 503 / 471 / 339 / 307 / 175 / 144 / 28 / 14；
- scanpy（教材阈值 2,000）：8 群，534 / 643 / 340 / 497 / 432 / 142 / 35 / 12；
- scanpy（对齐阈值 2,500）：9 群，559 / 615 / 340 / 209 / 155 / 281 / 430 / 36 / 13。

读法：先看最小的两群（DC 28-36 个、血小板 12-14 个）是否在三种方案里都稳定存在——在，说明它们是真实细胞类型而不是噪声；再看两个大 T 腔室是否被拆分——Seurat 拆了 naive/memory CD4，scanpy 拆的是"核糖体高"的一个 CD4 群。

### 2.4 `04_annotation/`

图：`marker_Features.png`（9 个谱系 marker）、`marker_Features_2.png`（12 个细分 marker）、`Marker_Heatmap.png`（原样：top10 只看 log2FC）与 `Marker_Heatmap_corrected.png`（先卡 `p_val_adj < 0.05 & log2FC > 1` 再取 top10）、`cell annotations.png`（Seurat 9 群带标签 UMAP）、`umapscanpy_final_umap.png`（scanpy 8 群带标签 UMAP）。

表：这一节真正的价值在表里。

- `marker_panel_mean_expression_seurat.csv` / `..._scanpy.csv`：一个 40 个 marker × 各群均值的面板，一眼能看出每个群到底是什么。例如 Seurat 的 `Platelet` 群 PPBP 5.85、PF4 4.93，而 RPS12 只有 1.34（血小板核糖体含量低）。
- `cluster_label_audit.csv`：教材的两份硬编码映射（cell 39 的表 vs cell 43 的向量）并排，加上观测到的程序与该群的 top5 marker。本次运行里 cluster 1 是 CD4、cluster 2 是单核，所以**向量对、表错**。
- `cluster_label_audit_py.csv`：scanpy 侧的同类审计，8 个群各自的 top5 marker 与细胞数。
- `cluster_correspondence_counts.csv` 及两份 share 表：Seurat 标签 × scanpy 群的列联表。`CD8 T` 这一列 78% 是 Seurat 的 naive CD4 T，是教材位置式命名出错最直接的证据。
- `markers_all_clusters.csv`（11,986 行）、`markers_top10_per_cluster.csv`（90 行）、`cluster0_markers.csv` … `cluster7_markers.csv`（scanpy，每群前 200 个）：教材算了但从未落盘的 marker 统计。
- `appendix_mincells3_*.csv`：基因全集改成 13,714 后的同一批表，用来看"改动 30 个 HVG 会让什么跟着动"。

### 2.5 `05_deg/`

`Naive_Memory_CD4_DEG.png` 是教材本想画、但原样跑崩溃的那张火山图（`library(ggrepel)` 缺失）；横轴 `avg_log2FC`、纵轴 `-log10(p_adj)`，灰色是有差异但不显著的基因，红色是显著项。

`DEG_NaiveCD4_vs_rest.csv`（724 个基因）与 `DEG_MemoryCD4_vs_NaiveCD4.csv`（650 个基因）是两张明细表，含 `p_val`、`p_val_adj`、`avg_log2FC`、`pct.1`、`pct.2`。scanpy 侧的 `rank_genes_groups_leidenscanpy_rank_genes_groups.png`（8 群各列前 25）、`..._violin.png`（群 0 的前 10 个基因表达分布）、`dotplot_scanpy_rank_genes_groups_dotplot.png`（17 个经典 marker × 8 群）是同一件事的另一种算法。

### 2.6 `06_integration/`

`UMAP_without_harmony.png` 与 `UMAP_with_harmony.png` 是教材唯一教的两张图，按 `orig.ident`（ctrl / stim）上色。原样跑只到这里；**看图得到的结论是"去批次成功了"，但图本身不可跨运行复现**（同一输入两次运行结构相关只有 0.563），所以本节配了数字：

- `integration_metrics.csv`：ctrl 与 stim 细胞数（QC 后 12,932 / 12,779）、合并后 25,711、跨样本近邻比例（PCA 空间 0.0287 → harmony 空间 0.4051）、harmony 聚类数 12、单样本占比 >90% 的群数 0。
- `harmony_cluster_by_sample.csv`：12 个群里 ctrl/stim 各多少细胞，最大群 6,065（3,096 / 2,969），最小群 92（38 / 54）。

读法：跨样本近邻比例从 3% 升到 41%，说明两批细胞被拉到一起了；同时每个群内部两样本接近 1:1（偏离 0.5 的平均绝对值 0.030），说明没有把某一种细胞"整合没了"。

### 2.7 `07_objects/`

`pbmc3k.h5ad` 是 Part 8 的交付对象：2,635 个细胞、1,833 个高变基因（`X` 是 scale 后的矩阵，`raw` 保留全部 13,714 个基因的 log 标准化值），`obs` 里有 `leiden`（已改名为 8 个细胞类型）与全部 scanpy 质控指标。`pbmc3k_obs.csv` 是同一份 `obs` 的纯文本版，方便不装 anndata 的人读。

`pbmc3k_final.rds` 没有入库：Seurat 为全部 32,738 个基因保留稠密 `scale.data`，文件 331 MiB；Part 6 里那行 `saveRDS()` 一行即可重建。

### 2.8 `08_figure_comparison/`

四张表 + 一份书面报告：

- `compare_aswritten_vs_reference.tsv`：原样跑 vs 28 张参考图；
- `compare_corrected_vs_reference.tsv`：修正版 vs 参考图；
- `compare_appendix_mincells3_vs_reference.tsv`：两个附录运行 vs 参考图；
- `compare_aswritten_vs_corrected.tsv`：我们自己的两版之间的逐字节比较。

每行给出像素尺寸、结构相关（1.0000 = 降采样后结构一致）、着色数量与各色占比。`figure_comparison.md` 是同一批数据的解读，包括"为什么热图的 0.7250 是长宽比造成的"和"为什么两张 harmony UMAP 的 0.56 等于噪声"。

### 2.9 `09_run_logs/` 与 `scripts/`

六次运行各有一份 `*_cell_log.txt`（逐单元的代码、stdout、警告、报错、对象状态）与 `*_cell_summary.tsv`（每单元一行）。`R_tutorial_*` 是原样跑，`R_corrected_mincells0_*` 是交付版修正跑，`R_corrected_mincells3_*` 是基因全集附录，`Py_asis_*` / `Py_fixed_2000_*` / `Py_fixed_2500_*` 同理。

`scripts/` 里 13 个文件分三类：教材单元文件（`*_cells_aswritten.*`、`*_cells_corrected.*`）、驱动器（`run_cells.R`、`run_cells.py`）、分析脚本（`table_qc_clusters.R`、`marker_panel_audit.*`、`annotation_audit.R`、`cluster_correspondence.*`、`compare_figs.py`）。所有脚本都用绝对路径或命令行参数，不依赖当前目录。

## 3. 手把手练习

**练习 1：证明原样跑确实少了 7 张图。**
```bash
grep -c "^  --status-- ERROR" results/09_run_logs/R_tutorial_cell_log.txt   # 1
grep -c "^  --status-- ERROR" results/09_run_logs/Py_asis_cell_log.txt      # 4
```
R 侧 13 张图 + scanpy 侧 5 张 = 18 张；`25 - 18 = 7`。再看 `Py_asis_cell_log.txt` 里 cell 67 的报错行，注意报错类型是 `SyntaxError`——整个单元格一行都没执行。

**练习 2：算一下少的那张火山图值不值得补。**
从 `05_deg/DEG_NaiveCD4_vs_rest.csv` 里数出 `p_val_adj < 0.05 & |avg_log2FC| > 0.5` 的基因数（199），再数只满足 `p_val_adj < 0.05` 的（353）。差值说明"显著"到底被 log2FC 门槛砍掉了多少。

**练习 3：验证"CD8 T"那个标签是错的。**
```bash
awk -F, 'NR==1 || $1=="CD8A" || $1=="CD8B" || $1=="IL7R"' \
    results/04_annotation/marker_panel_mean_expression_scanpy.csv
```
scanpy 侧名为 `CD8 T` 的群里 CD8A = 0.14、CD8B = 0.24，而 IL7R = 1.61——它是 CD4 T 群。再对照 `cluster_correspondence_counts.csv`：这一列 78.1% 是 Seurat 的 Naive CD4 T。

**练习 4：把两份硬编码映射的冲突找出来。**
```bash
column -s, -t results/04_annotation/cluster_label_audit.csv | cut -c1-120
```
看 `table_label` 与 `rename_label` 两列在第 1、2 行（cluster 1、2）的差别，再看 `observed_programme` 判定谁对。然后打开教材 cell 39 与 cell 43，确认这两份清单确实同时存在。

**练习 5：复现"30 个 HVG 改变 36% 的标签"。**
比较 `results/04_annotation/cluster_label_audit.csv` 与 `appendix_mincells3_cluster_label_audit.csv` 的 `top5_markers` 列：对照 run 里 cluster 1 是 CD4（CCR7/LEF1），附录 run 里 cluster 1 是单核（S100A8/FCN1）。教材的固定向量在后者会把 957 个细胞标反。

**练习 6：理解为什么 harmony 的图不能比像素。**
```bash
awk -F'\t' 'NR==1 || $1=="UMAP_without_harmony.png" || $1=="UMAP_with_harmony.png"' \
    results/08_figure_comparison/compare_aswritten_vs_corrected.tsv
```
我们自己的两次运行（同代码同输入）相关只有 0.5631 / 0.5676，与对参考图的 0.5664 / 0.5689 同量级。改用 `integration_metrics.csv` 的跨样本近邻比例做判据。

**练习 7：把 scanpy 的阈值改成 2,500，看群数。**
```bash
mkdir -p /tmp/t2500 && cd /tmp/t2500 && N_GENES_MAX=2500 MPLBACKEND=Agg \
  /lustre1/share/miniconda3/envs/scRNA_seq/bin/python \
  /path/to/results/scripts/run_cells.py \
  /path/to/results/scripts/scRNA_seq_py_cells_corrected.py /tmp/t2500 Py2500
```
结果：2,638 个细胞（与 Seurat 完全相同）、9 个 leiden 群，8 个硬编码名字不再够用——修正版会保留数字标签并写出审计表。

**练习 8：给自己的数据定阈值。**
读一组新的 10x 数据，先画 `nFeature` vs `nCount` 与 `nFeature` vs `mt%` 两张散点，再决定上下限；然后用 `qc_rederivation_rawmatrix.csv` 的列结构记下每个阈值丢了多少细胞，别只记留下的数量。

## 4. 一句话总结

这 25 张图里有 18 张是教材代码原样能产出的，另 7 张死在两个字符级的错误上（少一行 `library(ggrepel)`、多一个逗号）；而即使产出的那 18 张，也有两张注释的映射互相矛盾、8 个 scanpy 标签错 2 个——所以本目录里每个"看起来对"的图，都配了一张能查的表。
