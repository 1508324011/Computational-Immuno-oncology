# 单细胞 RNA-seq（PBMC 3k）· 矫正增强版

**第五节课 · scRNA-seq 分析**
教材来源：*Single Cell RNA-seq Analyze Tutorial*，秦博文、曾泽贤，北京大学生命科学中心 / 定量生物学中心 / 前沿交叉学科研究院（本目录下 `scRNA_seq_analyze_totorial_2026.html`）。

---

## 0. 关于本版

本文件是第五节课教材的矫正增强版，沿用本仓库 WES、RNA-seq、RNA-seq-application、TBCR-analysis 四个模块的既定流程：

1. 教材先**按原样执行**——每个代码单元逐字照抄，按教材顺序，不做优化、不在运行中修错；
2. 所有遇到或预见的问题**只记录不修改**（清单见 §1，具体解释放在对应 Part）；
3. 随后**另跑一份修正版**，让每条修正的影响可被测量而不是被断言（§13）；
4. 文档、脚本、日志、图、表一并入库。

本模块与前四节课有两处不同。

**第一次同一份数据上有两套工具链。** Part 1-7 是 R/Seurat，Part 8 是 Python/scanpy，都从同一个 `matrix.mtx` 出发。教材从不交代两侧的关系，而它们确实不一致：R 侧保留 32,738 个基因、得到 9 个群，Python 侧保留 13,714 个基因、得到 8 个群。§12 把这一差异逐细胞量化。

**第一次教材在两门语言里都跑不到结尾。** 少一行 `library()` 让 R 的差异表达单元崩掉，多一个逗号让整个 Python 单元崩掉——合起来就是 25 张图里的 7 张，以及两个交付对象全部丢失（§3）。

本版相对原教材增加的内容：

- 完整的一次运行：逐单元耗时、警告、以及准确的报错原文（§3）；
- 二十条已验证的坑（`SC-01` … `SC-20`），每条给出修正代码，重要处给出实测影响；
- 每张重生成的图与课程随附 28 张参考 PNG 的逐图像素级比对（§3.4，`results/08_figure_comparison/`）；
- 两套工具链的对账：所谓"9 群 vs 8 群"到底是什么（§12）；
- 一份修正版，原样复现了全部数字，并补出教材只留在内存里的对象与表格（§13）。

全文统一四档严重度：**致命**（原样代码根本产不出结果）、**改变结果**（能跑，但标签/数字不可信）、**稳健性**（这里能跑，别处会崩）、**代码卫生**（风格、命名、结果不落盘）。

---

## 1. 问题清单速览

| 编号 | Part | 问题 | 严重度 | 修正 |
|:-------|:------|:---------------------------------------------|:-----------|:---------------------------------|
| SC-01 | 4 | `ScaleData(features = rownames(seurat_obj))` 对全部 32,738 个基因做 scale，而正文写的是"用选出的特征"（2,000 个 HVG）。这不是错误：`DoHeatmap()` 从 `scale.data` 取基因，只 scale HVG 会让任意 marker 热图直接报错。代价是：对象 0.73 GB、耗时 4.2 s，而 13,714 基因版是 0.35 GB / 2.2 s。 | 稳健性 | 保留，但写明理由；内存敏感时 scale `unique(c(VariableFeatures, 需要的 marker))`。 |
| SC-02 | 6 | cluster → 细胞类型的映射在教材里**写了两份且互相矛盾**：cell 39 的表说 cluster 1 = CD14+ Mono、cluster 2 = Memory CD4 T；cell 43 的 `new.cluster.ids` 说 cluster 1 = Memory CD4 T、cluster 2 = CD14+ Mono。本次运行里向量对、表错（cluster 1 是 CCR7/LEF1/NOSIP 即 CD4，cluster 2 是 S100A8/FCN1/CD14 即单核）。而 Louvain 编号并不稳定：§13.2 里仅 30 个 HVG 的扰动，就让同一个固定向量标错 957 个细胞（36%）。 | 改变结果 | 由 marker 反推标签、断言群数、输出审计表。 |
| SC-03 | 6 | `seurat_obj.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)` 是一条**值被丢弃的管道**；随后 `top10` 用 `slice_max(avg_log2FC)` 取自**未过滤**的表，完全不看 `p_val_adj`。过滤从未生效。 | 改变结果 | 先赋值过滤结果、同时卡 `p_val_adj < 0.05`，再 `slice_max()`。 |
| SC-04 | 6 | `naive_cd4_genes`、`memory_vs_naive_genes` 算完再也不使用；火山图纵轴 `-log10(p_val_adj)`，当 `p_val_adj` 下溢为 0 时是 `Inf`（本数据未触发——最小 1.24e-115——但测试精度再少一位就会触发）；两个 DEG 表都不落盘。 | 稳健性 | 用 `pmax(p_val_adj, .Machine$double.xmin)`；导出两张表；真正用上基因列表。 |
| SC-05 | 7 | `data_pbmc_multi <- load("pbmc_stim.RData")` 拿到的是被载入对象的**名字**，不是对象；`ctrl.sparse`/`stim.sparse` 被倒进全局环境，再用 `get()` 取回。放进函数或 knit 文档里会静默泄漏或直接失败。 | 稳健性 | `load(..., envir = e)` 载入私有环境。 |
| SC-06 | 7 | Part 7 直接把**未过滤的原始矩阵**送去整合——没有 `percent.mt`、没有 `nFeature` 阈值，与 Part 2 完全脱节；`RunHarmony()` 之后只画一张 UMAP，不做 `FindNeighbors()`/`FindClusters()`，去批次效果从未被量化。 | 改变结果 | 两侧都做 QC（13,019 → 12,932、12,875 → 12,779），并用量化指标替代肉眼判断（§10）。 |
| SC-07 | 4 | `PCHeatmap()` 是已废弃别名（Seurat 5.1.0 里仍在、且不报警告）；`ElbowPlot()` 默认画 20 个 PC，而下游全部用 `dims = 1:10`。 | 稳健性 | 换 `DimHeatmap(dims = 1:6)`；显式 `ndims = 20` 并说明为何取 10。 |
| SC-08 | 3 | `p <- p + plot2` 中 `plot2` 本身是基于 `p` 的完整 ggplot，结果是**两联图**（左 `p`、右 `plot2`），而不是"同一个图上叠了一层"。随附的参考图 `VariableFeaturePlot.png` 就带着这个痕迹：左右两半是同一张散点、标注不同。 | 代码卫生 | `p <- LabelPoints(plot = p, points = top_features, repel = TRUE)`，单联。 |
| SC-09 | 8 | `sc.pl.pca_variance_ratio(adata, log=True,,save="scanpy_PCA_ratio.png")` 的双逗号是 `SyntaxError`。Python 以整单元格编译，于是该单元里**什么都没有执行**：不只是比值图，连 `sc.tl.pca()`、`sc.pp.neighbors()`、`sc.tl.umap()` 全都没跑。后续单元连锁失败于 `KeyError: 'No "neighbors" in .uns'`。 | 致命 | 去掉一个逗号。 |
| SC-10 | 2, 8 | 两侧用不同阈值过滤同一批 2,700 个细胞（`nFeature_RNA < 2500` vs `n_genes_by_counts < 2000`），于是分别描述 2,638 与 2,635 个细胞、9 与 8 个群。正文与代码都没提这个不一致。 | 改变结果 | 统一阈值；§12 显示 Python 侧改用 2,500 恰好得到 2,638 个细胞。 |
| SC-11 | 8 | `sc.pp.scale(adata, max_value = 10)` 在 cell 64 末尾调一次、cell 67 开头又调一次。第二次是对已截断数据的重复截断（本数据上是空操作，但一旦 `max_value` 不同就会静默改变数值）。 | 代码卫生 | 只保留一次。 |
| SC-12 | 8 | `rename_categories("leiden", new_cluster_names)` 要求每个类别正好一个名字，数量不符就抛难以理解的错误；8 个名字对 9 个类别时整单元崩掉。更关键的是这些名字是**位置式**的——按类别顺序而非 marker 证据分配——而且**是错的**：名叫 `CD8 T` 的群有 78.1% 是 Seurat 的 Naive CD4 T（CD8A 0.14、CD8B 0.24、IL7R 1.61），名叫 `NK` 的群是 61.8% CD8 T 加 33.3% NK。 | 改变结果 | 校验类别数，并由 `rank_genes_groups` 结果反推映射（§11.4）。 |
| SC-13 | 6, 8 | 两侧算完 marker 都不落盘；scanpy 的 `save=` 是"函数名 + 字符串"拼接，于是产生 `highest_expr_genesscanpy_highest_expr_genes.png`、`violinscanpy_violin_plot_QC.png`、`umapscanpy_final_umap.png`。课程随附的图之所以要改名，正是这个原因。 | 代码卫生 | 导出 `markers_all_clusters.csv`、`clusterN_markers.csv`、`cluster_label_audit*.csv`。 |
| SC-14 | 1, 8 | `genes.tsv` 有 95 个重复行、涉及 91 个 symbol（32,738 行、32,643 个唯一名，`7SK` 出现 4 次），另有 21 个 symbol 含 `_`。Seurat 把 `_` 改写成 `-`（"Feature names cannot have underscores, replacing with dashes"），scanpy 的 `var_names_make_unique()` 把重复项改成 `ACE` → `ACE-1`。两侧基因名双向不一致。 | 稳健性 | 记录两份映射；跨工具关联用 `gene_ids`。 |
| SC-15 | 1 | 数据目录是**旧版 Cell Ranger 1.1.0 结构**：`barcodes.tsv`、`genes.tsv`、`matrix.mtx`，没有 `features.tsv`、没有 `.gz`。`Read10X()` 与 `read_10x_mtx()` 都能读，但任何期待 `features.tsv.gz` 的新工具都读不了。 | 稳健性 | 本模块无需处理，记录在案。 |
| SC-16 | 6 | 调用了 `geom_text_repel()`，却从未 `library(ggrepel)`（只载入了 Seurat、dplyr、ggplot2）。整个 DEG 单元以 `could not find function "geom_text_repel"` 崩溃——**连带最后一行 `saveRDS(seurat_obj, "pbmc3k_final.rds")` 也不执行**。 | 致命 | Part 1 加 `library(ggrepel)`。 |
| SC-17 | 1 | `CreateSeuratObject()` 用默认 `min.cells = 0`，保留全部 32,738 个基因（其中 19,024 个 QC 后全零），而 scanpy 的 `filter_genes(min_cells = 3)` 只留 13,714 个。这不是无害差异：2,000 个 HVG 里有 30 个变化，Louvain 解随之移动（模块度 0.8584 vs 0.8728，群大小差 ±5%），教材那两份硬编码注释表也不再自洽（§13.2）。 | 改变结果 | 显式声明基因全集，不要静默改动。 |
| SC-18 | 7 | Part 7 的两张 UMAP **不可跨运行复现**：同一输入、同一代码跑两次，结构相关只有 0.5631（`UMAP_without_harmony`）与 0.5676（`UMAP_with_harmony`），与"对参考图的相关"（0.5664、0.5689）同量级。uwot 是多线程的，`RunUMAP()` 没有线程无关的随机种子。 | 稳健性 | 比对群体结构或嵌入矩阵，永不比像素；需要逐位可复现就设 `n.threads = 1`。 |
| SC-19 | 6 | 随附参考图 `Marker_Heatmap.png` 是 3,000 × 6,000 px，即 300 dpi 下 10 × 20 英寸；而教材 `ggsave(..., width = 10, height = 10)` 只会得到 3,000 × 3,000。参考图不是教材代码产出的。 | 代码卫生 | 无需处理，记录课程材料内部不一致。 |
| SC-20 | 7 | `RunHarmony(merged_seurat, group.by.vars = "orig.ident")` 省略了 `dims`，于是对 `RunPCA()` 得到的全部 50 个 PC 做整合（默认 `npcs = 50`）；而随后的 `RunUMAP(reduction = "harmony", dims = 1:10)` 只用前 10 个 harmony 分量。去批次与可视化跑在不同维度上。 | 稳健性 | 显式 `dims = 1:10`（或 `1:20`），并让 UMAP 使用同一组。 |

另有两点不是代码 bug，但必须记录：

- 教材的 HTML 里从未出现随机种子，而整个模块处处随机（Louvain、UMAP、harmony）；
- 本目录下 28 张参考 PNG 是教师 2025 年的输出，其中 13 张与 2024 版逐字节相同；另有 3 张（`UMAP_noPCA.png`、`UMAP_three_stages.png`、`immune_markers.png`）不是教材里任何代码产出的（§14）。

---

## 2. 环境与数据

### 2.1 环境

本模块不新建环境。教师共享环境里两套工具链都在：

```bash
source /lustre1/share/miniconda3/etc/profile.d/conda.sh
conda activate /lustre1/share/miniconda3/envs/scRNA_seq
```

| 组件 | 版本 | 备注 |
|---|---|---|
| R | 4.4.1 | |
| Seurat | 5.1.0 | SeuratObject 5.0.2；v5 分层对象 |
| harmony | 1.2.1 | |
| ggrepel | 0.9.6 | 装了，但教材从未载入（SC-16） |
| dplyr / ggplot2 / patchwork | 1.1.4 / 3.5.1 / 1.3.0 | |
| Matrix / uwot / leiden / igraph | 1.7.1 / 0.2.2 / 0.4.3.1 / 2.1.1 | |
| FNN / RANN | 1.1.4.1 / 2.6.2 | 仅供矫正版的混合度指标使用 |
| presto | **缺失** | 缺它 `FindAllMarkers()` 退回慢速 Wilcoxon 路径：一个单元 330 s |
| Python | 3.12.7 | |
| scanpy / anndata | 1.10.3 / 0.10.9 | |
| leidenalg / igraph | 0.11.6 | |
| numpy / pandas / seaborn | 2.0.2 / 2.2.3 / 0.13.2 | |
| Pillow | 11.0.0 | 供图像比对脚本使用 |

本模块在登录节点串行运行，与第 3、4、5 节课一致：总数据量 2,700 个细胞（Part 7 为 25,894 个），峰值内存约 4 GB，最长的单个单元是 `FindAllMarkers()` 的 5.5 分钟。没有比对步骤、不涉及调度器，故未加 `sbatch` 包装——与 RNA-seq、TBCR 两个模块记录的理由相同。

### 2.2 输入数据

| 路径 | 内容 |
|---|---|
| `/lustre1/share/scRNA_seq/hg19/` | `barcodes.tsv`（2,700）、`genes.tsv`（32,738）、`matrix.mtx`（2,286,884 个非零），Cell Ranger 1.1.0 / hg19 |
| `/lustre1/share/qinbowen/pbmc_stim.RData` | `ctrl.sparse`（35,635 × 13,019）、`stim.sparse`（35,635 × 12,875），稀疏 `dgCMatrix` |
| `/lustre1/share/scRNA_seq_student/qinbowen_2025/` | 教师输出，仅作比对参照 |

矩阵是标准 10x PBMC 3k：2,700 个细胞，中位 817 个基因、2,197 个 UMI，线粒体占比中位 2.03%。`genes.tsv` 有 95 个重复 symbol、21 个含 `_` 的 symbol（SC-14）。

### 2.3 运行目录

```
/lustre1/user/bjx131_pkuhpc/scrna_run_20261026/
├── hg19/                     计数矩阵副本（按教材 cell 7 的 rsync 指令）
├── pbmc_stim.RData           刺激组数据副本
├── rr/                       R 原样版（图直接写在该目录）
├── rr_fixed/                 R 矫正版（MIN_CELLS = 0）
├── rr_mincells3/             R 矫正版、基因全集对齐 scanpy（§13.2）
├── py/  py_asis/             Python 原样版
├── py_fixed/                 Python 矫正版（N_GENES_MAX = 2000，教材原值）
├── py_align2500/             Python 矫正版、阈值对齐 Seurat（§12.3）
├── logs/                     每次运行一份日志 + 一份逐单元状态表
└── scripts/                  单元文件 + 驱动器 + 分析脚本
```

教材的每个代码单元都按 HTML 中的原始编号抽出（R 为 11、13、16…，Python 为 54、55、57…），交给两个小驱动器执行：

- `scripts/run_cells.R`：在**同一个常驻 R 会话**里逐单元求值，逐单元捕获错误、警告、消息，记录耗时，并在单元结束后记录各对象状态；
- `scripts/run_cells.py`：同一思路，Python 常驻命名空间。

两个驱动器精确复刻 notebook 语义：报错的单元被记录，下一个单元继续执行；不自动修正、不重试、不重排。这正是把"教材有 bug"变成"25 张图里只有 18 张存在，另外 7 张死在这些单元"的关键。

### 2.4 仓库内的产物布局

```
scRNA-seq/
├── scRNA_seq_analyze_totorial_2026.html   课程材料（未改动）
├── scRNA_seq_lecture_2026.pdf             课程幻灯片（未改动）
├── *.png                                  28 张参考图（未改动）
├── ref/                                   本文档、中文版、结果导读
└── results/
    ├── 01_qc/                  质控图 + qc_summary.csv + qc_rederivation_rawmatrix.csv
    ├── 02_hvg_pca/             高变基因与 PCA 图（两侧、两版）
    ├── 03_clustering/          UMAP 图 + cluster_sizes.csv
    ├── 04_annotation/          marker/注释图、marker 表、审计表与对账表
    ├── 05_deg/                 火山图 + DEG 表 + scanpy 排名/点图
    ├── 06_integration/         harmony 图 + integration_metrics.csv + 分样本群占比表
    ├── 07_objects/             pbmc3k.h5ad + pbmc3k_obs.csv
    ├── 08_figure_comparison/   三张比对表 + 书面比对报告
    ├── 09_run_logs/            六次运行的逐单元日志与状态表
    └── scripts/                重新生成以上全部内容的脚本
```

`pbmc3k_final.rds` **不入库**：Seurat 为全部 32,738 个基因保留稠密 `scale.data`，文件 331 MiB。Part 6 里那行 `saveRDS()` 即可重建。

---

## 3. 原样运行：教材实际产出了什么

### 3.1 R 侧（Part 1-7）

| 单元 | 步骤 | 秒 | 结果 |
|---|---|---|---|
| 11 | `Read10X` + `CreateSeuratObject` | 7.2 | 32,738 × 2,700；警告：特征名 `_` → `-` |
| 13 | 线粒体 QC + 小提琴图 | 1.9 | 2,700 → 2,638 个细胞 |
| 16 | `NormalizeData` | 0.8 | |
| 18 | `FindVariableFeatures`（vst，2,000） | 4.2 | 警告：零离散度取 log-10 |
| 21 | `ScaleData`（全部 32,738 基因） | 5.8 | 对象 0.73 GB |
| 23 | `RunPCA` + 4 张图 | 8.7 | 50 个 PC |
| 27 | `FindNeighbors(1:10)` + `FindClusters(0.5)` | 1.3 | 9 个社区，模块度 0.8584 |
| 29 | `RunUMAP(1:10)` | 6.7 | 65 条消息（uwot 进度） |
| 33 | 9 个 marker 的 `FeaturePlot` | 3.0 | |
| 36 | `FindAllMarkers(only.pos)` + 热图 | **330.3** | 11,986 行；无 presto |
| 39 | 注释表（仅注释） | 0.0 | |
| 40 | 12 个 marker 的 `FeaturePlot` | 4.3 | |
| 43 | `RenameIdents` + 带注释 UMAP | 0.6 | 9 个标签 |
| 46 | DEG + 火山图 + `saveRDS` | 11.3 | **报错：`could not find function "geom_text_repel"`** |
| 48 | merge + 未去批次 UMAP | 63.1 | 35,635 × 25,894 |
| 51 | harmony + UMAP | 60.1 | 50 维 harmony 嵌入 |

合计 509 s，其中 `FindAllMarkers()` 占 65%。cell 46 在崩溃前已经算完两次 `FindMarkers()`（724 行与 650 行），所以 DEG 数字留在日志里，但火山图与 `pbmc3k_final.rds` 都没有。cell 51 的 harmony 嵌入是 50 维（SC-20）：`RunHarmony()` 默认用 `pca` 的全部主成分，而合并对象经 `RunPCA()` 得到 50 个。

### 3.2 Python 侧（Part 8）

| 单元 | 步骤 | 秒 | 结果 |
|---|---|---|---|
| 54 | import | 0.0 | |
| 55 | `read_10x_mtx` + `var_names_make_unique` | 0.5 | 2,700 × 32,738 |
| 57 | `highest_expr_genes` | 0.8 | 出图 |
| 59 | `filter_cells(200)` + `filter_genes(3)` + QC 小提琴图 | 1.1 | 过滤基因后 2,700 × 13,714 |
| 61 | 两张散点 + 阈值 + normalize + log1p | 0.4 | 2,635 个细胞 |
| 64 | HVG + 取子集 + scale | 0.9 | 1,833 个 HVG |
| 67 | PCA + 比值图 + neighbors + UMAP | 0.0 | **SyntaxError（SC-09）——该单元一行未执行** |
| 70 | leiden + UMAP + rank_genes_groups | 0.0 | `KeyError: 'No "neighbors" in .uns'` |
| 72 | dotplot | 0.0 | `ValueError: groupby has to be a valid observation` |
| 75 | 改名 + 终版 UMAP + 写 h5ad | 0.0 | `ValueError: leiden is neither in .obs nor in .var` |

合计 3.8 s、5 张图。没有 `pbmc3k.h5ad`。

### 3.3 结论

| | 产出 | 缺失 |
|---|---|---|
| Seurat 图 | 14 张里的 13 张 | `Naive_Memory_CD4_DEG.png` |
| scanpy 图 | 11 张里的 5 张 | PCA、PCA 比值、markers、rank、rank 小提琴、dotplot、终版 UMAP |
| 对象 | 2 个里的 0 个 | `pbmc3k_final.rds`、`pbmc3k.h5ad` |

**25 张图里的 18 张、2 个对象里的 0 个。** 修正版（§13）产出 25/25，外加教材本想画却从未画出的比值图，以及两个对象。

### 3.4 与 28 张参考图的逐图像素比对

每张重生成的图都比对了像素尺寸、着色数量与 64 × 64 归一化结构相关（`results/08_figure_comparison/compare_aswritten_vs_reference.tsv`）。

| 图 | 像素（我们 / 参考） | 相关度 |
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
| scanpy：highest-expr、HVG、2 张散点、QC 小提琴 | 尺寸完全一致 | 0.9992 - 1.0000 |

PC 热图上的 1.0000 不是巧合：同一份 PCA、每个主成分同样取 100 个平衡细胞、完全确定的绘图路径，PNG 逐字节级可复现。热图的 0.7250 纯粹是长宽比效应（SC-19）。Part 7 的两张 UMAP 是唯一"真实差异"，而这点差异落在跨运行噪声里（SC-18：我们自己两次运行之间是 0.5631，对参考图是 0.5664）。

---

## 4. Part 1 —— 载入数据

### 4.1 教材做了什么

`Read10X()` 把旧版 Cell Ranger 的三个文件读成 32,738 × 2,700 的稀疏矩阵，`CreateSeuratObject()` 包成对象。出现两条警告，两条都有信息量：

```
Feature names cannot have underscores ('_'), replacing with dashes ('-')
```

`genes.tsv` 里有 21 个 symbol 含 `_`（`RP11-442N24__B.1`、`Metazoa_SRP` 等）。Seurat 静默改写，此后用原符号查找必然失败。Python 侧则对重复项做相反的处理：95 行重复、涉及 91 个符号（32,738 行、32,643 个唯一名），`var_names_make_unique()` 把第二个 `ACE` 变成 `ACE-1`。因此跨工具链关联基因必须用 `gene_ids`，不能用 symbol。

```
Default search for "data" layer in "RNA" assay yielded no results; utilizing "counts" layer instead.
```

这是 Seurat 5 的消息：`PercentageFeatureSet()` 要 `data` 层，但此刻还没 `NormalizeData()`。功能上无害。

### 4.2 实测

| | |
|---|---|
| 矩阵 | 32,738 个特征 × 2,700 个细胞，2,286,884 个非零 |
| 基因数 / UMI / 线粒体占比中位 | 817 / 2,197 / 2.03% |
| 基因数范围 | 212 - 3,422 |
| 线粒体占比范围 | 0 - 22.6% |

### 4.3 坑与调参

SC-14（符号改写、重复）与 SC-15（旧版布局、无 `features.tsv`）都在这一步。`Read10X()` 能读旧布局，`scanpy.read_10x_mtx()` 靠 `gex_only` 默认值也能读。

调参：`CreateSeuratObject(min.cells = 3)` 是把 R 的基因全集对齐 scanpy `filter_genes(min_cells = 3)` 的自然做法，但**不要随手改**——SC-17 量化了它的代价。

### 4.4 扩展

若数据来自 Cell Ranger 3+，用 `Read10X_h5()`；其它来源用 `ReadMtx()`；多模态数据要按 `features.tsv` 的类型列区分基因表达与其它特征。

---

## 5. Part 2 —— 质量控制

### 5.1 做什么、为什么

三条过滤，对应三种失败模式：基因数太少 = 空滴或死细胞；基因数太多 = 双细胞；线粒体占比高 = 应激或濒死细胞，胞质 mRNA 已外漏。

### 5.2 实测

| | |
|---|---|
| 规则 | `nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5` |
| 保留细胞 | 2,700 中的 **2,638**（97.7%） |
| 因基因数过少丢弃 | 0 |
| 因基因数过多丢弃 | 5 |
| 因线粒体占比 ≥ 5% 丢弃 | 57 |
| 质控后中位基因数 / UMI / 线粒体占比 | 820 / 2,214 / 2.01%（线粒体最高 4.99%） |

scanpy 侧用 `n_genes < 2000` 保留 **2,635** 个细胞（SC-10）。在原始矩阵上按同样数字独立重算得到 2,634——差的这一个细胞，是线粒体占比在 scanpy 里按**过滤后**的 13,714 基因矩阵计算、而 Seurat 用全部 32,738 个特征，边界细胞的占比因此落在 5% 两侧。记录在 `results/01_qc/qc_rederivation_rawmatrix.csv`。

### 5.3 坑

实质性的是阈值不一致（SC-10）：它就是 §12 里"9 群 vs 8 群"的源头。另需注意 `^MT-` 只匹配实际存在的 13 个基因；换成 `^MT\.`（MitoCarta 风格）或小鼠数据（`^mt-`）要改正则，而一旦不匹配，`PercentageFeatureSet()` 的表现会被误读成数据问题。

### 5.4 调参

`nFeature` 上限与数据本身强相关：PBMC 3k 是 2016 年的浅测序（中位 817 个基因），2,500 很宽松；现代 10x 数据中位常超 2,000，同样阈值会删掉大部分细胞。定阈值前先看 `nFeature` vs `nCount`、`nFeature` vs `mt%`（教材 cell 61 的两张散点正是为此）。

---

## 6. Part 3 —— 高变基因

`FindVariableFeatures(selection.method = "vst", nfeatures = 2,000)` 拟合"方差 ~ 平均表达"的 loess 曲线，按标准化残差排序；从 32,738 个特征里选出 2,000 个，出 `VariableFeaturePlot.png`。

两点说明：

- `log-10 transformation introduced infinite values` 是预期警告：`vst` 对离散度取 `log10`，零方差基因给出 `-Inf`，Seurat 会丢弃，消息仅为提示。
- `LabelPoints()` 标注方差最高的 20 个基因：PPBP、LYZ、S100A9、IGLL5、FTL、GNLY、PF4、FTH1、GNG11、S100A8、FCER1A、HLA-DRA、CD74、CLU、NKG7、CST3、GZMB、CCL4、C1QA、HLA-DPB1，全是白细胞 marker，这本身是"方法有效"的体检。

SC-08 在此：`p <- p + plot2` 画出的是**两联图**，参考图带着这个痕迹。`results/02_hvg_pca/` 两个版本都在：`VariableFeaturePlot.png`（原样、两联、对参考 0.9751）与 `VariableFeaturePlot_corrected.png`（修正、单联）。

---

## 7. Part 4 —— 标准化与 PCA

`ScaleData(features = rownames(seurat_obj))` 对全部 32,738 个特征做中心化与缩放（SC-01），随后 `RunPCA(features = VariableFeatures(...))` 只用 2,000 个 HVG，返回 50 个主成分。

| | |
|---|---|
| scale 的特征数 | 32,738（其中 13,714 非空） |
| 对象体积 | 0.73 GB |
| `ScaleData` 耗时 | 4.2 s |
| 主成分 | 50 |
| 出图 | `PCA_VizDimLoadings.png`、`PCA.png`（肘图，20 维）、`PCA_dimplot.png`、`PCA_heatmap_ALL.png`（6 张热图 × 100 细胞） |

`PCHeatmap()` 是 `DimHeatmap()` 的废弃别名；Seurat 5.1.0 里仍可用且不报警告，很容易一用多年（SC-07）。肘图呈 PBMC 3k 的经典形状——拐点在 PC 7-10 之间——这正是 `dims = 1:10` 的依据；但肘图画的是 20 维，而正文从未解释为什么取 10。

`PCA_heatmap_ALL.png` 同时也是"我们的 PCA 与教师完全相同"的证据：结构相关 1.0000。

---

## 8. Part 5 —— 聚类与 UMAP

```
FindNeighbors(seurat_obj, dims = 1:10)      # 10 个 PC 上的 SNN 图
FindClusters(seurat_obj, resolution = 0.5)  # Louvain
RunUMAP(seurat_obj, dims = 1:10)
```

| | |
|---|---|
| 图 | 2,638 个节点、95,650 条边 |
| Louvain | 模块度 0.8584，**9 个社区** |
| 群大小 | 657 / 503 / 471 / 339 / 307 / 175 / 144 / 28 / 14 |
| UMAP | `UMAP.png`；关于 uwot 取代 Python UMAP 的警告 |

群大小本身已经回答了"这是不是真实生物学拆分"：657 与 503 是两个 CD4 腔室，471 单核、339 B、307 CD8、175 FCGR3A+ 单核、144 NK，最后两个小群——28 个树突状细胞、14 个血小板。血小板群只占 0.5%，在 `resolution = 0.5` 下能保留，0.3 就没了。

---

## 9. Part 6 —— 注释与差异表达

### 9.1 两份硬编码映射互相矛盾（SC-02）

教材把 cluster → 细胞类型映射写了两遍：

| cluster | cell 39 的表 | cell 43 的 `new.cluster.ids` | marker 说的是什么 |
|:--------|:---------------------|:--------------------------------------|:------------------------------------------|
| 0 | Naive CD4 T | Naive CD4 T | CD4（LDHB、IL7R、LTB、CD3D、IL32） |
| 1 | **CD14+ Mono** | **Memory CD4 T** | **CD4**（CCR7、NOSIP、LEF1、FHIT、NELL2） |
| 2 | **Memory CD4 T** | **CD14+ Mono** | **单核**（LGALS2、S100A8、FCN1、S100A9、CD14） |
| 3 | B | B | B（CD79A、MS4A1、CD79B） |
| 4 | CD8 T | CD8 T | CD8（GZMK、CCL5、GZMA、CST7） |
| 5 | FCGR3A+ Mono | FCGR3A+ Mono | FCGR3A+ 单核（CDKN1C、HES4、MS4A7、FCGR3A） |
| 6 | NK | NK | NK（GZMB、FGFBP2、SPON2、PRF1） |
| 7 | DC | DC | DC（SERPINF1、FCER1A、CLEC10A） |
| 8 | Platelet | Platelet | 血小板（GP9、ITGA2B、LY6G6F） |

本次运行里**向量对、表错**：把 cell 39 的表抄进 `new.cluster.ids` 的学生，会把 503 个 CD4 T 细胞标成单核、把 471 个单核标成 CD4 T。两份输出都在 `results/04_annotation/cluster_label_audit.csv`，其中"观测到的程序"由 marker 面板（`marker_panel_mean_expression_seurat.csv`）给出。

用同样的方法检查最终对象的标签，9 个全对：

| 标签 | 细胞数 | 证据（log 标准化均值） |
|---|---|---|
| Naive CD4 T | 657 | IL7R 1.70、CCR7 0.58、LEF1 0.57、CD14 0.03 |
| Memory CD4 T | 503 | IL7R 1.19、CCR7 0.84、CD14 0.03 |
| CD14+ Mono | 471 | CD14 1.53、LYZ 5.09、S100A8 3.91 |
| B | 339 | MS4A1 2.17、CD79A 2.80 |
| CD8 T | 307 | CD8A 1.01、CD8B 0.79、NKG7 3.36 |
| FCGR3A+ Mono | 175 | FCGR3A 2.75、MS4A7 1.71 |
| NK | 144 | GNLY 4.17、GZMB 3.68、KLRF1 0.52 |
| DC | 28 | FCER1A 2.11、CLEC10A 1.28 |
| Platelet | 14 | PPBP 5.85、PF4 4.93 |

对得**靠运气**：§13.2 显示，只要 2,000 个 HVG 里变 30 个（1.5%），cluster 1 与 2 就互换编号，同一个固定向量立刻变成 957 个细胞标错。

### 9.2 Marker 与热图

| | |
|---|---|
| `FindAllMarkers(only.pos = TRUE)` | 9 个群共 11,986 行，330 s（无 presto） |
| 通过 `p_val_adj < 0.05 & avg_log2FC > 1` | 3,483 行 |
| 进热图 | 90 个（每群 10 个） |
| 每群显著数 | 89（cluster 0）… 753（cluster 5） |

SC-03：过滤算完就丢，`top10` 来自未过滤表，`slice_max(avg_log2FC)` 完全可以把 `p_val_adj = 0.2` 的基因选进去。修正版先过滤（3,483 行）再每群取 10 个；由此得到的热图与参考图相关 0.9865。

### 9.3 DEG 单元崩了（SC-16）

`geom_text_repel()` 需要 `library(ggrepel)`，而教材只载入了 Seurat、dplyr、ggplot2。于是这个单元在**最贵的部分之后**失败：

| | |
|---|---|
| `FindMarkers(ident.1 = "Naive CD4 T")` | 724 个基因，353 个 `p_val_adj < 0.05`，其中 199 个同时 `|log2FC| > 0.5` |
| `FindMarkers("Memory CD4 T", "Naive CD4 T")` | 650 个基因，129 个 `p_val_adj < 0.05` |
| 报错 | `could not find function "geom_text_repel"` |
| 后果 | 没有 `Naive_Memory_CD4_DEG.png`，**也没有 `pbmc3k_final.rds`** |

火山图按 `p_val_adj` 排前五的是 LDHB（1.2e-115）、IL7R（1.8e-102）、LTB（1.9e-97）、CD3D（2.0e-91）、IL32（7.7e-86）——即"naive CD4 T 相对其余细胞"主要是 naive marker 更高、单核/B/NK 基因更低，这正是 `avg_log2FC > 0` 的含义。

SC-04：`-log10(p_val_adj)` 在本数据上是有界值（最小 1.24e-115），但离 `Inf` 只差一次舍入；两个基因列表算完不用；两张表都不落盘。

---

## 10. Part 7 —— 去批次整合

### 10.1 教材做了什么

载入 `pbmc_stim.RData`，建两个 Seurat 对象，设 `orig.ident`，`merge()`，然后**完全不质控**地 `NormalizeData` → `FindVariableFeatures` → `ScaleData` → `RunPCA` → `RunUMAP(1:10)` → `DimPlot`（SC-06），再 `RunHarmony(group.by.vars = "orig.ident")` → `RunUMAP(reduction = "harmony", dims = 1:10)` → `DimPlot`。

### 10.2 实测

| | |
|---|---|
| ctrl / stim 细胞数 | 13,019 / 12,875（35,635 个特征） |
| 合并后 | 25,894 个细胞 |
| 教材省略的质控 | 12,932 / 12,779 → **25,711**（−183，−0.7%） |
| 去批次前 UMAP | `UMAP_without_harmony.png` |
| 去批次后 UMAP | `UMAP_with_harmony.png` |
| harmony 输入 | **全部 50 个 PC**（SC-20），UMAP 只用前 10 个 |

### 10.3 把"只能看图"的部分量化

修正版补了一个不需要聚类的指标：对每个细胞，取 30 个最近邻，算其中来自**另一样本**的比例（0 = 完全分家，0.5 = 完全混合）。

| 空间 | 跨样本近邻比例 |
|---|---|
| PCA（10 个 PC），未去批次 | **0.0287** |
| harmony（10 维） | **0.4051** |

再把 harmony 嵌入按分辨率 0.5 聚类：

| | |
|---|---|
| 群数 | 12 |
| 单样本占比 > 90% 的群 | **0** |
| ctrl 占比偏离 0.5 的平均绝对值 | 0.030 |

最大的群 6,065 个细胞（3,096 ctrl / 2,969 stim），最小的 92 个（38 / 54）。这就是"harmony 起作用了"的数字版：跨样本近邻比例从 2.9% 升到 40.5%（去批次前，每个细胞的 30 个近邻里有 97% 来自同一样本），而每个群内部两个样本仍然平衡。

SC-18：两张 UMAP 不可跨运行复现（同一输入两次运行之间 0.5631 / 0.5676），所以本节比较的是嵌入与群组成，不是像素。

---

## 11. Part 8 —— scanpy

### 11.1 哪些活下来、哪些没有

cell 54-64 用 3.8 s 跑完，产出 5 张图。cell 67 里是那个双逗号（SC-09）。由于 Python 以整单元格编译，`SyntaxError` 抹掉了整个单元：PCA 从未计算，`neighbors()` 无法建图，`umap()` 无法运行，cell 70、72、75 连锁失败（`KeyError: 'No "neighbors" in .uns'`、`groupby has to be a valid observation`、`leiden is neither in .obs nor in .var`）。

这也解释了为什么课程参考图里没有 `pca_variance_ratio`：教师自己那一遍撞在同一堵墙上，被缩短的文件名清单也证明只产出了 5 张 scanpy 图。

### 11.2 Python 侧确实走到的数字

| | |
|---|---|
| `filter_cells(200)` + `filter_genes(3)` 之后 | 2,700 × **13,714** 个基因 |
| `n_genes < 2000 & pct_mt < 5` 之后 | **2,635** 个细胞 |
| HVG（`min_mean 0.0125, max_mean 3, min_disp 0.5`） | **1,833** 个 |
| PCA | 50 个主成分，`arpack` |
| 近邻 / UMAP | 10 / 40 个 PC |

修正版把剩下的补完：分辨率 0.9 的 8 个 leiden 社区（534 / 643 / 340 / 497 / 432 / 142 / 35 / 12 个细胞）、`rank_genes_groups`（Wilcoxon）、点图，以及 `pbmc3k.h5ad`。

### 11.3 位置式命名是错的（SC-12）

按教材阈值，8 个名字正好对得上——其中两个落在了错的群上。证据来自 `rank_genes_groups` 与跨工具对账表（§12.2）：

| scanpy 群 | 细胞数 | 分配的名字 | 主要 marker | Seurat 组成 | 判定 |
|:----------|:-------|:--------------|:---------------------------|:------------------------------|:----------------|
| 0 | 534 | CD4 T | RPS12、RPS27、RPL32 | 70.6% Memory CD4 T、26.2% Naive CD4 T | 可接受（确为 CD4 群） |
| 1 | 643 | **CD8 T** | LDHB、LTB、CD3D、IL32、TPT1 | **78.1% Naive CD4 T**、18.2% Memory CD4 T | **错**（CD8A 0.14、CD8B 0.24） |
| 2 | 340 | B | CD74、CD79A、HLA-DRA | 99.1% B | 对 |
| 3 | 497 | CD14+ Monocytes | LYZ、S100A9、S100A8 | 92.4% CD14+ Mono | 对 |
| 4 | 432 | NK | NKG7、CST7、GZMA、CTSW | **61.8% CD8 T**、33.3% NK | 部分错 |
| 5 | 142 | FCGR3A+ Monocytes | LST1、FCGR3A、FCER1G | 97.2% FCGR3A+ Mono | 对 |
| 6 | 35 | Dendritic | HLA-DPA1、HLA-DPB1、HLA-DRA | 77.1% DC | 对 |
| 7 | 12 | Platelet | PF4、PPBP、GNG11 | 100% Platelet | 对 |

`rename_categories()` 按类别顺序（即"Louvain 先编号了哪个社区"）分配名字，与生物学无关。于是 8 个标签里有 2 个实质性错误，其中最大的一个（`CD8 T`，占 24%）正是终版 UMAP 里最大的群。

### 11.4 修正版怎么做

先断言类别数，再由 `rank_genes_groups` 的主要 marker 反推映射，打印审计表，最后才改名。按教材阈值会得到 8 个类别、改名照常进行；按对齐阈值（§12.3）会得到 9 个，修正版拒绝猜测，保留数字标签并写出 `cluster_label_audit_py.csv` 交给人判断。

---

## 12. 两套工具链对账

### 12.1 "9 群 vs 8 群"之谜

| | Seurat | scanpy（教材） | scanpy（对齐） |
|---|---|---|---|
| 细胞数 | 2,638 | 2,635 | **2,638** |
| 基因数 | 32,738 | 13,714 | 13,714 |
| HVG | 2,000（vst） | 1,833（离散度） | 1,838 |
| 使用的主成分 | 10 | 40 | 40 |
| 分辨率 | 0.5（Louvain） | 0.9（leiden） | 0.9（leiden） |
| 群数 | **9** | **8** | **9** |

"R 出 9 群、Python 出 8 群"常被解释成"算法不同"。这不是全部：把 Python 侧阈值从 2,000 提到 2,500，得到**与 Seurat 完全相同的细胞数**（2,638）和 **9 个群**。群数主要由质控阈值决定，而不是 Louvain 与 leiden 的差别。

两个 9 群解在组成上仍然不同：Seurat 把 CD4 拆成 naive/memory 并把 NK 与 CD8 分开，而对齐后的 scanpy 解里有一个核糖体基因偏高的 CD4 群和一个合并的 CD8/NK 群。

### 12.2 跨工具列联表

两侧共享 barcode，可以直接交叉（两个对象共有 2,635 个细胞；`results/04_annotation/cluster_correspondence_counts.csv`）：

| Seurat 标签 | CD4 T | CD8 T | B | CD14+ Mono | NK | FCGR3A+ Mono | Dendritic | Platelet |
|---|---|---|---|---|---|---|---|---|
| Naive CD4 T | 140 | 502 | 1 | 0 | 13 | 0 | 0 | 0 |
| Memory CD4 T | 377 | 117 | 1 | 0 | 8 | 0 | 0 | 0 |
| CD14+ Mono | 0 | 0 | 0 | 459 | 0 | 4 | 8 | 0 |
| B | 1 | 1 | 337 | 0 | 0 | 0 | 0 | 0 |
| CD8 T | 16 | 23 | 1 | 0 | 267 | 0 | 0 | 0 |
| FCGR3A+ Mono | 0 | 0 | 0 | 37 | 0 | 138 | 0 | 0 |
| NK | 0 | 0 | 0 | 0 | 144 | 0 | 0 | 0 |
| DC | 0 | 0 | 0 | 0 | 0 | 0 | 27 | 0 |
| Platelet | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 12 |

列是教材命名后的 scanpy 群。B、单核、DC、血小板近乎对角，两个大的 T 腔室则完全错位：scanpy 的"CD8 T"列里 78% 是 Seurat 的 naive CD4 T。

### 12.3 哪些能比、哪些不能比

**能比、且本模块比了**：质控后细胞数、群数、按 marker 面板给出的群组成、marker 排名、DEG 数量、整合指标。

**不能比**：群的**编号**（两种方法各有自己的编号规则）、基因全集（HVG 之前 32,738 vs 13,714）、分辨率参数（0.5 Louvain 不等于 0.9 leiden——这里恰好给出 9 与 8，是巧合而非标定）。

让两侧真正对齐有两个开关：`CreateSeuratObject(min.cells = 3)`（13,714 个基因，但见 SC-17）与 scanpy 侧 `n_genes < 2500`（2,638 个细胞、9 个群）。两个都做，两侧就是 2,638 个细胞；本模块没有这么做，因为教材的图是按它自己的取值标定的，而且 SC-17 证明改动基因全集并不免费。

---

## 13. 修正版

### 13.1 改了什么，就只改了什么

`results/scripts/scRNA_seq_R_cells_corrected.R` 与 `..._py_cells_corrected.py` 就是教材原文，每处改动都标 `# --- FIX ... (SC-xx)`、每处新增都标 `# --- ADD ...`。真正影响产出的改动：

| 改动 | 效果 |
|---|---|
| `library(ggrepel)`（SC-16） | DEG 单元跑通：火山图、DEG 表、`pbmc3k_final.rds` |
| 比值图那一个逗号（SC-09） | PCA/近邻/UMAP 得以运行；生成 `scanpy_PCA_ratio.png`（参考图里没有这一张） |
| `pmax(p_val_adj, .Machine$double.xmin)`（SC-04） | 火山图纵轴构造上有界 |
| 过滤后的 marker 表赋值（SC-03） | 热图用 90 个显著基因，而不是"只看 log2FC 的前 10" |
| 群数守卫 + 审计表（SC-02、SC-12） | 映射由 marker 反推并打印；编号漂移变成警告，而不是静默错标 |
| 整合前质控 + harmony 显式 `dims` + 之后聚类（SC-06、SC-20） | 25,711 个细胞、12 个群、混合度指标 |
| `load()` 载入私有环境（SC-05） | 不再污染全局环境 |
| 导出表格（SC-13） | 10 个 CSV，教材只留在内存里 |

修正版 R 运行**精确复现**原样数字：2,638 个细胞、9 个群、11,986 行 marker、724 与 650 行 DEG、合并 25,894 个细胞、模块度 0.8584。两者共有的 13 张图里有 9 张**逐字节相同**（`results/08_figure_comparison/compare_aswritten_vs_corrected.tsv`）；不同的 4 张是两处有意改动（`VariableFeaturePlot` 单联、`Marker_Heatmap` 加了显著性过滤）与 Part 7 的两张 UMAP（SC-18 的不可复现性）。

### 13.2 基因全集的对照实验（SC-17）

同一份修正脚本改用 `MIN_CELLS=3`（即 scanpy 的基因全集），细胞数不变，但：

| | `MIN_CELLS = 0`（教材） | `MIN_CELLS = 3` |
|---|---|---|
| 特征数 | 32,738 | 13,714 |
| 与另一设置共享的 HVG | — | 2,000 个里的 1,970 个（98.5%） |
| Louvain 模块度 | 0.8584 | 0.8728 |
| 群大小 | 657/503/471/339/307/175/144/28/14 | 684/481/476/344/291/162/155/32/13 |
| marker 行数 | 11,986 | 11,629 |
| naive-CD4 DEG 行数 | 724 | 733 |
| memory-vs-naive DEG 行数 | 650 | 887 |
| 最终标签 | 正确 | **cluster 1 与 2 互换** |

最后一行是关键。`MIN_CELLS = 3` 时 Louvain 把单核社区编号为 1、memory-CD4 编号为 2；教材那个固定 `new.cluster.ids` 于是把 "Memory CD4 T" 贴到 CD14 1.51 / LYZ 5.07 / IL7R 0.19 的群上，把 "CD14+ Mono" 贴到 CD14 0.03 / IL7R 1.69 的群上——**957 个细胞（占 36%）标签错误，而所有图看起来都正常**。只有审计表能抓到这件事。

### 13.3 修正版 Python 补了什么

| | |
|---|---|
| 细胞数 / 基因数 | 2,635 / 13,714（教材阈值）；2,638（对齐阈值，附录） |
| HVG | 1,833 / 1,838 |
| 群数 | 8（教材阈值）/ 9（对齐阈值） |
| 出图 | 12 张（参考集的 11 张 + `pca_variance_ratio`） |
| 表格 | 8 个 `clusterN_markers.csv`、`cluster_label_audit_py.csv`、`pbmc3k_obs.csv` |
| 对象 | `pbmc3k.h5ad`（2,635 个细胞，19.5 MiB，gzip） |

### 13.4 有意不做的修正

- **不给 Louvain/UMAP/harmony 设种子。** 教材没有种子、参考图也没有，设种子只会让我们的图更难对齐。
- **不统一分辨率。** 0.5 Louvain 与 0.9 leiden 就是教材的取值；把它们"等价化"既不可能，也会掩盖 §12 的结论。
- **不安装 presto。** `FindAllMarkers()` 无 presto 要 330 s；装上能降到秒级，但会改变 Seurat 的数值路径，因此这里记录代价而不改变计算。
- **不改写交付的 Seurat 对象标签。** 它的标签是对的（§9.1），改写只会破坏与参考图的可比性。

---

## 14. 教材范围之外的参考图

28 张参考 PNG 里有 3 张不是教材任何代码的产物，因此本模块没有尝试重画：

- `UMAP_noPCA.png` 与 `UMAP_PCA_vs_noPCA.png`——教师 2025 幻灯片里"有没有 PCA"的对比；
- `UMAP_three_stages.png`——讲义里用来说明流程三阶段的三联示意图；
- `immune_markers.png`——CST 免疫 marker 速查表，属于外部参考图，不是分析产物。

教师 2025 目录里还有 HIV/EBV 病毒定量素材（`virus_star_ref/`、`ebv.fasta`、`HIV.fasta`、`virus.fa/gtf`）。教材完全未涉及，故不在本模块范围。

另有一处不一致只记录、不修正：参考图 `Marker_Heatmap.png` 是 3,000 × 6,000 px，而教材 `ggsave()` 要求 10 × 10 英寸、300 dpi，只能得到 3,000 × 3,000（SC-19）。

---

## 15. 复现本模块

```bash
# 0. 环境（教师共享环境，无需新建）
source /lustre1/share/miniconda3/etc/profile.d/conda.sh
E=/lustre1/share/miniconda3/envs/scRNA_seq

# 1. 按教材 cell 7 的指令准备数据
RUN=/lustre1/user/$USER/scrna_run
mkdir -p $RUN/{rr,py,logs,scripts}
rsync -r /lustre1/share/scRNA_seq/hg19 $RUN/
rsync /lustre1/share/qinbowen/pbmc_stim.RData $RUN/

# 2. 原样运行（25 张图里的 18 张，两个单元致命）
cp results/scripts/*.R results/scripts/*.py $RUN/scripts/
cd $RUN/rr && $E/bin/Rscript $RUN/scripts/run_cells.R \
    $RUN/scripts/scRNA_seq_R_cells_aswritten.R $RUN/logs R_aswritten
cd $RUN/py && MPLBACKEND=Agg $E/bin/python $RUN/scripts/run_cells.py \
    $RUN/scripts/scRNA_seq_py_cells_aswritten.py $RUN/logs Py_aswritten

# 3. 修正运行（25/25 张图 + 两个对象）
mkdir -p $RUN/rr_fixed && cd $RUN/rr_fixed && MIN_CELLS=0 $E/bin/Rscript \
    $RUN/scripts/run_cells.R $RUN/scripts/scRNA_seq_R_cells_corrected.R $RUN/logs R_corrected
mkdir -p $RUN/py_fixed && cd $RUN/py_fixed && N_GENES_MAX=2000 MPLBACKEND=Agg \
    $E/bin/python $RUN/scripts/run_cells.py \
    $RUN/scripts/scRNA_seq_py_cells_corrected.py $RUN/logs Py_corrected

# 4. 附录运行（基因全集与阈值对齐）
mkdir -p $RUN/rr_mincells3 && cd $RUN/rr_mincells3 && MIN_CELLS=3 $E/bin/Rscript \
    $RUN/scripts/run_cells.R $RUN/scripts/scRNA_seq_R_cells_corrected.R $RUN/logs R_mincells3
mkdir -p $RUN/py_align2500 && cd $RUN/py_align2500 && N_GENES_MAX=2500 MPLBACKEND=Agg \
    $E/bin/python $RUN/scripts/run_cells.py \
    $RUN/scripts/scRNA_seq_py_cells_corrected.py $RUN/logs Py_align2500

# 5. 表格与图像比对
cd $RUN && $E/bin/Rscript scripts/table_qc_clusters.R $RUN ../class/scRNA-seq/results
$E/bin/python scripts/compare_figs.py $RUN ../class/scRNA-seq \
    results/08_figure_comparison/compare_corrected_vs_reference.tsv rr_fixed py_fixed
```

教材自己的占位符（`data_dir <- "yourfolder/hg19"`、`load("pbmc_stim.RData")`、`read_10x_mtx("yourfolder/hg19")`）是原样版单元文件里唯一的改动，每个文件头部都用 `#@@NOTE` 标注。其余部分逐字照抄，包括那个双逗号。

---

## 16. 文件清单

| 路径 | 内容 |
|---|---|
| `ref/scRNA_seq_analyze_corrected.md` | 本文档（英文） |
| `ref/scRNA_seq_analyze_corrected_zh.md` | 中文版 |
| `ref/scRNA_seq_analyze_results_guide.md` | `results/` 导读 + 练习 |
| `results/scripts/scRNA_seq_R_cells_aswritten.R` | 教材 R 单元，除两处路径外逐字照抄 |
| `results/scripts/scRNA_seq_py_cells_aswritten.py` | 教材 Python 单元，逐字照抄（含 `SyntaxError`） |
| `results/scripts/scRNA_seq_R_cells_corrected.R` | 修正版 R 单元，每处改动有标注 |
| `results/scripts/scRNA_seq_py_cells_corrected.py` | 修正版 Python 单元，每处改动有标注 |
| `results/scripts/run_cells.R`、`run_cells.py` | notebook 语义驱动器（单会话、逐单元状态） |
| `results/scripts/table_qc_clusters.R` | 质控账目与群大小表 |
| `results/scripts/marker_panel_audit.{R,py}` | 两侧的 marker 面板逐群均值 |
| `results/scripts/annotation_audit.R` | 两份硬编码映射 vs 观测程序 |
| `results/scripts/cluster_correspondence.{R,py}` | 跨工具列联表 |
| `results/scripts/compare_figs.py` | 与 28 张参考图的图像比对 |
| `results/09_run_logs/*_cell_log.txt` | 逐单元的代码、标准输出、警告、报错、对象状态 |
| `results/09_run_logs/*_cell_summary.tsv` | 每个单元一行：状态、秒数、警告/报错数 |
| `results/08_figure_comparison/*.tsv` | 三张机器可读比对表 |
| `results/08_figure_comparison/figure_comparison.md` | 同一批数据的书面比对与解读 |
| `results/07_objects/pbmc3k.h5ad` | Part 8 的交付对象（19.5 MiB） |
| `results/07_objects/pbmc3k_obs.csv` | 逐细胞元数据：分群、质控指标 |
