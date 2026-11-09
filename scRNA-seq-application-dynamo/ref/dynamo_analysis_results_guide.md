# Dynamo 分析结果导读

本文件导读 `scRNA-seq-application-dynamo/results/` 下的全部产物，并给出可以亲手复现的练习。
配套的矫正增强版文档在 `ref/scRNAseq_dynamo_corrected.md`（英文）与
`ref/scRNAseq_dynamo_corrected_zh.md`（中文）；教材原件的 cell 抽取结果在
`results/scripts/dynamo_cells_aswritten.py`。

## 1. 怎么读这个目录

```
results/
├── 01_preprocess/          步骤 1-2：对象自身审计表 + 每群细胞数 + 逐基因标记
├── 02_velocity/            步骤 3-4：流线图、逐细胞箭头图、γ/β/半衰期、速度模长
├── 03_vectorfield/         步骤 5：四联图、地形图、状态图；逐细胞场量；逐群统计与检验；
│                           状态转移矩阵（通量 + 平均时间 + 行和）
├── 04_jacobian/            步骤 6-7：四个 Jacobian 面板、基因级张量与其成对/逐群表、
│                           加速度排序 top20/top200
├── 05_enrichment/          步骤 8：GO BP 2021 全表（925 个 term）、barplot/dotplot、
│                           gseapy 原始报告（txt/pdf）
├── 06_corrected/           矫正链：07/08/09 面板（教材只显示不存的三张）、curl 图、
│                           GO 改标图、UMAP 重算对比、Ductal vs 其余统计
├── 07_figure_comparison/   与两套教师参考图的机器可读对拍表
├── 09_run_logs/            原样跑的完整控制台日志
└── scripts/                重新生成以上全部内容的六个脚本
```

约定：

- `01_`–`05_` 是**原样跑**（`dynamo_cells_aswritten.py`）的产物，文件名沿用教材单元里的
  `save_kwargs["prefix"]`：`01_velocity_streamline_dyn_savefig.png` 对应第 14 单元，
  `06_reg_Ngn3_activate_Pax4_dyn_savefig.png` 对应第 22 单元。凡形如
  `inline_cellNN_auto_MM.png` 的，是**教材只内联显示、没存盘**的那几张图，由 harness 在
  `plt.show()` 处截下（第 23 单元两张、第 24 单元一张、第 32 单元两张）。
- `06_corrected/` 是**矫正跑**（`dynamo_corrected.py`）的产物，文件名保留原前缀并加
  `_corrected` 或独立编号，以免与原样版混淆。
- `03_vectorfield/*.tsv` 与 `04_jacobian/jacobian_*` 来自两次**数值导出**运行
  （`dynamo_vf_percell_export.py`、`dynamo_jacobian_export.py`）：它们把教材只画不存的量落成表。
- 表格一律 UTF-8、TAB 分隔；第一列常是索引（细胞名或群名），不要当数据类型混用。
- 六张原样图与矫正跑**逐字节相同**（md5 一致），所以两张目录下的同名图不必比较,本来就是同一份。

## 2. 01_preprocess：对象本身

`adata_summary.tsv` 是这一讲的"户口本"，14 行，每行一个数：

| 项 | 值 | 意思 |
|---|---|---|
| `cells_raw` / `genes_raw` | 3696 / 27998 | 读进来时 |
| `genes_pass_basic_filter` | **6956** | `shared_count=20` 后仍被考虑的基因数 |
| `hvg_use_for_pca` | **4000** | dynamo 用 SVR 选出的、真正参与 PCA/速度的基因 |
| `legacy_hvg_true` | 4004 | 文件自带的 `highly_variable_genes` 里 True 的个数 |
| `legacy_hvg_nan` | **10656** | 该列里 NaN 的个数（= 27998 − 17342，即上游另一套过滤丢掉的基因） |
| `overlap_legacy_vs_dynamo` | 1339 | 两个集合的交集 |
| `jaccard_legacy_vs_dynamo` | **0.2009** | 1339 / 6665 |
| `pancreas_genes_absent` | 0 | 17 个锚基因全部在 |

`jaccard = 0.20` 是本目录最重要的一行：dynamo 的 4000 个特征与"这文件自称的高变基因"只共享三分之一。
任何"dynamo 的速度 vs scanpy 的聚类"式比较，都是在两套特征空间之间做。

`var_flags.tsv`（27998 行 × 5 列）是逐基因的四个布尔标记：`use_for_pca`、`pass_basic_filter`、
`legacy_hvg`（文件自带）、`pancreas_gene`。想查某个基因为什么没进速度模型，grep 这一行。

`cluster_counts.tsv` 是八群规模：Ductal 916、Ngn3 high EP 642、Pre-endocrine 592、Beta 591、
Alpha 481、Ngn3 low EP 262、Epsilon 142、Delta 70。

**练习**：`awk -F'\t' '$2=="False" && $3=="True"' var_flags.tsv | wc -l` —— 有多少基因通过了基础过滤却没被选进 PCA？
（答：6956 − 4000 + 17 的重叠，实际输出会把答案给你。）

## 3. 02_velocity：速度与场

两张图是全课的门面：

- `01_velocity_streamline_dyn_savefig.png`（cell 14，1464×1023）：核平滑后的**流线**。
- `02_cell_wise_velocity_dyn_savefig.png`（cell 15，1222×816）：逐细胞**箭头**。

两者用的是同一个投影速度，差别只在"把箭头平均成场"这一步。对拍时这两张的墨迹 IoU 是 0.996 与 1.000 ——
它们也是能证明"速度场本身一致"的那两张：场若有差异，会先在这里暴露。

`gamma_summary.tsv` 记录动力学拟合：γ 中位 **0.382**（均值 0.953、q25 0.128、q75 0.908、max 64.58），
`beta` 的 min = max = **1.0**（`model="stochastic"` 固定 β，不是估计出来的），`half_life` = `ln2/γ`
（中位 1.814）。γ 有 3999 个有效值——4000 个基因里有一个拟合失败。

`cell_velocity_stats.tsv` 是 `velocity_umap_norm` 的分布（中位 0.375、q25 0.184、q75 0.597、max 1.12）。
模长是投影的结果，没有物理量纲；能解读的只有方向。

## 4. 03_vectorfield：场的三个导数与状态图

三张图：

- `03_combined_plots.png`（cell 18）：手工拼的四联图（潜势/速度/散度/加速度）。**注意这张没给 dpi**，
  1000×600；教师 2026 重跑版同面板是 3099×1834。
- `04_cell_state_dyn_savefig.png`（cell 19）：`topography`，潜势的地形图。
- `05_cell_state_transition_map_dyn_savefig.png`（cell 20）：`pl.state_graph`，群间转移画在 UMAP 上。

四张统计表里，`vf_per_cell.tsv` 是根：

| 列 | 含义 |
|---|---|
| `speed_pca` | ‖f(x)‖，单位 PCA 空间的场速 |
| `divergence_pca` / `divergence_umap` | ∇·f，**两个基各一份**，量级不同（±0.14 vs ±1.3），别混用 |
| `acceleration_pca` | ‖(f·∇)f‖ |
| `curl_umap` | ∇×f（教材只算不画，矫正版补了图） |
| `pca_ddhodge_potential` | ddhodge 潜势（类伪时间，低=早） |
| `velocity_umap_norm` | 投影速度模长，可与 §3 的 `cell_velocity_stats.tsv` 对照 |

`vf_stats.tsv` 是上表各列的全域分布（min/q25/median/mean/q75/max）；
`vf_by_cluster.tsv`（长表：cluster × quantity × mean/median/n）与
`vf_cluster_sign_summary.tsv`（宽表：只挑散度/加速度/旋度的均值与"正比例"）是逐群汇总。

`vf_cluster_tests.tsv` 是**推断**表，一行一个量：

| 列 | 含义 |
|---|---|
| `kruskal_H` / `kruskal_p` | 八群 Kruskal–Wallis |
| `ductal_mean` / `rest_mean` | Ductal(916) vs 其余(2780) 的均值 |
| `mannwhitney_p_ductal_vs_rest` | 两组 Mann–Whitney p |
| `frac_negative_ductal` / `frac_negative_rest` | 负值细胞比例 |

结论速记：五个量的 KW p 全 < 1e-90；但 Ductal 的散度**显著更不负**（−0.0140 vs −0.0263，p = 6.5e-06），
所以"Ductal 是鞍点"这个说法只在"7/8 个群都为负"的平凡意义上成立。真正干净的是 Ngn3 high EP：
唯一平均散度为正（+0.0229）、85 % 细胞为正、加速度最高（0.0158）；以及终末群的反面：
Beta −0.0682、Alpha −0.0502、Delta −0.0272，三者的正散度细胞比例都是 **0 %**。

状态图的两张矩阵：

- `clusters_graph_group_graph.tsv`（8×8，行=起点，列=终点）：**不是行随机矩阵**
  （Ductal 行和 1.76、Ngn3 low 0.94、Epsilon 1.57），对角项是驻留项，读作**通量**而不是转移概率。
  有意义的非对角项：Ductal→Alpha 0.68、Ductal→Ngn3 high 0.08、Pre-endocrine→Alpha 0.22、
  Delta→Pre-endocrine 0.27；Beta/Alpha 是完全吸收态（自 1.00、零流出）。
- `clusters_graph_group_avg_time.tsv`（8×8）：平均转移时间，单位是积分场上的**伪时间**，
  只能表内比较，不是小时也不是胚胎天数。最短：Delta→Pre-endocrine 16.0、Ductal→Alpha 24.6。
- 两个 `*_rowsum_*.tsv` 是行和，用来一眼看出"哪一行的通量溢出了"。

**练习**：把 `vf_per_cell.tsv` 按 `clusters` 分组，算 `pca_ddhodge_potential` 的均值并排序，
你会得到 Ductal(0.21) → Ngn3 low(0.26) → Ngn3 high(0.98) → Pre-endocrine(1.81) → Epsilon(1.99) →
Delta(2.12) → Alpha(2.13) → Beta(2.26)：**不用任何时间信息**恢复出的发育顺序。

## 5. 04_jacobian：谁调控谁

四个面板（`06_` 存盘，`07_/08_/09_` 是教材只内联显示的三张，原样跑时由 harness 截下）：

| 图 | 边 | 读法 |
|---|---|---|
| `06_reg_Ngn3_activate_Pax4_dyn_savefig.png` | Neurog3 → Pax4 | 红 = 激活；红色集中在 Ngn3-high/pre-endocrine |
| `inline_cell23_auto_01.png` | Arx → Pax4 | 与下一张互为镜像 |
| `inline_cell23_auto_02.png` | Pax4 → Arx | 同上 |
| `inline_cell24_auto_01.png` | Pax4 → Ins2 | β 细胞那条边 |

`jacobian_*` 系列（矫正链导出）是这些颜色的数字版：

| 文件 | 内容 |
|---|---|
| `jacobian_metadata.json` | 张量形状、17 个 regulator/effector 名单、轴序说明 |
| `jacobian_mean_all_cells.tsv` | 17×17，细胞平均；**行=regulator，列=effector**（表头已转置成人类习惯） |
| `jacobian_mean_ductal_cells.tsv` | 同上，仅 Ductal 916 细胞 |
| `jacobian_pairs.tsv` | 每条边一行：min/q25/median/mean/q75/max/`frac_positive` |
| `jacobian_pairs_by_cluster.tsv` | 边 × 群：均值/中位/正比例 |
| `jac_<reg>_<eff>.tsv` | 逐细胞值（cell, clusters, jacobian），可自己画分布或做检验 |

这几张表最值得先看的两条结论：

- **`Pdx1 → Ins2` 是全表最强的一条边**：在 Beta 群里均值 **+1.5e-3、99.3 % 细胞为正**，比其它任何边大 6 倍；
  而教材实际画的 `Pax4 → Ins2` 在 Beta 里是 **−4e-5（46 % 为正）**。教材 cell 24 的注释写的是
  「Pdx1 激活 Ins2」，代码画的却是 Pax4——**数据站在注释这边，代码画错了基因**（DYN-07）。
  若要用这张图讲课，建议改画 `regulators=["Pdx1"]`。
- **Pax4 ⤙ Arx 的互抑只发生在 Ngn3 high EP**：两个交叉项同时为负仅在该群
  （Arx→Pax4 −4e-5、40 % 为正；Pax4→Arx −1.2e-4、17 % 为正）；到 Pre-endocrine 已翻回正值，
  在 Alpha/Beta 两个交叉项都为正（互相激活）。也就是说教材说的「分叉点上的互相抑制」成立，
  但它是**局部**开关，不是全局开关。

**两个必须知道的实现事实**（否则数字会读反）：

1. 存进 `adata.uns['jacobian_pca']['jacobian']` 的是 **PCA 空间的 30×30×3696**。想要基因级结果，
   必须用 `jacobian_gene`（dynamo 用 `Qi·J·Qjᵀ` 反变换得到）。
2. `jacobian_gene` 的轴序是 **(effector, regulator, cell)** —— 本目录的导出按
   `Jg[effector_idx, regulator_idx]` 取值，表格再转置成"行=regulator"以便阅读。

`rank_acceleration_*.tsv`：`top200.tsv` 是教材 `adata.uns['rank_acceleration']` 的原貌（逐群宽表，
每个群一列、200 行），`top20_per_cluster.tsv` 是摘要。Ductal 前 10 是 Hmgb2/Hist1h2bc/Ccnb2/Cdc20/
Dynll1/Birc5/Cks2/Ccnb1/Cenpf/Hn1（9 个有丝分裂基因），Beta 前二是 Ins1/Ins2，Ngn3 high EP 是
Chga/Chgb/Tm4sf4/Akr1c19/Cryba2/Cpe/Hmgn3/Pcsk1n/Cystm1/Glud1。Beta 第三位是 `Xist`——这份数据混合性别，
那一项反映的是雌/雄差异。

## 6. 05_enrichment：GO 富集

`GO_BP_2021_Ductal_top100.tsv` 是 Enrichr 的完整返回：925 行 × 10 列
（`Gene_set, Term, Overlap, P-value, Adjusted P-value, Old P-value, Old Adjusted P-value, Odds Ratio,
Combined Score, Genes`）。**107 行**的 `Adjusted P-value < 0.05`；最小的校正 p 是 1.74e-21
（mitotic spindle organization，21/157）。

两张图：`inline_cell32_auto_01.png`（barplot）与 `inline_cell32_auto_02.png`（dotplot）。
它们的**标题是错的**：分别写着 `GO_Biological_Process_2018` 与 `KEGG_2016`，而查询的库是
`GO_Biological_Process_2021`，且从未查过 KEGG——标题是绘图调用里硬编码的字符串。
矫正版 `06_corrected/10_..._corrected.png` 把库名从 `res2d['Gene_set']` 读出，条形像素不变、标签正确。

`*.enrichr.reports.txt/pdf` 是 gseapy 自己倒出来的报告（教材用 `outdir='./figure'` 把它们丢进图目录，
这里收进 `05_enrichment/`）。带 `aswritten_` 前缀的是原样跑那份，不带的是矫正跑那份——两份内容一致，
留着是为了说明"为什么图目录里会长出 txt 和 pdf"。

## 7. 06_corrected：矫正链新增了什么

| 文件 | 修正项 | 内容 |
|---|---|---|
| `07_reg_Arx_Pax4_dyn_savefig.png` | C1 | 教材只内联显示的第二张 Jacobian 面板，按规矩存盘（1580×1048） |
| `08_reg_Pax4_Arx_dyn_savefig.png` | C1 | 第三张 |
| `09_reg_Pax4_Ins2_dyn_savefig.png` | C1 | 第四张 |
| `12_curl_nullcline_dyn_savefig.png` | C6 | 旋度图——教材算了 curl 却从没画过 |
| `10_Ductal_GO_enrichment_barplot_corrected.png` | C2 | 标题改成 `GO_Biological_Process_2021`，条形与原版像素一致 |
| `11_Ductal_GO_enrichment_dotplot_corrected.png` | C2 | 同上（dotplot） |
| `12_umap_shipped_vs_recomputed.png` | C4 | 左：文件自带 UMAP；右：`reduceDimension(enforce=True)` 的结果 |
| `umap_recompute_comparison.tsv` | C4 | 成对距离 Spearman **0.9393**、平均位移 9.74、展布 5.11、位移/展布 **1.905** |
| `ductal_subset_stats.tsv` | C7 | Ductal(916) vs 其余 的场量均值（让教材里死掉的 `progenitor` 变量派上用场） |

`umap_recompute_comparison.tsv` 值得单独说：它是"教材的降维到底算没算"这个问题的答案。
`reduceDimension` 因文件自带 `X_umap` 而**跳过**（控制台 `[UMAP] completed [0.0012s]`），
而强制重算会把每个细胞挪走 1.9 个展布单位——所以矫正跑**不**在原对象上改 `X_umap`，
只在副本上算、单独出一张对比图，从而保证图 01–09 逐字节不变。

## 8. 07_figure_comparison：与教师参考图对拍

`figure_comparison_aswritten.tsv` 一行一张图，列分三组：

- **对照 2026 脚本重跑图**（`ref2026_*`）：宽高、长宽比、`structural_r`（墨迹掩膜形状相关，
  对配色不敏感）、`ink_iou`、双方墨迹占比、配色数与 Top6 配色；
- **对照 2025 HTML 内嵌图**（`vshtml_*`）：同上一套指标；
- 元信息：该图出自教材哪个单元、用的是什么调用、图类型（stream/quiver/composite/graph/jacobian/bar/dot）。

两套参考**不能互换**：2025 是 notebook 内联显示的那 6 张，2026 是升级脚本重跑的 11 张。
结论速记：01/02 几乎完全一致（IoU 0.996 / 1.000），03–09 在 0.64–0.97 之间且都能解释
（dpi 差、等高线重涂、积分轨迹），**10/11 的 2026 版是纯白空画布**（墨迹 0、颜色 0、12–15 kB），
所以 GO 那两张只能拿 2025 内嵌图当真值（与本目录的 bar/dot 分别是 IoU 0.892 / 0.991）。

**练习**：打开表，按 `vs2026_ink_iou` 升序看前两行，再去 `03_vectorfield/` 与 `06_corrected/` 找对应图，
亲自确认"低 IoU ≠ 分析不同"。04 的 0.639 来自填充等高线（色阶挪一档就重涂大片），
05 的 0.792 来自积分轨迹的控制点落点。

## 9. 09_run_logs 与 scripts

`09_run_logs/aswritten_console.txt` 是原样跑的完整控制台（含每一步耗时）。想看证据时优先用它：

- `6956 genes passed basic filters` → §2 的过滤数；
- `applying PCA … <insert> X_pca to obsm` → PCA 是 recipe 算的（不是 no-op）；
- `[UMAP] completed [0.0012s]` → UMAP 是**跳过**的（DYN-03）；
- `incomplete neighbor graph info detected … indices not in adata.uns.neighbors` → 文件里没有 kNN 图，
  是 dynamo 自己建的；
- `[ddhodge completed] completed [21.0884s]` / `[VectorField] completed [24.1135s]` → 各步成本。

`scripts/` 六个脚本：

| 脚本 | 作用 |
|---|---|
| `extract_dynamo_cells.py` | 从 HTML 导出里抽出 50 个 cell（标准库 `html.parser`，环境无 bs4） |
| `dynamo_cells_aswritten.py` | 34 个代码单元按序还原并逐格执行（只替换教师硬编码路径） |
| `dynamo_corrected.py` | 同管线 + C1–C10 |
| `dynamo_vf_percell_export.py` | 导出逐细胞场量、逐群统计、检验、状态图矩阵 |
| `dynamo_jacobian_export.py` | 导出基因级 Jacobian（含 Pdx1→Ins2 等 10 条边）与张量元信息 |
| `compare_figs_dynamo.py` | 与两套参考图对拍，产出 `07_figure_comparison/*.tsv` |

## 10. 一句话记住这一讲

教材给出的四句定性判断里，**三句站得住**（Ngn3 high EP 是分叉点、早期祖细胞以高加速度出口、终末群散度为负），
**一句需要改写**（Ductal 不是散度意义上的鞍点，它是**最慢**且**潜势最低**的根）。
把颜色变成数字之后，这一讲最有价值的产物不是那五张图，而是
`03_vectorfield/vf_per_cell.tsv` —— 3696 行 × 8 列，任何后续统计离它都只有一个 `groupby`。
