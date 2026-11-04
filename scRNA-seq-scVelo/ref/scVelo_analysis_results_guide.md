# scVelo 分析结果讲解

本文件导读 `scRNA-seq-scVelo/results/` 下的全部产物，并给出可以亲手复现的练习。配套的矫正增强版文档在
`ref/scVelo_analysis_corrected.md`（英文）与 `ref/scVelo_analysis_corrected_zh.md`（中文）。

## 1. 怎么读这个目录

```
results/
├── 01_preprocess/          步骤 1：占比面板 + 每群细胞数 + 输入对象自身的审计表
├── 02_velocity/            步骤 2-3：stream/箭头面板、每细胞速率、速率图统计
├── 03_phase_portraits/     步骤 4：Cpe/Gnao1/Ins2/Adk 面板、Cpe 群面板、Hells/Top2a 面板
├── 04_gene_dynamics/       步骤 5：两张 top5 面板 + 排名表 + 每基因动力学
├── 05_cellcycle/           步骤 6：细胞周期渐变面板 + 相位 marker 面板 + 分数表
├── 06_corrected/           矫正链：13 张矫正面板 + 17 张表（dynamical 拟合、隐时间、
│                           参数敏感性、两模型排名并置、运行摘要）
├── 07_figure_comparison/   与课程随附 11 张参考 PNG 的机器可读对拍（原样版 / 矫正版）
├── 09_run_logs/            三次运行的控制台记录（原样、矫正、两次审计）
└── scripts/                重新生成以上全部内容的五个脚本
```

约定：

- `01_`–`05_` 目录里的图**原样跑**产出，文件名就是教材单元编号（`02_velocity_embedding_stream.png`
  对应第 63 单元）。它们与课程随附的 11 张参考图一一对应，可以直接 `diff` 眼睛。
- `06_corrected/` 里带 `corrected_` 前缀的是**矫正链**产出；不带前缀的 `sensitivity_*`、
  `rank_velocity_genes_*`、`latent_time_*` 是同一链条落下的表。
- 每个 `09_run_logs/*_console.txt` 里，`[describe]` 段打印对象形状与非零数，`[state]` 段打印该步骤
  结束后的关键量，`[VELO-nn]` 段是问题清单对应位置的现场取证。
- 表格一律 UTF-8 CSV，索引列优先；`proportions_by_cluster.csv` 有两个口径列，**不要**混用（见 §2）。

## 2. 01_preprocess：对象本身与占比

`01_proportions.png` 是教材第 50 单元那张饼图/柱图。柱高取自 `obs['initial_size_*']`，也就是过滤**之前**
的文库大小：Ductal 12.3%、Ngn3-low 12.4%、Ngn3-high 15.0%、Pre-endocrine 20.5%、Beta 20.7%、
Alpha 19.5%、Delta 19.3%、Epsilon 17.9%，整体 16.3%。

`proportions_by_cluster.csv` 把两套口径并排放：`unspliced_fraction_initial_size`（图上的）与
`unspliced_fraction_current_layers`（模型实际看到的）。后者整体 23.0%，逐群 19.9–25.3%。差 3–9 个
百分点，不是误差，是两个问题：前者说的是"文库里有大多比例是内含子"，后者说的是"速率模型手上有多少
unspliced 信号"。

`obs_cell_metadata.csv`（3,696 行）与 `var_gene_annotation.csv`（2,000 行）是**速率计算所用对象**的
逐细胞/逐基因元数据，包含 `velocity_gamma`（逐基因稳态斜率）、`velocity_r2`、`spearmans_score`、
`fit_*` 等。`cluster_cell_counts.csv` 是八群的规模（916 / 642 / 592 / 591 / 481 / 262 / 142 / 70）。

## 3. 02_velocity：场

`02_velocity_embedding_stream.png` 与 `03_velocity_embedding.png` 是全课的主图：前者是核平滑后的流线，
后者是逐细胞箭头。`11_velocity_embedding_stream_repeat.png` 与前者**逐字节相同**（同 md5）——教材第 84
单元重复了第 63 单元的调用，这是问题清单里的 VELO-10，保留它是为了让读者自己验证。

`velocity_per_cell.csv`（3,696 行）给每个细胞的速率模长 `velocity_norm`、所在群、以及
`velocity_confidence`；`velocity_summary_by_cluster.csv` 是逐群的模长中位数与置信度中位数。
`velocity_graph_aswritten.csv` 记录速率图的非零数（593,854）、密度（4.3%）与 k（30）——
图本身是 3,696 × 3,696 的稀疏矩阵，不入库。

`top_abs_velocity_genes_by_cluster.csv` 是按"速率绝对值均值"排的每群前若干基因，与教材的
`rank_velocity_genes`（按速率–表达相关排）不是同一个量，放在一起可以看出两种排法的差别。

## 4. 03_phase_portraits：相位图

三张面板，对应教材第 69、71、83 单元：`04_phaseportrait_markers.png`（Cpe/Gnao1/Ins2/Adk）、
`05_cpe_groups_velocity.png`（Cpe 按群+速率）、`10_phaseportrait_hells_top2a.png`（Hells/Top2a）。

`phase_portrait_gene_stats.csv` 逐基因记录相位图的关键量：`spliced_max`、`unspliced_max`、细胞点云的
相关系数、以及稳态斜率 `velocity_gamma`。读法：`Ins2` 的 spliced 上到 120、点云近乎一条直线，说明成熟
转录本主导、相位图信息最少；`Cpe` 有干净的环，是教材选它讲"速率方向"的原因。

## 5. 04_gene_dynamics：速率基因与每基因动力学

`06_top5_ngn3_high_ep.png` 与 `07_top5_pre_endocrine.png` 是教材第 75 单元的两张 top5 面板（各 5 个
基因的 UMAP+速率）。`rank_velocity_genes_top50.csv` 是每群前 50 名（列=群，行=名次），
`rank_velocity_genes_top5_summary.csv` 是每群前五名的拼接字符串。

两个必须在正文里说清的事实：`min_corr=.3` 筛的是 `var['spearmans_score']` 而**不是**排名相关性，而
真正的截断是 `n_genes=100`（默认值）。实测：阈值 0.1/0.3/0.5/0.7 分别留下 1,861/1,469/973/430 个
合格基因，而 Ductal 与 Beta 的可见 top-3 一模一样——阈值在 0.9 之前从未生效。

`gene_kinetics_scores.csv`（2,000 行）是逐基因的 `velocity_gamma`（中位 0.343、IQR [0.090, 0.966]）、
`velocity_r2`（中位 0.052）、`velocity_qreg_ratio`、`spearmans_score`。中位 r² 只有 5%，是所有
"速率基因"结论该被如何看待的前提。

## 6. 05_cellcycle：细胞周期

`08_cellcycle_gradients.png` 是教材第 79 单元的 S/G2M 渐变图（**唯一**与参考图对不上的那张，结构相关
0.589），`09_phase_marker_genes.png` 是第 81 单元的相位 marker 面板。

`cellcycle_scores_per_cell.csv` 是**重算**后的分数（教材的 `score_genes_cell_cycle` 覆盖了数据自带的
那一对）；`cellcycle_summary_by_cluster.csv` 是逐群的相位计数（G1 2,537 / S 878 / G2M 281）与均值；
`phase_marker_gene_ranking.csv` 与 `phase_marker_top5.csv` 是逐群的 S/G2M marker 排名
（S：Mcm6、Hells、Tyms、Pola1、Uhrf1；G2M：Hn1、Top2a、Kif23、Ckap2l、Kif2c）。

## 7. 06_corrected：矫正链做了什么

13 张图里，`corrected_01`–`corrected_03`、`06`–`08`、`10` 与教材同位置的面板同源（用来对拍），
`corrected_04`/`05` 是相位图的两种画法（稳态 / 只画 dynamical 拟合），`corrected_09` 用数据自带分数
重画了细胞周期，`corrected_11` 用 dynamical 拟合重画了 Hells/Top2a，`corrected_12`/`13` 是教材从来
没画过的两张（隐时间+置信度、以隐时间上色的 stream）。

关键表：

| 文件 | 内容 | 值得注意的数字 |
|---|---|---|
| `run_summary_corrected.csv` | 矫正对象自检 | 3,696 × 2,000，`X` 非零 1,740,903，`Ms` 非零合计 5,778,956.5 |
| `latent_time_summary_by_cluster.csv` | 逐群隐时间 | Ductal 0.157 → Ngn3-high 0.465 → Beta 0.774 |
| `latent_time_and_confidence_per_cell.csv` | 每细胞隐时间+两版置信度 | 两模型置信度相关仅 0.207 |
| `rank_velocity_genes_top5_dynamical.csv` | dynamical 排名 | 与 stochastic 的 top-20 Jaccard 0.132 |
| `sensitivity_graph_and_hvg.csv` | 六条完整流水线 | 1000 HVG → Jaccard 0.50；k=15 → 0.58；重算 PCA → 0.78 |
| `sensitivity_min_corr.csv` | min_corr 扫描 | 0.1→1,861 个合格、0.9→27 个 |
| `sensitivity_n_genes_cap.csv` | n_genes 上限 | 10/50/100/300 就是表格行数 |
| `cellcycle_scores_original_vs_recomputed.csv` | 两版周期分数 | Pearson 0.964 / 0.973 |
| `proportions_by_cluster_corrected.csv` | 两套占比口径 | 16.3% vs 23.0% |

`corrected_summary.json` 是整条链条的机器可读摘要（版本、耗时 4.72 min、各项统计）。

## 8. 07_figure_comparison：和教师的 11 张图对拍

`figure_comparison_aswritten.tsv` 与 `figure_comparison_corrected.tsv` 各 11/13 行，列包括画布尺寸、
长宽比、结构相关 `struct_corr`、墨迹 IoU `ink_iou`、色彩分布。

原样跑：十张 0.972–0.993，第八张（细胞周期）0.589。矫正链：十张可比面板 0.904–0.986，细胞周期仍是
0.579。差异的原因不在数据：那张参考图与我们的画布（1205 × 648）和图例文字完全一致，但标记云内部配色
不同，且参考图集里存在三种 dpi（80 / 96 / 64），说明那 11 张是多次会话拼接的。

## 9. 练习

**练习 1（10 分钟，口径）。** 打开 `proportions_by_cluster.csv`，用 `unspliced_sum/(spliced_sum+unspliced_sum)`
自己算一遍逐群占比，然后与 `01_proportions.png` 上的柱高比较。差多少？哪一列对应图上的数？再看
`scv.pl.proportions` 的签名里 `use_raw` 默认值。

**练习 2（15 分钟，模型）。** 用 `corrected_adata.h5ad`（或重跑脚本）分别取
`adata.var['velocity_gamma']` 与 `adata.var['fit_gamma']`，画散点并算 Spearman。为什么两者秩相关是
负的（−0.45）？提示：稳态模型把 β 归一为 1，dynamical 学的是独立 β。

**练习 3（10 分钟，隐时间）。** 用 `latent_time_and_confidence_per_cell.csv` 的 `latent_time` 列按群画箱线图。Beta 的中位数比 Alpha 早
（0.742 vs 0.758），但均值更晚（0.774 vs 0.758）。画出两者的分布，解释为什么"全局隐时间排出单一顺序"
在这里会失真。

**练习 4（20 分钟，参数）。** 把 `scVelo_audit_run.py` 的 `run_variant` 复制出来，把 `n_top_genes` 设成
500，跑一遍，比较 top-20 Jaccard 与速率基因数。趋势与 `sensitivity_graph_and_hvg.csv` 里的 1,000 那一行
一致吗？

**练习 5（15 分钟，图的含义）。** 用 `corrected_04`/`corrected_05` 两张相位图对照：哪些基因的 dynamical
拟合曲线明显偏离稳态线？结合 `phase_portrait_gene_stats_corrected.csv` 里的 `velocity_gamma`，说说偏离
大的基因有什么共同点。

**练习 6（10 分钟，细胞周期）。** 用 `cellcycle_scores_original_vs_recomputed.csv` 算两版的绝对差，
找出差最大的 20 个细胞，看它们在 `clusters` 里属于哪群。为什么差异集中在某些群？

**练习 7（20 分钟，对拍）。** 打开 `figure_comparison_aswritten.tsv`，找出 `ref_px` 与 `ours_px` 比值
正好 0.8 与 1.2 的行。把 `scv.settings.set_figure_params(dpi=64)` 打开，重画其中任意一张，验证画布
变成与参考一致——然后解释为什么"画布一致"并不等于"分析一致"。

**练习 8（30 分钟，往下走）。** 用 `corrected_adata.h5ad` 跑 `ref/scVelo_analysis_corrected.md` 附录 B 的
三行 CellRank，得到每个细胞的命运概率，并与隐时间对照：内分泌终末群的命运概率是否与它们的隐时间排序
一致？哪一群最不一致？
