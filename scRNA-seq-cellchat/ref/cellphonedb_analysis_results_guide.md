# CellphoneDB 分析结果导读

本文件导读 `scRNA-seq-cellchat/results/` 下的全部产物，并给出可以亲手复现的练习。
配套的矫正增强版文档在 `ref/scRNAseq_cellphonedb_corrected.md`（英文）与
`ref/scRNAseq_cellphonedb_corrected_zh.md`（中文）；教材原件的 cell 抽取结果在
`results/scripts/cellchat_cells_aswritten.py`。

## 1. 怎么读这个目录

```
results/
├── 01_input/               数据侧事实：细胞计数、微环境、9x9 矩阵、每对/每互作汇总、
│                           四张统计表的行序审计
├── 02_simple/              cpdb_analysis_method 的 4 张表（不算统计）
├── 03_statistical/         cpdb_statistical_analysis_method 的 8 张表（1000 次置换）
├── 04_figures/             教程自己画出来的 3 张图（热图 + 2 张点图）
├── 05_corrected/           矫正链：24 列 vs 36 列面板、x 轴对数、检索结果、阈值扫描表
├── 06_figure_comparison/   与两套教师参考图的机器可读对拍表
├── 07_degs/                DEG 分支：8 张表 + 2 张点图 + 输入 DEG 表
├── 09_run_logs/            每条执行线的完整控制台、确定性核对、session info
└── scripts/                重新生成以上全部内容的脚本
```

约定：

- `02_`–`04_` 是**原样跑**（`cellchat_cells_aswritten.py`）的产物，文件名沿用 CellphoneDB 写盘时的
  `statistical_analysis_*.txt` 原名，只去掉了时间戳后缀，因此 `03_statistical/means.tsv` 就是那次调用的
  `statistical_analysis_means_<时间戳>.txt`。
- `05_`–`07_` 是矫正/取证/DEG 三线的产物，**不与原样跑混放**，因为它们的口径不同（C1 改了拼写、
  C2 换了线程、C10 换了判据）。
- 一切表都是 **tab 分隔**，带表头，第一列是 `id_cp_interaction`（`deconvoluted*` 是 `id_cp_interaction` +
  `gene_name`），可直接 `pd.read_csv(..., sep='\t')`。
- 教程给的"417 个显著互作"在 `03_statistical/significant_means.tsv` 里体现为 4248 个非 NaN 格；
  数"行"会得到 2440，数"列"会得到 95——**都对，但都不是 417**。

## 2. 01_input：数据侧的全部事实

| 文件 | 内容 | 怎么用 |
|---|---|---|
| `cell_type_counts.tsv` | 42 个标注的细胞数（`dS2` 588 … `eEVT` 28） | 看清"矩阵里 42 个类型"这一事实 |
| `microenvironments.tsv` | 9 个类型 → 微环境的对应 | 与上一行对照，就是 §2.6 的 33 个被丢弃类型 |
| `sig_matrix_9x9.tsv` | 81 个有序对的显著互作数（行和 411–565，列和 230–783） | 本文 §4.1 的全部内容 |
| `sig_per_cell_pair.tsv` | 每对的显著互作数（降序） | 排序与矩阵行和一致，用于核对 |
| `sig_top_interactions.tsv` | 每互作的显著对数（前 200） | §4.2 的榜单，**必须按 id 合并生成**（见下） |
| `table_row_order.tsv` | 四张统计表"同行是否同 id"的比例 | CPB-19 的证据：`significant_means` 只有 0.0004 |

`table_row_order.tsv` 是这一讲最该先看的文件。它说明 `pvalues` / `means` / `interaction_scores` 三张表
行序一致（比例 1.000），而 `significant_means` 与它们**完全不同**（0.0004）。于是

```python
# 错：按位置对齐，结果看起来像样
n = sm[PAIR_COLS].notna().sum(axis=1)
top = pd.DataFrame({'interaction': means['interacting_pair'], 'n': n})   # 名字与计数错位

# 对：按 id 合并
top = sm[['id_cp_interaction']].assign(n=n).merge(
        means[['id_cp_interaction', 'interacting_pair']], on='id_cp_interaction')
```

我的第一版 Top 表就踩了这个坑：把 `COL4A1_integrin_a1b1_complex` 的 76 挂到了
`IAPP_AMY3_receptor` 头上（后者真实值是 0）。`scripts/evidence_recheck_order.py` 是修正版，
`09_run_logs/evidence_log.txt` 里有当时的运行记录。

## 3. 02_simple：不算统计的那四张表

| 文件 | 形状 | 读法 |
|---|---|---|
| `means_result.tsv` | 2440 × 94 | 13 列元信息 + 81 对；**NaN = 没过 10 % 掩码** |
| `deconvoluted.tsv` | 6825 × 49 | 长表：一（互作，基因）一行，复合物占多行 |
| `deconvoluted_percents.tsv` | 6825 × 49 | 该基因在该簇的表达细胞比例 |
| `interaction_scores.tsv` | 2440 × 94 | 参与者表达乘积打分，点图的透明度来源 |

练习：在 `deconvoluted_percents` 里查 `IAPP`、`AMY3`、`PTPRS` 三行，你会看到
`IAPP` 只在 eEVT 有 7 % 细胞表达、其余簇 0 %；而 `PTPRS` 在 iEVT 54 %、GC 58 %、eEVT 57 %。
这解释了为什么"胰腺激素"类互作在这份蜕膜数据里能"显著"（0 对）或不能——**掩码是逐基因的**。

## 4. 03_statistical：1000 次置换之后的 8 张表

| 文件 | 形状 | 关键事实 |
|---|---|---|
| `means.tsv` | 2440 × 94 | 与 `02_simple/means_result.tsv` 同内容 |
| `pvalues.tsv` | 2440 × 94 | **96 % 的格是 1.0**（= 没测）；最小非零 p 是 0.004 |
| `significant_means.tsv` | 2440 × 95 | 4248 个非 NaN 格；行序与其他表不同（CPB-19） |
| `deconvoluted.tsv` / `_percents.tsv` | 6825 × 49 | 长表重算 |
| `interaction_scores.tsv` | 2440 × 94 | |
| `CellSign_active_interactions.tsv` | 6 × 95 | 6 个互作：BTC_EGFR、LGALS3_MERTK、DLK1/JAG1/JAG2_NOTCH3、EPO_EPOR |
| `CellSign_active_interactions_deconvoluted.tsv` | 28 × 9 | 长表形式 |

练习（三分钟，能改写你对这一讲的理解）：

```python
import pandas as pd
pv = pd.read_csv('03_statistical/pvalues.tsv', sep='\t')
P  = [c for c in pv.columns if '|' in c]
v  = pv[P].values.ravel()
print((v == 1).mean())        # 0.96 -> 96% 没测
print((v == 0).sum())         # 2946 -> 被截断在 1/1000
print((v <= 0.05).sum())      # 4248 -> 教程的"显著"
```

## 5. 04_figures：教程画出来的三张图

- `inline_cell30_auto_01.png` — `plot_cpdb_heatmap` 的 "Sum of significant interactions"，
  实为 81 列 p 值表在 0.05 下的**行和**（9 个数），丢掉了方向。
- `inline_cell31_gg.png` — 两基因点图，**24 列**：`"PV MMPP11"` 拼错被静默丢弃（CPB-04）。
- `inline_cell32_gg.png` — 同一张点图按 `classification` 分面。

三张图对 2026 教师重跑的对拍是 r = 0.999 / 0.997 / 0.999，墨迹 IoU 0.987 / 0.899 / 0.937（见第 7 节）。

## 6. 05_corrected：矫正链新增了什么

| 文件 | 是什么 | 一句话 |
|---|---|---|
| `31_PV_trophoblast_aswritten_typo.png` | 原样（拼错）面板 | **24 列** |
| `31_PV_trophoblast_corrected.png` | C1 矫正面板 | **36 列**，第三个 PV 类型回来了 |
| `31_xaxis_pairs_*.tsv` | 两张面板的 x 轴对数 | 24 vs 36，可机器核 |
| `search_results.tsv` | 第 28 格检索的 58 行 | 全是 TGF-β；`CSF1_CSF1R` 被分类过滤删掉（CPB-12） |
| `threshold_sensitivity.tsv` | 0.05 / 0.10 / 0.20 三档 | 7007/708、4251/417、2413/218 |
| `table_stat_significant_means.tsv` | 单线程参考运行的显著表 | md5 `fcf4c9dc…`，与第二次单线程跑逐字节相同 |

注意 `search_results.tsv` 的行数（58）**不是**"检索到 58 个显著互作"：它是 5 个互作 × 若干细胞对的
长表。要数互作用有多少个，得 `nunique`。

## 7. 06_figure_comparison：与教师参考图对拍

`results_figcmp.tsv` 每行一张图，两套参考集各一列组：

| 图 | vs 2026 重跑 (r / IoU) | vs 2025 HTML (r / IoU) |
|---|---|---|
| 第 30 格热图 | 0.999 / 0.987 | 0.999 / 0.979 |
| 第 31 格点图 | 0.997 / 0.899 | 0.882 / 0.521 |
| 第 32 格分面图 | 0.999 / 0.937 | 0.861 / 0.462 |
| C1 矫正面板 | 0.977 / 0.835 | 0.880 / 0.545 |

2025 HTML 一列的 IoU 偏低不是"图画错了"，而是**渲染不同**：内嵌图的包围盒与 dpi 与我们不同
（2200×800 vs 2160×760；2600×1800 vs 2552×1751）。表里另有 `*_wh`、`*_aspect_ratio` 两列可自行核对。
把 IoU 当"相似度"读会误判；它是"渲染相似度"。

## 8. 07_degs：教程关掉的那一支

`degs_summary.tsv` 是总表，`degs_A_asis_*.tsv` 是 8 张标准表，`A_asis/degs_used.tsv` 是实际喂进去的
`cluster / gene` 两列表。

| 数字 | 值 |
|---|---|
| 附带 DEG 表 | 13 809 行、**6 个簇**、100 % 的行 `adj.P.Val < 0.05`（已是显著清单，再过滤是空操作） |
| 覆盖的微环境类型 | **6 / 9**——三个 PV 类型一个都没有 |
| DEG 分支结果 | 4191 格 / 419 互作 / 81 对 |
| 与置换法的重合 | 格 Jaccard **0.469**、互作 Jaccard **0.742** |

反直觉但重要：DEG 分支**不是**更严格的过滤器。`permutation only` 有 1552 格，`DEG only` 有 1495 格，
两者差不多大——它们问的是不同问题。

## 9. 09_run_logs 与 scripts

- `aswritten_console.txt`：原样跑全文。里面能看到 1000 次置换的进度条（26 次/秒，41 s）、
  `Time taken` 之外的每个阶段起止时间戳。
- `corrected_console.txt` / `corrected_summary.json`：矫正跑的控制台与机器可读汇总
  （阈值扫描、`active_TFs` 有无、100 次迭代、DEG 失败原因）。
- `determinism.json` / `determinism_console.txt` / `determinism.tsv`：单线程两次 + 5 线程一次的 md5 对拍。
  读法：`a`、`b` 两行 md5 相同 → 单线程可复现；`c` 行不同 → `threads=5` 不可复现。
- `evidence_log.txt`：矩阵、行序审计、p 值结构、检索复算的完整输出。
- `degs_console.txt`：DEG 分支；第一版把附带 limma 表直接喂进去抛的
  `AttributeError: id_cp_interaction` 原文就在里面。
- `session_info.txt`：cellphonedb 5.0.1 / ktplotspy 0.3.3 / plotnine 0.14.5 / pandas 2.2.0 等。
- `patches_applied.txt`：原样跑施加的三处偏离（P1、P2、H1）及其理由。

`scripts/` 与本文各节一一对应：`cellchat_cells_aswritten.py`（原样）、`cellchat_corrected.py`（矫正）、
`cellchat_evidence.py`（取证）、`cellchat_degs.py`（DEG）、`cellchat_plots.py`（五张分析图）、
`compare_figs_cellphonedb.py`（对拍），以及生成原样脚本的 `make_driver.py`、`extract_cellchat_cells.py`
和修正行序 bug 的 `evidence_recheck_order.py`。

## 10. 一句话记住这一讲

**4248 个"显著互作"里，96 % 的矩阵根本没被检验；被检验的部分没有做多重校正；把阈值从 0.1 改到 0.05
或 0.2，这个数字就在 7007 和 2413 之间；把线程从 5 改成 1，它又是另一个数。这一讲真正的价值不是
"417 个互作"，而是"如何把一个通讯推断结果读到可以写进论文的程度"。**

练习题（都能在本目录内完成，不需要重跑）：

1. 用 `03_statistical/pvalues.tsv` 对 81 列做 BH 校正，报告 q ≤ 0.05 的格数；再与 4248 对比。
   （答案在 `05_corrected/bh_fdr.tsv`：**2946 格 / 311 互作**——而且正好是 p = 0 的那些格。）
2. 用 `01_input/sig_matrix_9x9.tsv` 算"发出−接收"的对角化程度，找出最不对称的一对细胞类型。
3. 用 `07_degs/degs_A_asis_significant_means.tsv` 与 `03_statistical/significant_means.tsv` 按
   `id_cp_interaction`（**不是**按行位置）求交集，验证 2696 与 356 两个数。
4. 在 `05_corrected/search_results.tsv` 里加回 `CSF1_CSF1R`：把 `query_classifications` 去掉重跑
   `search_utils`（30 秒），看 6 行长表长什么样。
