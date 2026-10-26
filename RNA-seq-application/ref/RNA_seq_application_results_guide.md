# RNA-seq 下游应用分析结果讲解（results guide）

: 2026-10-12 · Lesson 3 (RNA-seq-application) 的实跑产物导读。配套矫正增强版见 `RNA_seq_application_corrected.{md,html,pdf}`（双语）。

## 1. 这批结果是什么

这是第三节课（RNA-seq 下游应用分析）在登录节点按教材原样实跑的全部产物：四个 Part（基因相关性、生存分析、免疫浸润、TIDE），共 31 个文件。**不含任何原始数据**——表达矩阵和随访表在共享只读路径 `/lustre1/share/rnaimmune/`，未复制进仓库；这里的都是分析产出（图、去卷积结果表、TIDE 输出、脚本、日志）。

## 2. 逐类导读

### 02_gene_correlation/ —— Part 1
- `GZMA_cor_PRF1.png`：GZMA vs PRF1 散点 + Pearson/Spearman 标注。两个杀伤效应基因强共表达（Pearson r=0.71，Spearman ρ=0.93）。
- `CTLA4_cor_CYT.png`：CTLA4 vs CYT 签名散点。注意 Pearson–Spearman 分裂（r=0.16 p=0.32 vs ρ=0.80 p≈0）——离群/非线性单调的教科书特征。

![图 1 — `CTLA4_cor_CYT.png`：CTLA4 对 CYT 签名。](./figures/CTLA4_cor_CYT.png){width=70%}

把这张图当作"为什么两个系数都要报"的教具：云在秩上几乎完全单调（ρ=0.80），但幅度上被少数高 CTLA4 样本压歪（r=0.16）。只报 Pearson 会得出"没关系"，只报 Spearman 会高估线性可分性——两个一起报，再指着散点说清楚是哪一种关系。

### 03_survival_analysis/ —— Part 2
- `survival_CYT.png`：KM 曲线，CYT high/low 各 21 例，带 log-rank p（**p = 0.027，显著**）、置信区间、风险表。
- `cox.png`：多变量 Cox 森林图（Age/Gender/Stage/CYT/CYT_level/CTLA4），6 协变量均不显著（队列小+CYT 重复入模共线），Concordance=0.679。

![图 2 — `survival_CYT.png`：KM 曲线（CYT high/low 各 21 例，带风险表）。](./figures/survival_CYT.png){width=70%}

先看两条数字：log-rank p = 0.027（刚好显著），而风险人数在一两年内就减半、之后每条曲线只剩个位数病人——所以"显著"是 42 个病人给出的显著，不是稳定估计。这张图的用途是展示流程与方向，不是证明疗效。

![图 3 — `cox.png`：多变量 Cox 森林图。](./figures/cox.png){width=70%}

所有 CI 都跨过 1，这正是"29 个事件上放 6 个协变量"的结果（约 5 事件/变量，低于 10 EPV）。变量清单本身就暴露了结构问题：`CYT` 与 `CYT_level` 是同一变量的两份，互相抬高标准误。要引用就引 Concordance 0.679。

### 04_immune_infiltration/ —— Part 3（重点）
- `qts_barplot.png`：quanTIseq 堆叠柱（**绝对分数法才能这样画**），前 10 样本，按细胞型着色。
- `qts_scatterplot.png`：quanTIseq 点/分面图（每细胞型自由 y 轴，对任何方法都安全）。
- `qts_heatmap.png` / `qts_heatmap_row_scaled.png`：quanTIseq 热图，未缩放 vs 按行 z-score。后者是公平跨样本比较。
- `deconv_methods/*.tsv`：七种去卷积方法的结果表（含 MCPcounter 矫正版 11 细胞型 × 42 样本）。`res_quantiseq.tsv` 是首次调用，`res_quantiseq_rmgenes.tsv` 是删 6 个 B 细胞基因后的覆盖版（下游图用此版）。
- `CIBERSORT-Results.txt`：CIBERSORT.R 源码向 cwd 写的原始输出（42 样本 × 22 型 + P/Correlation/RMSE 列）。
- `mcp_counter_blocked.txt`：APP-01 阻塞与矫正记录（conda MCPcounter 1.1.0 http 硬编码 → 登录节点 reset；随后用本地 https 标志集矫正跑出，11 细胞型 × 42 样本）。
- `log_part3.txt`：完整运行日志。

![图 4 — `qts_barplot.png`：quanTIseq 的绝对比例，堆叠柱画法在这里是正确的。](./figures/qts_barplot.png){width=70%}

读这张图的第一件事不是看免疫细胞，而是看浅色的"未分类细胞"有多大——它占了每根柱子的六到九成，免疫细胞只是右端那些细条。这是方法的正确输出（肿瘤里确实有大量细胞无 signature 可解释），但也是"堆叠柱不适合看低丰度群体"的直接演示。

![图 5 — `qts_scatterplot.png`：点/分面图，每个细胞类型独立 y 轴。](./figures/qts_scatterplot.png){width=70%}

任何方法都安全的选择。当某个群体永远不超过 5% 时，堆叠柱里它是看不见的，而这里它有自己的一格。要比"同一种细胞在不同样本间"的高低，读这张。

![图 6 — `qts_heatmap_row_scaled.png`：按细胞类型 z-score 后的热图。](./figures/qts_heatmap_row_scaled.png){width=70%}

公平的跨样本比较，本节最该看的一张。按行缩放后每种细胞类型各有一条自己的色阶，于是 B 细胞这种数值极小的群体也能看出"谁高谁低"。两个提醒：z 分数是相对的（红 = 相对本队列偏高，不是绝对高）；42 个样本的行缩放会让单点离群看起来夸张，任何惊人的行都回到图 5 核对。

### 05_tide/ —— Part 4
- `res_tide.txt`：我们自己跑 `tidepy` 的输出（42 样本 × 14 列，字节尺寸与教师预跑副本一致）。
- `tide_heatmap.png`：TIDE 各指标热图（按行 z-score），列注释含真实/预测应答、CTL.flag、No.benefits、Age/Gender/Stage。
- `tide_predicted_responder.png` / `tide_responder.png`：TIDE 分柱图，按 CTL.flag 分面，分别按**预测**与**真实**应答着色——对比即读预测 vs 真实。
- `roc.png` / `roc_smooth.png`：ROC 曲线，AUC=0.7872（原始）/ 0.7718（平滑）。
- `log_tidepy.txt`：tidepy 运行日志（含三行预期 `[WARN]`）。

![图 7 — `tide_heatmap.png`：TIDE 各指标按行 z-score 的热图，列注释含真实/预测应答与 CTL.flag。](./figures/tide_heatmap.png){width=70%}

对着注释条竖着读：CTL/Dysfunction/Exclusion 这块就是 TIDE 的两轴逻辑。注意它把真实应答（`Responder`）也画在同一条注释带上，所以这张图可以直接回答"预测和真实在哪些病人上打架"。

![图 8 — `tide_predicted_responder.png` 与 `tide_responder.png`：同一批柱图，分别按预测与真实着色。](./figures/tide_predicted_responder.png){width=70%}

并排看就是读"预测 vs 真实"最快的办法：多数柱子两边一致，少数不一致——那些不一致就是 ROC 概括掉的错误。TIDE 的规则是 `TIDE < 0 → 预测响应`，所以贴近 0 的病人判定会随 `--vthres` 翻转，值得单独列出来。

![图 9 — `roc.png` / `roc_smooth.png`：AUC 0.7872（原始）与 0.7718（平滑）。](./figures/roc.png){width=70%}

阶梯状不是模型的性质，而是 42 个病人只有 42 个阈值。能读的只有它相对虚线的位置（AUC 0.787 ≈ 随机抽一对时 79% 排序正确），以及平滑后掉的那 0.015——原始 AUC 偏乐观，因为它是在同一批 42 个点上拟合出来的。报数时用平滑值 + `ci.auc()`。

### */scripts/
每个 Part 一个脚本（`.R` 或 `.sh`），代码逐字取自教材，文件头登记所有偏差（交互式→Rscript、setwd、LC_ALL=C、跨 Part 对象重派生、MCPcounter tryCatch 等）。

## 3. 动手练习（全部仓库内可完成）

**练习 1 —— MCPcounter 阻塞与矫正复现**
```bash
# 1a. 读阻塞与矫正记录
cat /gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/04_immune_infiltration/mcp_counter_blocked.txt
# 1b. 亲手验证 http vs https（原样跑为何阻塞）：
curl -sS -o /dev/null -w "http: %{http_code}\n" --max-time 15 \
  http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt
curl -sS -o /dev/null -w "https: %{http_code}\n" --max-time 15 \
  https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt
# 预期：http 000（reset），https 200
# 1c. 矫正脚本用本地标志集跑（离线）：
source /lustre1/share/miniconda3/etc/profile.d/conda.sh && conda activate rnaimmune
export LC_ALL=C LANG=C
Rscript /gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/04_immune_infiltration/scripts/part3_mcp_counter_fix.R
# 预期：11 cell types x 42 samples，写入 deconv_methods/res_mcp_counter.tsv
```

**练习 2 —— 七法维度核对 + MCPcounter 重复映射验证**
```bash
cd /gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/04_immune_infiltration/deconv_methods
for f in res_*.tsv; do echo -n "$f: "; head -1 "$f" | awk -F'\t' '{print NF-1" 样本"}'; echo -n "  细胞型="; echo "$(($(wc -l < "$f")-1))"; done
# 预期：epic 8型、quantiseq 11型、xcell 39型、cibersort 22型、cibersort_abs 22型、timer 6型、mcp_counter 11型
# MCPcounter 重复映射验证（Monocyte 与 Macrophage/Monocyte 应字段完全相同）：
m=$(awk -F'\t' '$1=="Monocyte"{print substr($0,length($1)+2)}' res_mcp_counter.tsv)
ma=$(awk -F'\t' '$1=="Macrophage/Monocyte"{print substr($0,length($1)+2)}' res_mcp_counter.tsv)
[ "$m" == "$ma" ] && echo "Monocyte ≡ Macrophage/Monocyte（immunedeconv cell_type_map 把 Monocytic lineage 映射成两个名）" || echo "不同"
```

**练习 3 —— TIDE 输出与教师副本字节比对**
```bash
wc -c /gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/05_tide/res_tide.txt \
      /lustre1/share/rnaimmune/res_tide.txt
# 预期：两者都 10259 字节
diff <(head -1 .../05_tide/res_tide.txt) <(head -1 /lustre1/share/rnaimmune/res_tide.txt)  # 表头应一致
```

**练习 4 —— Pearson–Spearman 分裂亲眼确认**
打开 `02_gene_correlation/CTLA4_cor_CYT.png`，对照矫正版 §3 的数字表：Pearson r=0.16（p=0.3247，不显著）与 Spearman ρ=0.80（p≈0）并存——散点上能看到秩关系强但线性关系被离群点拉偏。

**练习 5 —— TIDE 分数与应答方向**
```bash
# TIDE<0 预测应答：看预测应答者的 TIDE 分符号
cd /gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/05_tide
awk -F'\t' 'NR>1{print $1, $3, $2}' res_tide.txt | sort -k2 -n | head -10
# 前几行（TIDE 最负）应多为 Responder=True
```

**练习 6 —— Cox 共线性亲验**
重跑 Part 2 脚本但去掉 `CYT_level`（公式改 `~ Age + Gender + Stage + CYT + CTLA4`），对比 `summary(cox_model)` 的标准误——CYT 的 se 应缩小（因不再与自身的二值化共线）。

**练习 7 —— 绝对分数 vs 相对口径可视化对照**
打开 `qts_barplot.png`（堆叠柱，和为 1，仅对 quanTIseq 这类绝对法正确）和 `qts_scatterplot.png`（点图，自由 y 轴，任何方法都安全）。对照矫正版 §5 的"输出口径"表，理解为何 MCP-counter/xCell/TIMER 不能画堆叠柱。

## 4. 未入库的大文件在哪

本模块**无大文件**——所有产出都是小图、小表、小日志，已全部入库。原始数据（表达矩阵、随访表）在共享只读路径 `/lustre1/share/rnaimmune/`，未复制（也无需复制，教材命令直接读该路径）。

## 5. 一句话总结

第三节课是纯下游 R 分析 + 一行 tidepy，登录节点 ~6 分钟跑完；MCPcounter 原样跑受阻（conda http 硬编码），已用本地 https 标志集矫正跑出，七法齐全；发现 MCPcounter 的 Monocyte 与 Macrophage/Monocyte 是 immunedeconv 映射表重复（10 独立群非 11）；TIDE 输出与教师副本字节一致，ROC AUC 0.79 复现原文。
