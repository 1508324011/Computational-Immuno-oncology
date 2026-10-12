# RNA-seq 下游应用分析结果讲解（results guide）

: 2026-10-12 · Lesson 3 (RNA-seq-application) 的实跑产物导读。配套矫正增强版见 `RNA_seq_application_corrected.{md,html,pdf}`（双语）。

## 1. 这批结果是什么

这是第三节课（RNA-seq 下游应用分析）在登录节点按教材原样实跑的全部产物：四个 Part（基因相关性、生存分析、免疫浸润、TIDE），共 31 个文件。**不含任何原始数据**——表达矩阵和随访表在共享只读路径 `/lustre1/share/rnaimmune/`，未复制进仓库；这里的都是分析产出（图、去卷积结果表、TIDE 输出、脚本、日志）。

## 2. 逐类导读

### 02_gene_correlation/ —— Part 1
- `GZMA_cor_PRF1.png`：GZMA vs PRF1 散点 + Pearson/Spearman 标注。两个杀伤效应基因强共表达（Pearson r=0.71，Spearman ρ=0.93）。
- `CTLA4_cor_CYT.png`：CTLA4 vs CYT 签名散点。注意 Pearson–Spearman 分裂（r=0.16 p=0.32 vs ρ=0.80 p≈0）——离群/非线性单调的教科书特征。

### 03_survival_analysis/ —— Part 2
- `survival_CYT.png`：KM 曲线，CYT high/low 各 21 例，带 log-rank p、置信区间、风险表。
- `cox.png`：多变量 Cox 森林图（Age/Gender/Stage/CYT/CYT_level/CTLA4），6 协变量均不显著（队列小+CYT 重复入模共线），Concordance=0.679。

### 04_immune_infiltration/ —— Part 3（重点）
- `qts_barplot.png`：quanTIseq 堆叠柱（**绝对分数法才能这样画**），前 10 样本，按细胞型着色。
- `qts_scatterplot.png`：quanTIseq 点/分面图（每细胞型自由 y 轴，对任何方法都安全）。
- `qts_heatmap.png` / `qts_heatmap_row_scaled.png`：quanTIseq 热图，未缩放 vs 按行 z-score。后者是公平跨样本比较。
- `deconv_methods/*.tsv`：七种去卷积方法的结果表（MCPcounter 为 NULL 占位+阻塞说明，其余六法完整）。`res_quantiseq.tsv` 是首次调用，`res_quantiseq_rmgenes.tsv` 是删 6 个 B 细胞基因后的覆盖版（下游图用此版）。
- `CIBERSORT-Results.txt`：CIBERSORT.R 源码向 cwd 写的原始输出（42 样本 × 22 型 + P/Correlation/RMSE 列）。
- `mcp_counter_blocked.txt`：APP-01 阻塞记录（conda MCPcounter 1.1.0 http 硬编码 → 登录节点 reset）。
- `log_part3.txt`：完整运行日志。

### 05_tide/ —— Part 4
- `res_tide.txt`：我们自己跑 `tidepy` 的输出（42 样本 × 14 列，字节尺寸与教师预跑副本一致）。
- `tide_heatmap.png`：TIDE 各指标热图（按行 z-score），列注释含真实/预测应答、CTL.flag、No.benefits、Age/Gender/Stage。
- `tide_predicted_responder.png` / `tide_responder.png`：TIDE 分柱图，按 CTL.flag 分面，分别按**预测**与**真实**应答着色——对比即读预测 vs 真实。
- `roc.png` / `roc_smooth.png`：ROC 曲线，AUC=0.7872（原始）/ 0.7718（平滑）。
- `log_tidepy.txt`：tidepy 运行日志（含三行预期 `[WARN]`）。

### */scripts/
每个 Part 一个脚本（`.R` 或 `.sh`），代码逐字取自教材，文件头登记所有偏差（交互式→Rscript、setwd、LC_ALL=C、跨 Part 对象重派生、MCPcounter tryCatch 等）。

## 3. 动手练习（全部仓库内可完成）

**练习 1 —— MCPcounter 阻塞复现**
```bash
cd /gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/04_immune_infiltration
cat mcp_counter_blocked.txt                       # 读阻塞记录
# 亲手验证 http vs https：
curl -sS -o /dev/null -w "http: %{http_code}\n" --max-time 15 \
  http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt
curl -sS -o /dev/null -w "https: %{http_code}\n" --max-time 15 \
  https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt
# 预期：http 000（reset），https 200
```

**练习 2 —— 六法维度核对**
```bash
cd /gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/04_immune_infiltration/deconv_methods
for f in res_*.tsv; do echo -n "$f: "; head -1 "$f" | awk -F'\t' '{print NF-1" 样本"}'; echo -n "  细胞型="; echo "$(($(wc -l < "$f")-1))"; done
# 预期：epic 8型、quantiseq 11型、xcell 39型、cibersort 22型、cibersort_abs 22型、timer 6型、mcp_counter 占位
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

第三节课是纯下游 R 分析 + 一行 tidepy，登录节点 ~6 分钟跑完，31 个产出文件全部入库；唯一阻塞是 conda MCPcounter 1.1.0 的 http 硬编码（已记录、未绕过），六法齐全，TIDE 输出与教师副本字节一致，ROC AUC 0.79 复现原文。
