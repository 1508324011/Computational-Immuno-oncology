# WES 分析结果 · Analysis Results

两轮实跑（demo 0.5× + 全量）的**分析产物摘录**——QC 报告、比对/去重指标、VQSR 模型与
tranches、变异集 VCF、体细胞检出、CNV 分段、HsMetrics 覆盖度、全部实跑脚本与作业日志。
**不含任何原始数据**（FASTQ / BAM / GVCF 未入库，仍在集群运行目录，见讲解文档 §5）。

- 📖 **先读讲解**：[WES_results_guide.html](WES_results_guide.html)（或
  [PDF](WES_results_guide.pdf) / [md](WES_results_guide.md)）——逐类文件怎么读、
  两轮对比、6 个动手练习（含亲手验证 VQSR 翻倍 bug 的三条命令）。
- `demo/`：0.5× 降采样轮（作业 3048），来源
  `/lustre1/user/bjx131_pkuhpc/wes_run_cluster_20260914`。
- `full/`：全量轮（作业 3049–3061，8 个依赖链作业），来源
  `/lustre1/user/bjx131_pkuhpc/wes_run_full_20260920`。

**目录编号即路线图步骤编号**（`0_fastq` → `5_somatic_tmb`），与
[矫正增强版路线图](../ref/WES_roadmap_corrected_zh.html)逐节对照。
