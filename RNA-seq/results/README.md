# RNA-seq/results — 两轮实跑分析结果（curated）

本目录是 RNA-seq 教程（mini 轮 2026-09-21 + homework/全量轮 2026-10-09）的分析产物入库，
与 `WES/results/` 同一策略：**只放我们自己的分析结果，不含任何原始数据**（无 FASTQ/BAM/索引；
`.gitignore` 第 9 节开了显式例外区并保留原始数据兜底封锁，pre-commit 钩子同步放行）。

逐文件导读与动手练习见 `../ref/RNA_seq_results_guide.md`（HTML/PDF 同目录）；
每一步的原理、命令原文、实测数字与坑位取证见 `../ref/RNA_seq_tutorial_corrected.md`（英文版）
与 `../ref/RNA_seq_tutorial_corrected_zh.md`（中文版）。

```
RNA-seq/results/
├── mini/                    # mini 轮（2026-09-21，10k reads/样本，chr22 索引，~4 分钟）
│   ├── 02_fastqc/           # 12 个 fastqc zip（zip 内含 html+数据，故未单收 html）
│   ├── 04_fastp/            # 12 对 fastp json+html
│   ├── 05_STAR/             # 6 份 Log.final.out + 6 份 SJ.out.tab
│   ├── 06_multiqc/          # multiqc 汇总报告（multiqc_report.html.gz + 解析数据表）
│   ├── 07_salmon_align/     # align 模式 6 个完整输出目录（小索引，quant.sf 仅数 KB）
│   ├── 08_salmon_direct/    # direct 模式 6 个完整输出目录
│   ├── 09_R/                # 3 张 PDF：火山 / GO dotplot / 热图（基因名正常版）
│   ├── scripts/             # 当轮运行脚本（命令原文 + 最小包装）
│   └── log/                 # 当轮 summary + deviations + R 分析日志
└── homework/                # homework/全量轮（2026-10-09，11.8万-22.7万读对/样本，全 GRCh38 索引，13 分钟）
    ├── 02_fastqc/           # 12 个 fastqc zip
    ├── 04_fastp/           # 12 对 fastp json+html
    ├── 05_STAR/             # 6 份 Log.final.out + 6 份 SJ.out.tab（太短率 47-53% 的出处）
    ├── 06_multiqc/          # multiqc 汇总报告（multiqc_report.html.gz + 解析数据表）
    ├── 07_salmon_align/     # align 模式 6 个目录——扣掉 10MB/份的 quant.sf（下游不用），保留
    │                        #   meta_info/偏差校正表/日志（映射率 100% 的证据在此）
    ├── 08_salmon_quasi/     # quasi 模式 6 个目录——quant.sf 转 .gz（10.6MB→2.1MB），其余原样
    │                        #   （映射率 30-35% 的证据；R 段实际输入）
    ├── 09_R/                # 4 张 PDF：火山 / GO dotplot / 热图（照原文：全 NA 证据）/
    │                        #   热图矫正版（0/50 NA）—— NA bug 证据对就在这里
    ├── scripts/             # 当轮运行脚本 + 热图矫正脚本 heatmap_corrected.R
    └── log/                 # 当轮 summary + deviations + R 分析日志
```

原始数据与重型中间产物（FASTQ、BAM、完整 salmon 目录、STAR 索引）保留在集群运行目录：
`/lustre1/user/bjx131_pkuhpc/rna_run_mini_20260921/`、`/lustre1/user/bjx131_pkuhpc/rna_run_homework_20261009/`。
