# Computational-Immuno-oncology

可复现的课堂内容仓库 / Reproducible course materials for computational immuno-oncology teaching.

> 当前内容：**WES（全外显子组测序）与 RNA-seq（转录组）两个模块的完整路线图与双轮实跑结果**。
> 后续将逐步加入免疫肿瘤学其他模块（TCR/BCR 分析、新抗原预测、肿瘤免疫微环境等）。

---

## 目录结构 / Layout

```text
.
├── WES/
│   ├── ref/                                  # 参考文档 (reference material)
│   │   ├── WES_roadmap_to_students.html      # 主路线图 (pandoc 导出, 图文完整)
    │   ├── WES_roadmap_to_students.pdf       # 同一文档的 PDF 版
    │   ├── WES_roadmap_english_embedded.html # 英文版 (图片内嵌)
    │   ├── WES_roadmap_corrected.html        # 矫正增强版：两轮实跑对照 + 逐坑修复说明
    │   ├── WES_roadmap_corrected.pdf         # 同一文档的 PDF 版
    │   ├── WES_roadmap_corrected_zh.html     # 同一矫正版的中文版 (正文全中文, 命令原样)
    │   ├── WES_roadmap_corrected_zh.pdf      # 中文版 PDF
    │   └── build.sh                           # pandoc+typst 构建脚本 (带表格/断行后处理与重叠自检)
    └── results/                              # 两轮实跑分析产物 (无任何原始数据: FASTQ/BAM/GVCF 不入库)
        ├── WES_results_guide.{md,html,pdf}  # 逐文件讲解: 怎么读、两轮对比、动手验证翻倍 bug 的命令
        ├── demo/                            # 0.5× demo 轮摘录 (fastp 报告/指标/VQSR 四件套/空 somatic/噪声 CNV/日志/脚本)
        └── full/                             # 全量轮摘录 (联合 55,123 + 顺序 VQSR 矫正/108 PASS somatic/amplicon vs hybrid CNV 对照/4 份 HsMetrics/8 个 sbatch/全部作业日志)
└── RNA-seq/                                 # RNA-seq 模块 (SEQC 双轮实跑完成)
    ├── rna_seq.pdf / RNA_seq_tutorial_mini.html  # 教程
    ├── RNA_seq_homework.html                # 作业
    ├── RNA_seq_code_explaination.pdf         # 代码讲解
    ├── ref/                                  # 矫正增强文档 (构建走 WES/ref/build.sh)
    │   ├── RNA_seq_tutorial_corrected{.md,.html,.pdf}      # 矫正增强版：两轮实跑对照 + 15 条坑位取证
    │   ├── RNA_seq_tutorial_corrected_zh{.md,.html,.pdf}   # 同一矫正版的中文版
    │   └── RNA_seq_results_guide{.md,.html,.pdf}           # 两轮结果逐类导读 + 7 个动手练习
    └── results/                              # 两轮实跑分析产物 (无任何原始数据: FASTQ/BAM 不入库)
        ├── mini/                             # mini 轮 (10k reads/样本, chr22 索引, ~4 分钟, 111/4989 转录本, DEG 45)
        └── homework/                        # 全量轮 (全 GRCh38 索引, 13 分钟, 831/231448 转录本, DEG 319, 含热图 NA bug 证据对与矫正版)
└── RNA-seq-application/                     # RNA-seq 下游应用模块 (第三节课, 实跑完成)
    ├── RNA-seq-application.html             # 教材原文 (Part 0-4, 33 个代码块)
    ├── ref/                                  # 矫正增强文档 (构建走 WES/ref/build.sh)
    │   ├── RNA_seq_application_corrected{.md,.html,.pdf}      # 矫正增强版：全流程实跑 + 7 条坑位取证
    │   ├── RNA_seq_application_corrected_zh{.md,.html,.pdf}   # 同一矫正版的中文版
    │   └── RNA_seq_application_results_guide{.md,.html,.pdf}  # 结果逐类导读 + 7 个动手练习
    └── results/                              # 实跑分析产物 (无任何原始数据)
        ├── 02_gene_correlation/             # Part 1: GZMA↔PRF1 / CTLA4↔CYT 相关散点
        ├── 03_survival_analysis/            # Part 2: KM 生存 + Cox 森林
        ├── 04_immune_infiltration/          # Part 3: 七法去卷积 (六法齐, MCPcounter 阻塞) + 堆叠柱/点/热图
        └── 05_tide/                         # Part 4: tidepy + 热图/柱图×2/ROC×2
└── TBCR-analysis/                          # TCR/BCR 组库下游分析模块 (第四节课, 实跑完成)
    ├── ref/
    │   ├── TBCR_analysis.html               # 教材原文 (Part 0-4, TRUST3 输出下游分析)
    │   ├── TBCR_analysis_corrected{.md,.html,.pdf}     # 矫正增强版：实跑 + 9 条坑位取证 + 修正版对拍
    │   ├── TBCR_analysis_corrected_zh{.md,.html,.pdf}  # 同一矫正版的中文版
    │   └── TBCR_analysis_results_guide{.md,.html,.pdf} # 结果逐类导读 + 7 个动手练习
    └── results/                              # 实跑分析产物 (无任何原始数据)
        ├── 01_clonotype/                    # Part 1: BCR consensus assembly 频率直方图
        ├── 02_gene_usage/                   # Part 2: IGHV/IGHJ top10 + V-J 热图
        ├── 03_cdr3_length/                  # Part 3: CDR3 nt/aa 长度分布
        ├── 04_diversity/                    # Part 4: 375 样本多样性 CSV + 三张相关图 (原样)
        └── 05_corrected/                    # 修正版: raw/L1/L2 三套指标 + TRB 扩展 + 对拍图
└── scRNA-seq/                              # 单细胞 RNA-seq 模块 (第五节课, 实跑完成)
    ├── scRNA_seq_analyze_totorial_2026.html  # 教材原文 (Part 0-8, 30 个代码单元, 无输出)
    ├── scRNA_seq_lecture_2026.pdf            # 课程幻灯片 (11 页)
    ├── *.png                                 # 28 张参考图 (教师 2025 输出, 未改动)
    ├── ref/                                  # 矫正增强文档 (构建走 WES/ref/build.sh)
    │   ├── scRNA_seq_analyze_corrected{.md,.html,.pdf}     # 矫正增强版: 原样跑 + 矫正跑 + 20 条坑位取证
    │   ├── scRNA_seq_analyze_corrected_zh{.md,.html,.pdf}  # 同一矫正版的中文版
    │   └── scRNA_seq_analyze_results_guide{.md,.html,.pdf} # 结果逐类导读 + 8 个动手练习
    └── results/                              # 实跑分析产物 (无任何原始数据)
        ├── 01_qc/                           # Part 2: 质控图 + 两份账目表 (运行值 / 原始矩阵重算值)
        ├── 02_hvg_pca/                      # Part 3-4: HVG 与 PCA 图 (Seurat/scanpy 两侧, 原样/修正两版)
        ├── 03_clustering/                   # Part 5: UMAP + 三种聚类方案的群大小表
        ├── 04_annotation/                   # Part 6: marker 图/表 + 两份硬编码映射审计 + 跨工具列联表
        ├── 05_deg/                          # Part 6: 火山图 + DEG 明细 + scanpy 排名/小提琴/点图
        ├── 06_integration/                  # Part 7: harmony 前后 UMAP + 混合度指标
        ├── 07_objects/                      # Part 8: pbmc3k.h5ad + 逐细胞元数据
        ├── 08_figure_comparison/            # 与 28 张参考图的机器可读/书面比对
        ├── 09_run_logs/                     # 六次运行的逐单元日志与状态表
        └── scripts/                          # 13 个脚本: 教材单元/驱动器/审计与比对
```

> **矫正增强版（WES_roadmap_corrected）**：把学生版路线图在集群上**原样跑两轮**（0.5× 降采样 demo 轮 + 全量轮）后，对照实测结果写成的修复与讲解版——每个步骤含原理/目的/实跑结果/坑与矫正/结果解读/调优扩展，附录带完整问题清单（12 条，含证据）与一页矫正版命令速查。已知问题速览：VQSR MergeVcfs 每条变异重复两遍（含老师参考产物）、ANNOVAR 富协议库名不存在、CNVkit amplicon 模式用于杂交捕获 panel、无靶区深度度量（实测全量 OC 仅 54.6× 可用）。
> **中文版（WES_roadmap_corrected_zh）**：与英文矫正版内容一致、互为对照；正文全中文，命令/路径/参数原样保留，数字逐一对拍核验。构建：`./build.sh WES_roadmap_corrected_zh.md "WES 全流程分析路线图 · 矫正增强版（中文版）"`（英文版：`./build.sh`）。
>
> **分析结果（WES/results/）**：两轮实跑的分析产物摘录（QC 报告、比对/去重指标、VQSR 模型与 tranches、变异集 VCF、体细胞检出、CNV 分段、HsMetrics、全部实跑脚本与作业日志）——**不含任何原始数据**。讲解文档 `WES/results/WES_results_guide.html/pdf` 逐类讲怎么读，含 6 个动手练习（比如用仓库内文件三条命令亲手验证 VQSR 翻倍 bug）。大文件（BAM/GVCF/富注释表）仍在集群运行目录，讲解文档 §5 有清单。

> **矫正增强版（RNA_seq_tutorial_corrected）**：把 mini 教程与 homework 作业在集群上**原样跑两轮**（mini 轮 2026-09-21 + 全量轮 2026-10-09，均按教程原生方式登录节点串行跑）后写成的修复与讲解版——每步含原理/实跑数字/坑与矫正/调优扩展，带 15 条问题登记表（含热图全 NA 沉默 bug、''' 伪注释、双模式映射率口径等，每条带证据）。中文版 `_zh` 与英文版数字逐一对拍一致。构建：`WES/ref/build.sh ../../RNA-seq/ref/RNA_seq_tutorial_corrected.md "RNA-seq 全流程教程 · 矫正增强版"`（中文版传 `_zh.md` 与中文标题）。
>
> **分析结果（RNA-seq/results/）**：两轮实跑的分析产物（fastqc/fastp 报告、STAR 指标与 SJ 表、multiqc、salmon 双模式定量与映射率证据、R 段四张 PDF 含**热图 NA bug 证据对**、全部脚本与日志）——**不含任何原始数据**。讲解文档 `RNA-seq/ref/RNA_seq_results_guide.html/pdf` 逐类导读 + 7 个动手练习（练习 1 三条命令在仓库内亲手复现热图 50×NA bug）。大文件（FASTQ/BAM/完整 salmon 目录）仍在集群运行目录，讲解文档 §4 有清单。

> **矫正增强版（RNA_seq_application_corrected）**：把第三节课教材（Van Allen 2015 黑色素瘤队列，基因相关性/生存/免疫浸润/TIDE 四个 Part）在登录节点**原样实跑一轮**（2026-10-12，纯 R + 一行 tidepy，~6 分钟）后写成的修复与讲解版——每 Part 含原理/实跑数字/坑与矫正/调优扩展，带 7 条问题登记表（含 MCPcounter conda 1.1.0 http 硬编码阻塞、EPIC 收敛警告、tidepy 三行预期 WARN 等，每条带证据）。中文版 `_zh` 与英文版数字逐一对拍一致。构建：`WES/ref/build.sh ../../RNA-seq-application/ref/RNA_seq_application_corrected.md "RNA-seq Application Analysis · Corrected & Annotated Edition"`（中文版传 `_zh.md` 与中文标题）。
>
> **分析结果（RNA-seq-application/results/）**：实跑分析产物（Part 1 相关散点、Part 2 KM+森林、Part 3 七法去卷积表+四张 quanTIseq 图、Part 4 TIDE 输出+热图+柱图×2+ROC×2、全部脚本与日志）——**不含任何原始数据**（表达矩阵/随访表在共享只读路径，未复制）。讲解文档 `RNA-seq-application/ref/RNA_seq_application_results_guide.html/pdf` 逐类导读 + 7 个动手练习（练习 1 亲手复现 MCPcounter http vs https 阻塞）。本模块无大文件，全部产出已入库。

> **矫正增强版（TBCR_analysis_corrected）**：把第四节课教材（TRUST3 输出的 TCR/BCR 组库下游分析：克隆型/基因使用/CDR3 长度/多样性指数）在登录节点**原样实跑一轮**（2026-10-19，纯 R，约 45 秒）后，再**补跑一份修正版**对拍写成的修复与讲解版——每 Part 含原理/实跑数字/坑与矫正/调优扩展，带 9 条问题登记表（核心三条：`nrow()` 冒充唯一 CDR3、`partial` 黑名单漏掉 `out_of_frame`、队列循环跳过 Part 3 的质量过滤；另有插补序列、单克隆 `NaN`、log2 vs 自然对数、`*` 占位符当基因等）。中文版 `_zh` 与英文版数字逐一对拍一致（218 个数字 token 出现次数完全相同）。构建：`WES/ref/build.sh ../../TBCR-analysis/ref/TBCR_analysis_corrected.md "TCR/BCR Repertoire Downstream Analysis · Corrected & Annotated Edition"`（中文版传 `_zh.md` 与中文标题）。
>
> **分析结果（TBCR-analysis/results/）**：原样四部分产物（克隆型频率直方图、IGHV/IGHJ top10、V-J 热图、CDR3 nt/aa 长度分布、375 样本多样性 CSV 与三张相关图、全部脚本与日志）+ `05_corrected/` 修正版（raw/L1/L2 三套逐样本指标、TRB 扩展 CSV、全部修正图与 `corrected_summary.txt`）——**不含任何原始数据**（TRUST3 产物在共享只读路径 `/lustre1/share/airr/data/`，未复制）。讲解文档 `TBCR-analysis/ref/TBCR_analysis_results_guide.html/pdf` 逐类导读 + 7 个动手练习（练习 2 五条命令现场复现 raw/L1/L2 三口径差异）。

> **矫正增强版（scRNA_seq_analyze_corrected）**：把第五节课教材（PBMC 3k：R/Seurat 的 Part 1-7 + Python/scanpy 的 Part 8，共 30 个代码单元）在登录节点**原样实跑一轮**（2026-10-26，一次 R 会话 509 s + 一次 Python 命名空间 4 s）后，再**跑一份修正版**对拍写成的修复与讲解版——每 Part 含原理/实跑数字/坑与矫正/调优扩展，带 20 条问题登记表。两条致命：R 侧从未 `library(ggrepel)`（DEG 单元连带 `saveRDS` 全灭）、Python 侧 `log=True,,save=` 双逗号 `SyntaxError` 抹掉整个单元（PCA/近邻/UMAP 全没跑）；**原样跑只产出 25 张图里的 18 张、2 个对象里的 0 个**。改变结果的几条：教材把 cluster→细胞类型映射写了两份且互相矛盾（cell 39 的表与 cell 43 的向量在 cluster 1/2 上互换，本次运行里向量对、表错）；scanpy 的位置式命名 8 个标签错 2 个（名为 `CD8 T` 的群 78.1% 是 Seurat 的 Naive CD4 T）；两侧 QC 阈值不一致（2500 vs 2000）正是"R 9 群 vs Python 8 群"的主因——把 Python 阈值提到 2500 就得到与 Seurat 完全相同的 2,638 个细胞和 9 个群。另记录：基因全集改动（`min_cells=3`，只有 30/2000 个 HVG 变化）会让同一份固定注释向量标错 957 个细胞（36%）；Part 7 的 UMAP 不可跨运行复现（同输入两次运行结构相关仅 0.563，与对参考图同量级）；参考图 `Marker_Heatmap.png` 是 3000×6000 px，而教材 `ggsave` 只能给 3000×3000。中文版 `_zh` 与英文版结构、数字逐一对拍。构建：`WES/ref/build.sh ../../scRNA-seq/ref/scRNA_seq_analyze_corrected.md "Single-cell RNA-seq (PBMC 3k) — Corrected & Annotated Edition"`（中文版传 `_zh.md` 与中文标题）。
>
> **分析结果（scRNA-seq/results/）**：原样跑与修正跑的图（质控、HVG/PCA、聚类 UMAP、marker 与注释、火山图与 DEG、harmony 前后 UMAP、scanpy 的 11 张）＋ 教材算了却从未落盘的表（11,986 行 marker、90 行每群 top10、两张 DEG 明细、每群前 200 的 scanpy marker）＋ 两份硬编码映射的审计表与跨工具列联表 ＋ 与 28 张参考图的逐图像素比对（`PCA_heatmap_ALL.png` 结构相关 **1.0000**，即逐像素可复现）＋ 六次运行的逐单元日志 ＋ `pbmc3k.h5ad` 交付对象——**不含任何原始数据**（计数矩阵在共享只读路径 `/lustre1/share/scRNA_seq/hg19`，未复制；331 MiB 的 `pbmc3k_final.rds` 与 120 MiB 的合并对象不入库，各一行 `saveRDS()` 可重建）。讲解文档 `scRNA-seq/ref/scRNA_seq_analyze_results_guide.html/pdf` 逐类导读 + 8 个动手练习（练习 3 一条 awk 命令现场证伪 `CD8 T` 标签）。

## WES 路线图内容一览 / What the roadmap covers

一条 13 步、完全可照抄执行的 WES 分析管线（ tumor: **OC** 卵巢癌组织 ｜ matched normal: **PBMC** 外周血 ）：

| # | 步骤 | 工具 | 核心概念 |
| --- | ------ | ------ | ---------- |
| 0 | 环境与约定 | conda (`wes` env) | 参考基因组一致性 (GRCh38/hg38)、目录约定 |
| 1 | FASTQ 格式 | zcat/awk | Phred 质量值 Q = −10·log₁₀(P_error) |
| 2 | 测序质量控制 | fastp | 接头检测 (PE overlap)、滑窗截取、Q30 报告 |
| 3 | 比对 | BWA-MEM | seed-and-extend、read group (SM 字段!) |
| 4 | 标记 PCR 重复 | GATK MarkDuplicates | 未剪切坐标判重、只标记不删除 |
| 5 | 碱基质量重校正 | GATK BQSR | 已知位点建模系统误差、调整"可信度"而非序列 |
| 6 | 胚系变异检测 | HaplotypeCaller | 局部单倍型重组装、de Bruijn 图、PairHMM、GVCF |
| 7 | 变异质量重校正 | VQSR | 注释空间高斯混合模型、tranche 敏感度 |
| 8 | 功能注释 | ANNOVAR | 基因结构、氨基酸改变、人群频率、致病预测 |
| 9 | 体细胞 TMB | Mutect2 | tumor-normal 配对、链取向偏倚校正；TMB = 合格体细胞变异 / 可靶向 Mb |
| 10 | 拷贝数变异 | CNVkit | bin 深度归一、log2 ratio、CBS 分段 |
| 11 | 结果解读清单 | — | Ti/Tv、覆盖率、基因组版本一致性等 9 项检查 |
| 12–13 | 文件与总结 | — | 每类产物的含义 |

每一步都先讲**算法为什么这样设计**，再给可直接运行的命令。

## 在教学集群上运行 / Running on the teaching cluster

> 本仓库只存放**文档与脚本**；测序数据与个人分析结果不入库。

```bash
ssh <teaching-node>
export SHARE=/lustre1/share
source "$SHARE/miniconda3/etc/profile.d/conda.sh"
conda activate wes            # fastp / bwa / samtools / bcftools / gatk / cnvkit

# 参考数据（集群共享、只读）
export REF=$SHARE/references                          # hg38.fa + .fai/.dict + BWA index
export TARGETS=$REF/S07604514_AllTracks_V6_60_hg38.bed  # 243,559 intervals, 60.51 Mb
export ANNOVAR=$SHARE/annovar_new                     # table_annovar.pl + humandb (hg38)

# 教学数据（降采样 demo 版，输出写到自己目录）
#   $SHARE/data_new/{OC,PBMC}             # 新降采样 demo（新版路线图指向这个）
#   $SHARE/data/OC_WES/{OC,PBMC}          # 旧降样本 (弃用，README 旧版指向这个)
#   $SHARE/data/OC_WES_all/{OC,PBMC}      # 全量 (OC 10.4 Gbp, PBMC 15.0 Gbp)

# 个人输出目录
export BASE=/gpfs1/home/$USER/wes_run_cluster_$(date +%Y%m%d)
mkdir -p "$BASE"/{0_fastq,1_alignment/OC,1_alignment/PBMC,2_variant_call,3_annotation,4_cnv,5_somatic_tmb}
```

按路线图第 1–10 节顺序执行即可；推荐用 SLURM (`sbatch`) 提交 BWA/GATK 等重计算步骤，不要在登录节点直接跑。

后四节课（RNA-seq、RNA-seq 下游应用、TCR/BCR 组库、单细胞 RNA-seq）都是**登录节点串行直跑**，与 WES 的 sbatch 纪律不同：`conda activate RNA-seq` / `rnaimmune` / `airr` / 直接用教师共享环境 `/lustre1/share/miniconda3/envs/scRNA_seq`（R+Seurat 与 Python+scanpy 同一环境）后，按各自 `ref/` 下的矫正版文档执行即可。

## 致谢 / Credits

- 路线图作者：**Jiahao Ma (马嘉昊)**, **Zexian Zeng (曾泽贤)** — 北京大学前沿交叉学科研究院 (Center of Life Sciences & Center for Quantitative Biology)
  原文面向课堂分发；本仓库仅作教学存档与复现配套，版权归原作者。
- 感谢 Yufeng He 此前的教学路线图工作。

## 仓库约定 / Repo conventions

- 大文件（HTML/PDF 教学文档 ≈16MB）有意直接入库，便于学生一次 clone 拿到全部讲义；其余大文件一律走 [Git LFS] 或集群共享目录。
- `.gitignore` 屏蔽所有测序数据与运行产物 (FASTQ/BAM/VCF/…) —— **任何原始数据、患者数据、凭据不入库**。
- 提交时 pre-commit 钩子会自动扫描 API key / 密钥模式（含 base64 感知，避免教学 HTML 中的嵌入图片误报）。
