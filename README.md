# Computational-Immuno-oncology

可复现的计算免疫肿瘤学课堂仓库 / Reproducible course materials for computational immuno-oncology teaching.

**八节课全部跑完并收口**（原样跑 → 问题登记 → 矫正链 → 双语矫正增强文档 → results 入库 → 与教师参考图逐图定量对拍），
累计 **117 条带证据的问题登记**、**24 份矫正版/结果导读**（英文矫正版 / 中文矫正版 / 结果导读，共 **501 页**，
另有 WES 主路线图 28 页）、**两套教师参考图的逐图对拍**。仓库只存文档、脚本与精选分析产物，**任何原始数据与凭据都不入库**。

最后更新：2026-11-10（第八节 CellphoneDB 收口，commit `f964b50`）。

---

## 快速上手 / Quick start

| 你想做什么 | 去哪儿 |
|---|---|
| 照着跑一遍某节课 | 该模块 `ref/` 下的**矫正增强版**（英文 `*_corrected.md/pdf`，中文 `*_corrected_zh.*`）——每步含原理、目的、命令讲解、结果解读、常见坑、参数调优、扩展用法 |
| 只想看结果怎么读 | 该模块 `ref/` 下的**结果导读**（`*_analysis_results_guide.pdf`）——逐目录导读 + 动手练习 |
| 想知道这节课哪里有问题 | 矫正增强版里的**问题登记表**（编号如 `DYN-01`、`CPB-19`），每条都指向证据文件 |
| 想复现某个数字 | `results/` 下对应目录的 TSV + `results/scripts/` 下的脚本 + `results/09_run_logs/` 下的控制台 |
| 想比对教师参考图 | `results/07_figure_comparison/`（结构相关 + 墨迹 IoU，逐图一行） |

---

## 模块总览 / Modules at a glance

| # | 模块目录 | 课程内容 | 教材单元 | 原样跑实测 | 问题登记 | 文档页数（EN / ZH / 导读） |
|---|---|---|---|---|---|---|
| 1 | `WES/` | 全外显子组 13 步路线图（OC 卵巢癌 + PBMC 配对） | 13 步 | 两轮（0.5× demo + 全量） | 12 | 35 / 30 / 9 |
| 2 | `RNA-seq/` | 转录组全流程（SEQC；mini + homework 双轮） | 教程 + 作业 | mini ≈4 min、全量 13 min | 15 | 20 / 19 / 8 |
| 3 | `RNA-seq-application/` | 转录组下游应用（Van Allen 2015 黑色素瘤，4 个 Part） | 33 代码块 | ≈6 min | 7 | 23 / 22 / 10 |
| 4 | `TBCR-analysis/` | TCR/BCR 组库下游分析（TRUST3 输出，4 个 Part） | 4 Part | ≈45 s | 9 | 28 / 25 / 12 |
| 5 | `scRNA-seq/` | 单细胞 RNA-seq（PBMC 3k；Seurat 7 个 Part + scanpy 1 个） | 30 单元 | R 509 s + Python 4 s | 20 | 38 / 36 / 13 |
| 6 | `scRNA-seq-scVelo/` | RNA velocity（小鼠胰腺 E15.5，3,696 × 27,998 基因） | 29 单元 | 42.6 s | 17 | 31 / 27 / 4 |
| 7 | `scRNA-seq-application-dynamo/` | Dynamo 矢量场（同一份胰腺数据，7 个 Part） | 34 单元 | ≈12 min | 18 | 31 / 27 / 6 |
| 8 | `scRNA-seq-cellchat/` | CellphoneDB 细胞通讯（蜕膜–胎盘界面，3,312 细胞） | 34 单元 | ≈100 s | 19 | 23 / 20 / 4 |

> 第八节课的目录叫 `cellchat`，内容是 **CellphoneDB**（全文无一行 R，`CellChat` 一次未出现）——改名会破坏与教师材料、已入库产物的对应关系，故保留原名并在文档 §1.1 说明。

---

## 目录结构 / Layout

```text
.
├── WES/                                        # 全外显子组路线图 (第一讲, 双轮实跑完成)
│   ├── ref/                                    #   参考文档 + 矫正增强版 + 构建脚本
│   │   ├── WES_roadmap_to_students.{html,pdf}  #   主路线图 (pandoc 导出, 图文完整)
│   │   ├── WES_roadmap_english_embedded.html   #   英文版 (图片内嵌)
│   │   ├── WES_roadmap_corrected{.md,.html,.pdf}      #   矫正增强版: 两轮实跑对照 + 12 条坑位取证
│   │   ├── WES_roadmap_corrected_zh{.md,.html,.pdf}   #   同一矫正版的中文版
│   │   └── build.sh                            #   pandoc+typst 构建脚本 (表格/断行后处理 + 重叠自检)
│   └── results/                                #   两轮实跑产物 (无原始数据: FASTQ/BAM/GVCF 不入库)
│       ├── WES_results_guide.{md,html,pdf}     #   逐文件讲解: 怎么读、两轮对比、动手验证翻倍 bug
│       ├── demo/                               #   0.5× demo 轮摘录 (fastp/指标/VQSR 四件套/空 somatic/噪声 CNV)
│       └── full/                               #   全量轮摘录 (联合 55,123 + 顺序 VQSR/108 PASS somatic/amplicon vs hyb)
├── RNA-seq/                                    # 转录组全流程 (第二讲, 双轮实跑完成)
│   ├── rna_seq.pdf / RNA_seq_tutorial_mini.html / RNA_seq_homework.html / RNA_seq_code_explaination.pdf
│   ├── ref/                                    #   RNA_seq_tutorial_corrected{,_zh} / RNA_seq_results_guide
│   └── results/                                #   mini/ (~4 min, DEG 45) 与 homework/ (13 min, DEG 319)
├── RNA-seq-application/                        # 转录组下游应用 (第三讲, 实跑完成)
│   ├── RNA-seq-application.html                #   教材原文 (Part 0-4, 33 个代码块)
│   ├── ref/                                    #   RNA_seq_application_corrected{,_zh} / RNA_seq_application_results_guide
│   └── results/                                #   02_gene_correlation/ 03_survival_analysis/ 04_immune_infiltration/ 05_tide/
├── TBCR-analysis/                              # TCR/BCR 组库下游 (第四讲, 实跑完成)
│   ├── ref/                                    #   教材原文 + TBCR_analysis_corrected{,_zh} / TBCR_analysis_results_guide
│   └── results/                                #   01_clonotype/ 02_gene_usage/ 03_cdr3_length/ 04_diversity/ 05_corrected/
├── scRNA-seq/                                  # 单细胞 RNA-seq (第五讲, 原样跑 + 矫正跑完成)
│   ├── scRNA_seq_analyze_totorial_2026.html    #   教材原文 (Part 0-8, 30 个代码单元, 无输出)
│   ├── scRNA_seq_lecture_2026.pdf              #   课程幻灯片
│   ├── *.png                                   #   28 张教师参考图 (未改动)
│   ├── ref/                                    #   scRNA_seq_analyze_corrected{,_zh} / scRNA_seq_analyze_results_guide
│   └── results/                                #   01_qc … 08_figure_comparison / 09_run_logs / scripts
├── scRNA-seq-scVelo/                           # RNA velocity (第六讲, 原样跑 + 矫正链完成)
│   ├── scRNAseq-scVelo.html                    #   教材原文 (88 单元含输出)
│   ├── 01_proportions.png … 11_*.png           #   11 张课程随附参考图 (未改动)
│   ├── ref/                                    #   教材 .ipynb + scVelo_analysis_corrected{,_zh} / guide / figures/(21)
│   └── results/                                #   01_preprocess … 07_figure_comparison / 09_run_logs / scripts
├── scRNA-seq-application-dynamo/                # Dynamo 矢量场 (第七讲, 原样跑 + 矫正链完成)
│   ├── scRNAseq-application_dynamo.html        #   教材原文 (50 单元 = 16 md + 34 code)
│   ├── ref/                                    #   自带 build.sh + corrected{,_zh} / guide / figures/(14+11+6)
│   └── results/                                #   01_preprocess … 07_figure_comparison / 09_run_logs / scripts
└── scRNA-seq-cellchat/                         # CellphoneDB 细胞通讯 (第八讲, 原样跑 + 矫正链完成)
    ├── scRNAseq-cellchat.html                  #   教材原文 (94 单元 = 60 md + 34 code)
    ├── ref/                                    #   自带 build.sh + corrected{,_zh} / guide / figures/(12+3+15)
    └── results/                                #   01_input … 07_degs / 09_run_logs / scripts
```

命名约定：`*_corrected.md` 英文矫正增强版，`*_corrected_zh.md` 中文版，`*_analysis_results_guide.md` / `*_results_guide.md` 结果导读；
`results/NN_topic/` 按分析阶段编号，`results/07_figure_comparison/`（WES 为 `08_`、`results/figures/`）放对拍表，
`results/09_run_logs/` 放各次运行的完整控制台，`results/scripts/` 放可重跑脚本。

---

## 八节课逐讲要点 / Lesson notes

每讲的完整论述在其 `ref/` 文档中；这里只留"这一讲到底发现了什么"。

### 1. WES 全外显子组路线图（第一讲）

- OC 卵巢癌组织 + 配对 PBMC，13 步：fastp → BWA-MEM → MarkDuplicates → BQSR → HaplotypeCaller → VQSR → ANNOVAR → Mutect2/TMB → CNVkit → 解读清单。
- 双轮实跑（0.5× demo + 全量）后发现的主要问题：**VQSR `MergeVcfs` 让每条变异重复两遍**（教师参考产物里也有）、
  ANNOVAR 富协议库名不存在、CNVkit **amplicon 模式用于杂交捕获 panel**、**缺靶区深度度量**（全量轮 OC 实测仅 54.6× 可用）。
- 全量轮：联合 55,123 变异、顺序 VQSR、108 个 PASS 体细胞变异。

### 2. RNA-seq 转录组全流程（第二讲）

- mini 教程（10k reads/样本、chr22 索引、≈4 min、111/4,989 转录本、DEG 45）与 homework 全量轮（GRCh38、13 min、831/231,448、DEG 319）都按教程原生方式在登录节点串行跑——两份教程 HTML 里没有任何 `sbatch`/`srun` 痕迹。
- 主要问题：**热图全 NA 的静默 bug**、`'''` 伪注释、STAR/salmon 双模式映射率口径不一致。矫正版给出可复现的对照与修法。

### 3. RNA-seq 下游应用（第三讲，Van Allen 2015）

- 四个 Part：基因相关性（GZMA↔PRF1、CTLA4↔CYT）、KM 生存 + Cox 森林、七法免疫浸润去卷积、TIDE 预测。
- 主要问题：**MCPcounter conda 1.1.0 的 http 硬编码**导致该方法整体阻塞（六法齐、一法缺）、EPIC 的收敛警告、tidepy 预期内的三行 WARN。
- 图文要点：quanTIseq 堆叠柱里"未分类细胞占六到九成"、CTLA4–CYT 的 Pearson/Spearman 分裂、TIDE 预测柱与真实柱并排会读漏判。

### 4. TCR/BCR 组库下游分析（第四讲，TRUST3）

- 四个 Part：克隆型频率、IGHV/IGHJ 使用与 V-J 热图、CDR3 nt/aa 长度、375 样本多样性指数。
- 原文三条核心 bug：`nrow()` 冒充唯一 CDR3、`partial` 黑名单漏掉 `out_of_frame`、队列循环跳过 Part 3 的质量过滤；
  另有插补序列、单克隆 `NaN`、`log2` 与自然对数混用、`*` 占位符当基因名。
- 修正版给出 raw/L1/L2 三套多样性口径与 TRB 扩展；EN–ZH 218 个数字 token 计数完全一致。

### 5. 单细胞 RNA-seq（第五讲，PBMC 3k）

- Seurat Part 1–7 + scanpy Part 8，30 个代码单元。**原样跑只产出 25 张图里的 18 张、2 个对象里的 0 个**：
  R 侧从未 `library(ggrepel)`（DEG 单元连带 `saveRDS` 全灭）、Python 侧 `log=True,,save=` 双逗号 `SyntaxError` 抹掉整个单元。
- 会改变结论的几条：cluster→细胞类型映射写了两份且互相矛盾；scanpy 的位置式命名 8 个标签错 2 个（名为 `CD8 T` 的群 **78.1 %** 是 Naive CD4 T）；
  两侧 QC 阈值不一致（2500 vs 2000）正是"R 9 群 vs Python 8 群"的主因——Python 阈值提到 2500 后得到与 Seurat 完全相同的 2,638 细胞 / 9 群。
  另：基因全集改动会让固定注释向量标错 **957 个细胞（36 %）**，Part 7 的 UMAP 不可跨运行复现（结构相关仅 0.563）。

### 6. RNA velocity（第六讲，scVelo / 小鼠胰腺 E15.5）

- **教材用 9 张幻灯讲 dynamical（EM）模型却从未调用 `recover_dynamics`**，只跑了稳态近似。补跑后同一批数据
  stochastic vs dynamical 的 top-20 速率基因 Jaccard 仅 **0.132**，置信度相关 **0.207**，逐基因 γ 秩相关 **−0.45**。
- 另有：`recover_dynamics` 默认 `var_names="velocity_genes"` 会把 EM 限制在上一次的 1,038 个基因（改 `"all"` 后 2,000 个）；
  `scv.pl.velocity(mode='dynamical')` 只控侧栏不选模型（两图 md5 相同）；`score_genes_cell_cycle` 静默覆盖自带 S/G2M 分数（Pearson 0.964/0.973）；
  `scv.pl.proportions` 默认画过滤前文库（16.3 % vs 23.0 %）；`min_corr=.3` 从未生效；预登记的"双重取对数"被**证伪**。
- 11 张参考图里十张结构相关 0.973–0.993，仅细胞周期那张 0.589（参考侧 dpi 分 80/96/64 三档）。

### 7. Dynamo 矢量场（第七讲，同一份胰腺数据）

- `reduceDimension` 是**静默 no-op**（`[UMAP] completed [0.0012s]`）——PCA 由 `recipe_monocle` 新算、UMAP 是文件自带的、kNN 图文件里根本没有；
  强制重算会把细胞挪走 1.9 倍 UMAP 展布，并能**翻转 Ductal 旋度符号**，故矫正只在副本上算。
- 三个 HVG 标记互不认账（`use_for_pca` 4,000 与文件自带 categorical 只重叠 1,339，Jaccard 0.201）；`model="stochastic"` 把 β 固定为 1.0。
- 教材 cell 24 的**注释对、代码错**：注释写 Pdx1→Ins2，代码画 Pax4→Ins2；实测 Pdx1→Ins2 在 Beta 群 **+1.5e-3、99.3 % 细胞为正**（全表最强），
  Pax4→Ins2 在 Beta 是 −4e-5（46 %）。Pax4⤙Arx 的"互抑"只发生在 Ngn3 high EP，到 Alpha/Beta 转为互激——是**局部开关**。
- 四句定性判断里"**Ductal 是散度意义上的鞍点**"不成立（Ductal 显著更不负，MWU p=6.5e-06；它真正是最慢 + 潜势最低）；
  "Ngn3 high EP 是分叉点"与"终末群收缩"完全成立。状态图矩阵**不是行随机**，只能读作通量。教师 2026 参考集里 GO 两张图是**纯白空画布**。
- 可复现性：两次运行相隔 2.5 h，**6/6 图 md5 逐字节相同**。

### 8. CellphoneDB 细胞通讯（第八讲，蜕膜–胎盘界面）

- 原样跑 ≈100 s（简单打分 27 s、1000 次置换 41 s、打分 22 s、两张点图 ≈65 s）。**2,440 互作 × 81 有序对 = 197,640 格**，
  显著 4,248 格 / 417 互作；矩阵里 **96 % 的格 p = 1**——不是"不显著"，而是**没过 10 % 可用性掩码、根本没测**。
- `"PV MMPP11"` 多一个 `P` → ktplotspy **静默丢类型**（24 列而非 36 列），而**教师 2026 参考图同样是 24 列**；
  `list(...).sort() == list(...).sort()` 是 `None == None`，barcode 校验永远 `True`。
- `threads=5` **不可复现**（单线程两次 md5 逐字节相同、4,246 格；5 线程 4,251 且每次不同）；
  197,640 次检验无多重校正（α=0.05 期望假阳性 9,882）；**BH–FDR 真算后 q ≤ 0.05 只剩 2,946 格 / 311 互作**，且正好是 p = 0 的那些格。
- 阈值是**覆盖度掩码**：0.05/0.10/0.20 → 7,007/708、4,251/417、2,413/218，而最小的簇只有 **7 个细胞**，10 % 等于"1 个细胞"。
- **`significant_means` 与 `means`/`pvalues` 行序完全不同**（同行同 id 仅 0.04 %）——按位置合并会把
  `COL4A1_integrin_a1b1` 的 76 挂到 `IAPP_AMY3_receptor`（真实 0）头上；审计表 `01_input/table_row_order.tsv` 把这一点钉成可核证据。
- DEG 分支（教材全程 `degs_analysis=False`）真跑得 4,191 格 / 419 互作，与置换法 Jaccard **0.469**（格）/ **0.742**（互作）；
  附带的 DEG 表只覆盖 9 个类型中的 6 个——三个 PV 类型一个都没有。

---

## 交付规范 / Delivery pipeline

每节课都走同一条六段管线，顺序不可颠倒：

1. **原样跑通**：把教材的代码单元逐格还原执行，参数一字不改；只允许最小偏离（shell magic 中和、硬编码路径替换），且逐条记录在 `patches_applied.txt`。
2. **问题登记**：执行期**不修**，只记录（编号 + 严重度 + 证据文件）。校正留到第 3 段。
3. **矫正链**：C1…Cn 逐项修复，每项产出可对拍的文件；文档写原理、目的、命令讲解、结果解读、常见坑、参数调优、扩展用法。
4. **results 入库**：按 `NN_topic/` 归档，配 `09_run_logs/` 与 `scripts/`；**不含原始数据**。
5. **仓库集成**：`.gitignore` 例外区 + pre-commit 放行 + `.pi-lens.json` + 本 README。
6. **图对拍**：与教师参考图逐图算**结构相关**与**墨迹 IoU**，两套参考（HTML 内嵌 vs 脚本重跑）分别记录。

**两条铁律**：原样跑不因"已知有错"而提前纠正；文档里的每个数字都必须能指到 `results/` 里的文件或 `09_run_logs/` 里的控制台。

---

## 在教学集群上运行 / Running on the teaching cluster

> 本仓库只存放**文档与脚本**；测序数据与个人分析结果不入库。

**WES（唯一走 sbatch 的模块）**：

```bash
ssh <teaching-node>
export SHARE=/lustre1/share
source "$SHARE/miniconda3/etc/profile.d/conda.sh"
conda activate wes            # fastp / bwa / samtools / bcftools / gatk / cnvkit

export REF=$SHARE/references                          # hg38.fa + .fai/.dict + BWA index
export TARGETS=$REF/S07604514_AllTracks_V6_60_hg38.bed  # 243,559 intervals, 60.51 Mb
export ANNOVAR=$SHARE/annovar_new                     # table_annovar.pl + humandb (hg38)

# 教学数据（降采样 demo 版；输出写到自己目录）
#   $SHARE/data_new/{OC,PBMC}             # 新版路线图指向这个
#   $SHARE/data/OC_WES/{OC,PBMC}          # 旧降样本 (弃用)
#   $SHARE/data/OC_WES_all/{OC,PBMC}      # 全量 (OC 10.4 Gbp, PBMC 15.0 Gbp)
```

BWA/GATK 等重计算一律 `sbatch` 提交到计算分区，登录节点只提交与查看。

**其余七讲（登录节点串行直跑）**：两份转录组教程 HTML 里没有任何 `sbatch`/`srun` 痕迹，教师 30 人并发负载测试也是在登录节点做的，故这些课按教程原生方式直跑。

| 讲 | 环境 | 说明 |
|---|---|---|
| 2 RNA-seq | `conda activate RNA-seq` | STAR / salmon / R |
| 3 RNA-seq-application | `conda activate rnaimmune` | 纯 R + 一行 tidepy |
| 4 TBCR-analysis | `conda activate airr` | 纯 R |
| 5 scRNA-seq | `/lustre1/share/miniconda3/envs/scRNA_seq` | R+Seurat 与 Python+scanpy 同一环境 |
| 6 scVelo | 新建 `velo` | scvelo 0.3.4 / scanpy 1.10.4；共享环境里没有 scVelo，硬塞会连坐升级 anndata/numba |
| 7 Dynamo | `/lustre1/share/scRNA_seq_application/lesson/dynamo` | 教师原环境（dynamo-release 1.3.2 / py3.11） |
| 8 CellphoneDB | `/lustre1/share/scRNA-adv-envs/scRNA-adv` | cellphonedb 5.0.1 / ktplotspy 0.3.3 |

共享只读数据（**不复制进仓库**）：WES 见上方 `$SHARE`；TCR/BCR `$SHARE/airr/data/`；单细胞 `$SHARE/scRNA_seq/hg19`；
胰腺数据 `/lustre1/share/scRNA-advancedTopics/data/Pancreas/endocrinogenesis_day15.h5ad`；
CellphoneDB `/lustre1/share/scRNA-adv-2026/data/data_cellphoneDB/`。其余数据位置见各模块文档的"复现"一节。

---

## 构建文档 / Building the documents

所有文档都走 pandoc 3.11 + typst 0.15.1 的同一套管线：表格跨页后处理、长 token 断行（真 U+200B）、
窄列缩字号、`pdftotext -bbox` 词框重叠自检（失败即非零退出）。`build.sh` 会 `cd` 到自己所在目录，故必须在该模块的 `ref/` 下运行。

| 模块 | 构建（在 `ref/` 下） | 产物页数（EN / ZH / 导读） |
|---|---|---|
| WES | `./build.sh` | 35 / 30 / 9（导读在 `results/`） |
| RNA-seq | `../WES/ref/build.sh RNA_seq_tutorial_corrected.md "<title>"` | 20 / 19 / 8 |
| RNA-seq-application | `../WES/ref/build.sh RNA_seq_application_corrected.md "<title>"` | 23 / 22 / 10 |
| TBCR-analysis | `../WES/ref/build.sh TBCR_analysis_corrected.md "<title>"` | 28 / 25 / 12 |
| scRNA-seq | `../WES/ref/build.sh scRNA_seq_analyze_corrected.md "<title>"` | 38 / 36 / 13 |
| scVelo | `../WES/ref/build.sh scVelo_analysis_corrected.md "<title>"` | 31 / 27 / 4 |
| Dynamo | `./build.sh scRNAseq_dynamo_corrected.md "<title>"` | 31 / 27 / 6 |
| CellphoneDB | `./build.sh scRNAseq_cellphonedb_corrected.md` | 23 / 20 / 4 |

中文版传 `*_zh.md` 与中文标题；Dynamo 与 CellphoneDB 各自复制了一份 `build.sh`（内容同 `WES/ref/build.sh`）。

---

## 仓库约定 / Repo conventions

- **大文件**：教学文档（HTML/PDF，单文件最大 5.5 MB）有意入库，便于学生一次 clone 拿到全部讲义；其余大文件一律走集群共享目录或 Git LFS。
- **文档插图**：图片放在各模块 `ref/figures/`（WES 的结果导读用 `results/figures/`），markdown 里写相对路径 `./figures/xxx.png`。
- **`.gitignore`**：屏蔽所有测序数据与运行产物（FASTQ/BAM/VCF/GVCF/二进制模型…）；各模块 `results/` 是**逐条列出的例外区**（第 8–15 节），
  例外区内仍硬性拦截原始数据扩展名与 `*.zip`。
- **pre-commit**：自动扫描 API key / 密钥模式（含 base64 感知，避免教学 HTML 内嵌图片误报）、大文件闸门、数据文件闸门；`results/` 例外区已放行。
- **未入库的本地素材**：`homework1/`（教师给的 CCI 作业说明）按约定暂不入库。

---

## 附录：WES 路线图 13 步一览

| # | 步骤 | 工具 | 核心概念 |
|---|---|---|---|
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

---

## 致谢 / Credits

- 路线图与课程内容作者：**Jiahao Ma (马嘉昊)**, **Zexian Zeng (曾泽贤)** — 北京大学前沿交叉学科研究院
  (Center of Life Sciences & Center for Quantitative Biology)。原文面向课堂分发；本仓库仅作教学存档与复现配套，版权归原作者。
- 感谢 Yufeng He 此前的教学路线图工作。
