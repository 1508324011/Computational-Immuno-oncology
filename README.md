# Computational-Immuno-oncology

可复现的课堂内容仓库 / Reproducible course materials for computational immuno-oncology teaching.

> 当前内容：**Whole-Exome Sequencing (WES) 数据分析完整路线图**。
> 后续将逐步加入免疫肿瘤学其他模块（TCR/BCR 分析、新抗原预测、肿瘤免疫微环境等）。

---

## 目录结构 / Layout

```text
.
└── WES/
    └── ref/                                  # 参考文档 (reference material)
        ├── WES_roadmap_to_students.html      # 主路线图 (pandoc 导出, 图文完整)
        ├── WES_roadmap_to_students.pdf       # 同一文档的 PDF 版
        ├── WES_roadmap_english_embedded.html # 英文版 (图片内嵌)
        ├── WES_roadmap_corrected.html        # 矫正增强版：两轮实跑对照 + 逐坑修复说明
        └── WES_roadmap_corrected.pdf         # 同一文档的 PDF 版
```

> **矫正增强版（WES_roadmap_corrected）**：把学生版路线图在集群上**原样跑两轮**（0.5× 降采样 demo 轮 + 全量轮）后，对照实测结果写成的修复与讲解版——每个步骤含原理/目的/实跑结果/坑与矫正/结果解读/调优扩展，附录带完整问题清单（12 条，含证据）与一页矫正版命令速查。已知问题速览：VQSR MergeVcfs 每条变异重复两遍（含老师参考产物）、ANNOVAR 富协议库名不存在、CNVkit amplicon 模式用于杂交捕获 panel、无靶区深度度量（实测全量 OC 仅 54.6× 可用）。

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

## 致谢 / Credits

- 路线图作者：**Jiahao Ma (马嘉昊)**, **Zexian Zeng (曾泽贤)** — 北京大学前沿交叉学科研究院 (Center of Life Sciences & Center for Quantitative Biology)
  原文面向课堂分发；本仓库仅作教学存档与复现配套，版权归原作者。
- 感谢 Yufeng He 此前的教学路线图工作。

## 仓库约定 / Repo conventions

- 大文件（HTML/PDF 教学文档 ≈16MB）有意直接入库，便于学生一次 clone 拿到全部讲义；其余大文件一律走 [Git LFS] 或集群共享目录。
- `.gitignore` 屏蔽所有测序数据与运行产物 (FASTQ/BAM/VCF/…) —— **任何原始数据、患者数据、凭据不入库**。
- 提交时 pre-commit 钩子会自动扫描 API key / 密钥模式（含 base64 感知，避免教学 HTML 中的嵌入图片误报）。
