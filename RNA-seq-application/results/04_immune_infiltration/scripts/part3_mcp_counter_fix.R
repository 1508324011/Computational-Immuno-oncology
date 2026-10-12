# part3_mcp_counter_fix.R — APP-01 矫正版：MCPcounter 本地标志基因集运行
# =============================================================================
# 背景（APP-01，详见 ref/RNA_seq_application_corrected*.md Part 3）：
#   conda 版 MCPcounter 1.1.0 把 genes/probesets 默认参数硬编码成
#   http://raw.githubusercontent.com/... 运行时拉取；登录节点到 GitHub 的
#   http 连接被重置（https 可达），包内又无本地副本 → 按原文原样跑必阻塞。
#   2026-10-12 按用户决定：原样跑记为阻塞（mcp_counter_blocked.txt）。
#
# 矫正（本脚本，2026-10-12 晚补跑，用户指示执行并记录）：
#   用 https 把两个标志文件下到 ref/（MCPcounter_genes.txt /
#   MCPcounter_probesets.txt），read.table 参数逐字复刻 MCPcounter.estimate
#   默认值（仅换数据源为本地文件），作为 genes=/probesets= 显式传入
#   immunedeconv::deconvolute —— 同一份数据，不改包、不改教程逻辑。
#
# 偏差登记：
#   - 本脚本是矫正版补跑，不是教程原样代码（原样代码见
#     part3_immune_infiltration.R 的 tryCatch 段与 mcp_counter_blocked.txt）
#   - 其余部分（读表达矩阵、CPM、deconvolute 调用形式）与 part3 脚本一致
# =============================================================================

setwd('/gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/04_immune_infiltration')
cat('=== APP-01 fix: MCPcounter with local marker files ===\n')
cat('time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n')

suppressMessages(library(immunedeconv))

# --- 读表达矩阵 + CPM（与 part3 脚本一致）---
exprsn <- read.table('/lustre1/share/rnaimmune/data/expression/VanAllen2015_CTLA4_Melanoma_RNASeq',
                     sep='\t', header=TRUE, row.names=1)
exprsn <- t(t(exprsn)/colSums(exprsn)*1e6)
cat('expression:', nrow(exprsn), 'genes x', ncol(exprsn), 'samples (CPM)\n')

# --- 本地标志基因集：read.table 参数逐字复刻包默认（仅换 http 为本地文件）---
# 包默认（args(MCPcounter.estimate)）：
#   probesets = read.table(curl("http://.../probesets.txt"),
#                          sep="\t", stringsAsFactors=FALSE, colClasses="character")
#   genes     = read.table(curl("http://.../genes.txt"), sep="\t",
#                          stringsAsFactors=FALSE, header=TRUE,
#                          colClasses="character", check.names=FALSE)
genes <- read.table('/gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/ref/MCPcounter_genes.txt',
                    sep='\t', stringsAsFactors=FALSE, header=TRUE,
                    colClasses='character', check.names=FALSE)
probesets <- read.table('/gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/ref/MCPcounter_probesets.txt',
                        sep='\t', stringsAsFactors=FALSE, colClasses='character')
cat('local genes:', nrow(genes), 'markers x', ncol(genes), 'cols;',
    nlevels(factor(genes[[2]])), 'cell populations\n')
cat('local probesets:', nrow(probesets), 'rows\n')

# --- deconvolute：与教程调用形式一致，仅显式传 genes=/probesets= ---
res_mcp <- deconvolute(exprsn, 'mcp_counter',
                       genes=genes, probesets=probesets)
cat('res_mcp:', nrow(res_mcp), 'cell types x', ncol(res_mcp)-1, 'samples\n')
print(res_mcp)

# --- 落盘：覆盖原 NULL 占位（NULL 占位内容存档于 mcp_counter_blocked.txt）---
write.table(res_mcp, 'deconv_methods/res_mcp_counter.tsv',
            sep='\t', quote=FALSE, row.names=FALSE)
cat('saved: deconv_methods/res_mcp_counter.tsv\n')
cat('=== APP-01 fix: OK ===\n')
