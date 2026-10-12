## =============================================================================
## Lesson 3 — RNA-seq application: Part 4 TIDE analysis (R part)
## Source: RNA-seq-application.html (code blocks 28-33, as written)
## Deviations (pre-registered, additive only):
##   - run via Rscript (tutorial: interactive R) under LC_ALL=C
##   - setwd to this step's directory so png() relative paths land here
##   - load the res_tide.txt produced by our own tidepy run (code block 27)
##     instead of the teacher's pre-baked copy. Tutorial code block 28 comment
##     explicitly documents this as the intended usage when you run tidepy
##     yourself: "If you run tidepy yourself, you can load your results via
##     res_tide <- read.table('/path_to_your_folder/res_tide.txt', ...)"
##   - re-derive metadata here (separate R process) for reordering + Response
##   - diagnostic prints: res_tide head, ROC auc
## =============================================================================
setwd("/gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/05_tide")
options(stringsAsFactors = FALSE)

## ---- re-derive metadata (Part 2 blocks 11,12,14) ----
library(tidyverse)
library(pheatmap)
exprsn <- read.table('/lustre1/share/rnaimmune/data/expression/VanAllen2015_CTLA4_Melanoma_RNASeq', sep='\t', header=TRUE, row.names=1)
exprsn <- t(t(exprsn)/colSums(exprsn)*1e6)
metadata <- read.table('/lustre1/share/rnaimmune/data/follow_up/VanAllen2015_CTLA4_Melanoma_RNASeq', sep='\t', header=TRUE, row.names=1)
metadata <- metadata[colnames(exprsn),]
metadata$CYT <- colMeans(exprsn[c('GZMA','PRF1'),])
metadata <- metadata %>% mutate(CYT_level = ifelse(CYT > median(CYT), "high", "low"))
metadata$CTLA4 <- exprsn['CTLA4',]
metadata$Gender <- factor(metadata$Gender, levels=c(0,1)); levels(metadata$Gender) <- c('Female','Male')
metadata$Response <- factor(metadata$Response, levels=c(0,1)); levels(metadata$Response) <- c('Non-responder', 'Responder')

## ---- code block 28: load TIDE results ----
res_tide <- read.table('res_tide.txt', sep='\t', header=TRUE, row.names=1)
# If you run tidepy yourself, you can load your results via
# res_tide <- read.table('/path_to_your_folder/res_tide.txt', sep='\t', header=TRUE, row.names=1)

colnames(res_tide)[2] <- 'Predicted Responder'
res_tide <- res_tide[rownames(metadata),] # reorder to be consistent with metadata
cat("== res_tide head ==\n"); print(head(res_tide))

# Process the TIDE results for proper visualization
res_tide$`Predicted Responder` <- factor(res_tide$`Predicted Responder`, levels=c('False', 'True'))
levels(res_tide$`Predicted Responder`) <- c('Non-responder', 'Responder')
res_tide$CTL.flag <- factor(res_tide$CTL.flag, levels=c('False', 'True'))
levels(res_tide$CTL.flag) <- c('CTL Low', 'CTL high')
res_tide$No.benefits <- factor(res_tide$No.benefits, levels=c('False', 'True'))

res_tide <- cbind(res_tide, metadata$Response)
colnames(res_tide)[15] <- 'Responder'

## ---- code block 29: TIDE heatmap ----
# annotations for columns
annotation <- res_tide[,c(15, 2, 8, 1)]
res_tide_heatmap <- res_tide[,-c(15, 2, 8, 1)]
annotation <- cbind(annotation, metadata[,c('Age', 'Gender', 'Stage')])
cat("== annotation head ==\n"); print(head(annotation))

# heatmap plot
png('tide_heatmap.png', width=800, height=1000)
pheatmap(t(res_tide_heatmap),
         scale = 'row',
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation
)
dev.off()

## ---- code block 30: TIDE bar (predicted responder) ----
# bar plot
p <- ggplot(res_tide, aes(x = reorder(rownames(res_tide), -TIDE), y = TIDE, fill = `Predicted Responder`))
p <- p + geom_bar(stat = 'identity', width = 0.7)
p <- p + scale_fill_manual(values = c('Responder' = 'red', 'Non-responder' = 'blue'))
p <- p + facet_grid(. ~ CTL.flag, space = 'free_x', scales = 'free_x')
p <- p + labs(x = '', y = 'TIDE Score')
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
png('tide_predicted_responder.png')
print(p)
dev.off()

## ---- code block 31: TIDE bar (true responder) ----
# bar plot
p <- ggplot(res_tide, aes(x = reorder(rownames(res_tide), -TIDE), y = TIDE, fill = Responder))
p <- p + geom_bar(stat = 'identity', width = 0.7)
p <- p + scale_fill_manual(values = c('Responder' = 'red', 'Non-responder' = 'blue'))
p <- p + facet_grid(. ~ CTL.flag, space = 'free_x', scales = 'free_x')
p <- p + labs(x = '', y = 'TIDE Score')
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
png('tide_responder.png')
print(p)
dev.off()

## ---- code block 32: ROC ----
library(pROC)

TIDE_roc <- roc(Responder ~ TIDE, data=res_tide, aur=TRUE)
cat("== TIDE_roc auc ==\n"); print(TIDE_roc$auc)

# ROC plot
p <- ggroc(TIDE_roc, legacy.axes = TRUE)
p <- p + geom_segment(x = 0, xend = 1, y = 0, yend = 1 ,color="darkgrey", linetype=4)
p <- p + theme_bw()
p <- p + annotate("text",x=0.75,y=0.125,label = paste("AUC =", round(TIDE_roc$auc,3)))
png('roc.png')
print(p)
dev.off()

## ---- code block 33: smoothed ROC ----
# smooth the roc curve
TIDE_roc_smooth <- roc(Responder ~ TIDE, data=res_tide, aur=TRUE, smooth=TRUE)
cat("== TIDE_roc_smooth auc ==\n"); print(TIDE_roc_smooth$auc)

# ROC plot
p <- ggroc(TIDE_roc_smooth, legacy.axes = TRUE)
p <- p + geom_segment(x = 0, xend = 1, y = 0, yend = 1 ,color="darkgrey", linetype=4)
p <- p + theme_bw()
p <- p + annotate("text",x=0.75,y=0.125,label = paste("AUC =", round(TIDE_roc_smooth$auc,3)))
png('roc_smooth.png')
print(p)
dev.off()

cat("== Part 4 done. outputs: tide_heatmap/tide_predicted_responder/tide_responder/roc/roc_smooth.png ==\n")
