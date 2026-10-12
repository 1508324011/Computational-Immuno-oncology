## =============================================================================
## Lesson 3 — RNA-seq application: Part 2 Survival analysis
## Source: RNA-seq-application.html (code blocks 11-16, as written)
## Deviations (pre-registered, additive only):
##   - run via Rscript (tutorial: interactive R) under LC_ALL=C
##   - setwd to this step's directory so png() relative paths land here
##   - re-derive exprsn/sig here because Part 1 ran in a separate R process
##     (tutorial: single interactive R session shared these objects)
##   - diagnostic prints: metadata head, survfit summary, cox summary
## =============================================================================
setwd("/gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/03_survival_analysis")

## ---- re-derive exprsn & CYT sig (from Part 1 code blocks 6,7,10) ----
exprsn <- read.table('/lustre1/share/rnaimmune/data/expression/VanAllen2015_CTLA4_Melanoma_RNASeq', sep='\t', header=TRUE, row.names=1)
exprsn <- t(t(exprsn)/colSums(exprsn)*1e6)            # CPM (block 7)
geneset <- c('GZMA', 'PRF1')
sig <- colMeans(exprsn[geneset,])                      # CYT (block 10)

## ---- code block 11: load metadata, reorder to exprsn ----
library(survival)
library(survminer)
library(dplyr)

metadata <- read.table('/lustre1/share/rnaimmune/data/follow_up/VanAllen2015_CTLA4_Melanoma_RNASeq', sep='\t', header=TRUE, row.names=1)
metadata <- metadata[colnames(exprsn),] # reorder to be consistent with exprsn
cat("== metadata head ==\n"); print(head(metadata))
cat("== metadata rownames match exprsn colnames? ==\n"); print(all(rownames(metadata)==colnames(exprsn)))

## ---- code block 12: add CYT, CYT_level ----
metadata$CYT <- sig
metadata <- metadata %>% mutate(CYT_level = ifelse(CYT > median(CYT), "high", "low"))
cat("== metadata with CYT/CYT_level head ==\n"); print(head(metadata))
cat("== CYT_level table ==\n"); print(table(metadata$CYT_level))

## ---- code block 13: KM survival plot ----
surv_obj <- Surv(time = metadata$OS, event = metadata$`OS.Event`)
fit <- survfit(surv_obj ~ CYT_level, data = metadata)
cat("== survfit summary ==\n"); print(summary(fit))

# K-M plot
p <- ggsurvplot(fit, data = metadata,
                pval = TRUE,
                conf.int = TRUE,
                risk.table = TRUE,
                ggtheme = theme_minimal())
p <- p + xlab('Time (days)')
png('survival_CYT.png')
print(p)
dev.off()

## ---- code block 14: Cox model ----
metadata$CTLA4 <- exprsn['CTLA4',]
metadata$Gender <- factor(metadata$Gender, levels=c(0,1))
levels(metadata$Gender) <- c('Female', 'Male')
cox_model <- coxph(Surv(OS, `OS.Event`) ~ Age + Gender + Stage + CYT + CYT_level + CTLA4, data = metadata)

## ---- code block 15: cox summary ----
cat("== cox_model summary ==\n"); print(summary(cox_model))

## ---- code block 16: forest plot ----
# forest plot
png('cox.png')
print(ggforest(cox_model, data = metadata,
         main = "Cox Proportional Hazards Model",
         fontsize = 1.2))
dev.off()

cat("== Part 2 done. outputs: survival_CYT.png, cox.png ==\n")
