###################################### Tian He
#####################
library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(patchwork)
library(ggplot2)
library(ggbump)
library(knitr)
library(data.table)
library(dplyr)
library(tidyr)
library(SCENIC)
library(SCopeLoomR)
library(muscat)
library(SingleCellExperiment)
library(RcisTarget)
library(reshape2)
library(geomtextpath)
library(ggpubr)
#library(scGSVA)
library(GSEABase)
library(GSVA)
library(ComplexHeatmap)
library(schard)
library(data.table)
#######################################
###############################
'%!in%' <- Negate('%in%') 
###################################
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE203160')
Data_matrix <- fread('./GSE203160-normalize_tpm_data.csv') %>% as.data.frame()
Data_meta <- fread('./Data_meta.csv') %>% as.data.frame()
rownames(Data_matrix) <- Data_matrix$V1
Data_matrix <- Data_matrix[,-1]
Genelist <- fread('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\WGCNA\\IMConn_rank2.csv') %>% as.data.frame()
Genelist <- Genelist$Gene
TurnGenelist <- mouse2human(Genelist)
###########################
Data_matrix <- Data_matrix[rownames(Data_matrix) %in% TurnGenelist$humanGene,]
##############
pbmc <- t(Data_matrix)
Gene_matrix <- as.data.frame(pbmc)
Gene_matrix <- cbind(rownames(Gene_matrix), Gene_matrix$CCNA2) %>% as.data.frame()

colnames(Gene_matrix) <- c('ID','Gene')
Gene_matrix$Group <- Data_meta$Disease
Gene_matrix$Gene <- as.numeric(Gene_matrix$Gene)
################################ Liner
newdata <- Gene_matrix
ko_color <- c('#3d4144', '#28aae2')
newdata$Group <- factor(newdata$Group, levels = c('Control','Case'))
my_comparisons <- list(c("Control", "Case"))
a <- ggboxplot(newdata, x="Group", y="Gene", fill = "Group",
            ylab="Transcripts Per Million (TMP)",add = "jitter", shape = 16,xlab="",bxp.errorbar = TRUE,outlier.shape = NA)
a <- a + theme_bw() + rotate_x_text(51) + scale_fill_manual(values = ko_color)
a <- a + stat_compare_means(comparisons = my_comparisons, symnum.args=list(cutpoints = c(0,0.001,0.01,0.05, 1), symbols = c("***", "**","*", "ns")),
                            label = "p.signif")
a
##########################
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE203160')
pdf(file = 'Mki67_compare.pdf', width = 3.4, height = 3.05)
a
dev.off()
###################################################
Gene_matrix <- as.data.frame(pbmc)
Gene_matrix <- cbind(rownames(Gene_matrix), Gene_matrix$CCNA2) %>% as.data.frame()

colnames(Gene_matrix) <- c('ID','Gene')
Gene_matrix$Group <- Data_meta$Disease
Gene_matrix$Gene <- as.numeric(Gene_matrix$Gene)
################################ Liner
newdata <- Gene_matrix
ko_color <- c('#3d4144', '#28aae2')
newdata$Group <- factor(newdata$Group, levels = c('Control','Case'))
my_comparisons <- list(c("Control", "Case"))
a <- ggboxplot(newdata, x="Group", y="Gene", fill = "Group",
            ylab="Transcripts Per Million (TMP)",add = "jitter", shape = 16,xlab="",bxp.errorbar = TRUE,outlier.shape = NA)
a <- a + theme_bw() + rotate_x_text(51) + scale_fill_manual(values = ko_color)
a <- a + stat_compare_means(comparisons = my_comparisons, symnum.args=list(cutpoints = c(0,0.001,0.01,0.05, 1), symbols = c("***", "**","*", "ns")),
                            label = "p.signif")
a
##########################
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE42955')
pdf(file = 'CCNA2_compare.pdf', width = 3.4, height = 3.05)
a
dev.off()
