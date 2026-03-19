###################################### 
####### WGCNA-PRAD
library(dplyr)
library(stringr)
library(patchwork)
library(knitr)
library(msigdbr) 
library(fgsea)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(limma)
library(clusterProfiler)
#library(DESeq2)
library(ggplotify)
library(cowplot)
library(edgeR)
library(tinyarray)
library(gridExtra)
library(Seurat)
library(WGCNA)
library(reshape2)
library('org.Mm.eg.db')
library('clusterProfiler')
library(ggpubr)
######################################
#####################################################
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\WGCNA')
pbmc <- read.csv('./TPM.gene.csv')
############################################## 基因转换
library(biomaRt)
GeneENSUM <- pbmc$geneID
ENTREZID <- bitr(GeneENSUM, fromType = "ENSEMBL", toType=c("SYMBOL","ENTREZID","GENETYPE"),OrgDb = org.Mm.eg.db)
ENTREZID <- ENTREZID[ENTREZID$GENETYPE == 'protein-coding',]
pbmc <- pbmc[pbmc$geneID %in% ENTREZID$ENSEMBL,]
pbmc$Match_Symbol <- ENTREZID$SYMBOL[match(pbmc$geneID, ENTREZID$ENSEMBL)]
Names <- pbmc$Match_Symbol
pbmc_matrix <- as.matrix(pbmc[,-1])
pbmc_matrix <- as.matrix(pbmc_matrix[,-34])
pbmc_matrix <- apply(pbmc_matrix,2,as.numeric)
rownames(pbmc_matrix) <- Names
pbmc_matrix <- avereps(pbmc_matrix)
pbmc_matrix <- pbmc_matrix[rowMeans(pbmc_matrix)>0,]
pbmc_matrix <- normalizeBetweenArrays(pbmc_matrix)

###############################################
pbmc_meta <- read.csv('./WGCNA_meta.csv')
rownames(pbmc_meta) <- pbmc_meta$Sample
########################################################
pbmc <- as.data.frame(pbmc_matrix)
pbmc <- t(pbmc)

Tcta_matrix <- as.data.frame(pbmc)
Tcta_matrix <- cbind(rownames(Tcta_matrix), Tcta_matrix$Tcta) %>% as.data.frame()

colnames(Tcta_matrix) <- c('ID','Tcta')
Tcta_matrix$Group <- pbmc_meta$subtype
Tcta_matrix$Tcta <- as.numeric(Tcta_matrix$Tcta)
################################ Liner
newdata <- Tcta_matrix
ko_color <- c('#3d4144', '#28aae2', '#1fb573', '#800040')
newdata$Group <- factor(newdata$Group, levels = c('SHAM','TAC','ART','ASIV80'))
my_comparisons <- list(c("ASIV80", "SHAM"),c("ASIV80", "TAC"),c("ASIV80", "ART"))
a <- ggboxplot(newdata, x="Group", y="Tcta", fill = "Group",
            ylab="Transcripts Per Million (TMP)",add = "jitter", shape = 16,xlab="",bxp.errorbar = TRUE,outlier.shape = NA)
a <- a + theme_bw() + rotate_x_text(51) + scale_fill_manual(values = ko_color)
a <- a + stat_compare_means(comparisons = my_comparisons, symnum.args=list(cutpoints = c(0,0.001,0.01,0.05, 1), symbols = c("***", "**","*", "ns")),
                            label = "p.signif")
a
#a <- a + stat_compare_means() + NoLegend()

setwd('G:\\WuXiang\\Experiment\\LOG')
pdf(file = 'Tcta_compare.pdf', width = 2.7, height = 2.95)
a
dev.off()



