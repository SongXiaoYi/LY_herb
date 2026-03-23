# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0, DESeq2 1.38.3
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
#library(umap)
library(dplyr) # 函数filter
library(tidyr) 
library(limma)
library(sva)
library(preprocessCore)
library(biomaRt)
# Download all cancer data
###############################
library(TCGAbiolinks)
library(MultiAssayExperiment)
#library(maftools)
library(ComplexHeatmap)
library(data.table)
library(hugene10sttranscriptcluster.db)

Sys.setenv ("VROOM_CONNECTION_SIZE"=99999999)

# load series and platform data from GEO

setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE57338')
################################# GPL11532
gset <- getGEO("GSE57338", GSEMatrix =TRUE, AnnotGPL = T,getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL11532", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

GPL11532_matrix <- exprs(gset)

############################ 获取注释文件
GPL11532 <-getGEO('GPL11532',destdir =".")# 读取对象
GPL11532_anno <- Table(GPL11532)
test_genes <- as.character(GPL11532_anno$ID)
############################ 转换ENTREZ_GENE_ID
gene_info <- select(hugene10sttranscriptcluster.db, test_genes, c("SYMBOL", "ENTREZID", "GENENAME","GENETYPE"))
###################### 转换
#colnames(GPL11532_anno)[2] <- 'MatchID'
GPL11532_anno$MatchGene <- gene_info$SYMBOL[match(GPL11532_anno$ID, gene_info$PROBEID)]
GPL11532_anno$MatchType <- gene_info$GENETYPE[match(GPL11532_anno$ID, gene_info$PROBEID)]
GPL11532_anno_select <- GPL11532_anno[GPL11532_anno$MatchType == 'protein-coding',]
GPL11532_anno_select <- na.omit(GPL11532_anno_select)
GPL11532_anno_select <- GPL11532_anno_select[GPL11532_anno_select$ID %in% rownames(GPL11532_matrix),]
rownames(GPL11532_matrix) <- GPL11532_anno$ID
GPL11532_matrix <- as.data.frame(GPL11532_matrix)
GPL11532_matrix$ID <- rownames(GPL11532_matrix)
GPL11532_anno_select$ID <- as.character(GPL11532_anno_select$ID)
GPL11532_matrix$MatchGene <- GPL11532_anno_select$MatchGene[match(GPL11532_matrix$ID, GPL11532_anno_select$ID)]
####################################################
GPL11532_matrix <- na.omit(GPL11532_matrix)
Names <- GPL11532_matrix$MatchGene
GPL11532_matrix <- GPL11532_matrix[,-c(30,31)]
GPL11532_matrix <- as.matrix(GPL11532_matrix)
rownames(GPL11532_matrix) <- Names
GPL11532_matrix <- avereps(GPL11532_matrix)
GPL11532_matrix <- GPL11532_matrix[rowMeans(GPL11532_matrix)>0,]

# log2 transform
qx <- as.numeric(quantile(GPL11532_matrix, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { GPL11532_matrix[which(GPL11532_matrix <= 0)] <- NaN
  GPL11532_matrix <- log2(GPL11532_matrix) }

#取交集
data <- GPL11532_matrix
################################# 筛选出膜性肾病的样本
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE42955')
#meta <- read.csv('./GSE108109_meta.csv')
#final_data <- data[,meta$ID]

write.csv(data,file = "GSE42955-normalize_data.csv")
