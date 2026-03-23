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

setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO')
################################# GPL6244
gset <- getGEO("GSE42955", GSEMatrix =TRUE, AnnotGPL = T,getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

GPL6244_matrix <- exprs(gset)

############################ 获取注释文件
GPL6244 <-getGEO('GPL6244',destdir =".")# 读取对象
GPL6244_anno <- Table(GPL6244)
test_genes <- as.character(GPL6244_anno$ID)
############################ 转换ENTREZ_GENE_ID
gene_info <- select(hugene10sttranscriptcluster.db, test_genes, c("SYMBOL", "ENTREZID", "GENENAME","GENETYPE"))
###################### 转换
#colnames(GPL6244_anno)[2] <- 'MatchID'
GPL6244_anno$MatchGene <- gene_info$SYMBOL[match(GPL6244_anno$ID, gene_info$PROBEID)]
GPL6244_anno$MatchType <- gene_info$GENETYPE[match(GPL6244_anno$ID, gene_info$PROBEID)]
GPL6244_anno_select <- GPL6244_anno[GPL6244_anno$MatchType == 'protein-coding',]
GPL6244_anno_select <- na.omit(GPL6244_anno_select)
GPL6244_anno_select <- GPL6244_anno_select[GPL6244_anno_select$ID %in% rownames(GPL6244_matrix),]
rownames(GPL6244_matrix) <- GPL6244_anno$ID
GPL6244_matrix <- as.data.frame(GPL6244_matrix)
GPL6244_matrix$ID <- rownames(GPL6244_matrix)
GPL6244_anno_select$ID <- as.character(GPL6244_anno_select$ID)
GPL6244_matrix$MatchGene <- GPL6244_anno_select$MatchGene[match(GPL6244_matrix$ID, GPL6244_anno_select$ID)]
####################################################
GPL6244_matrix <- na.omit(GPL6244_matrix)
Names <- GPL6244_matrix$MatchGene
GPL6244_matrix <- GPL6244_matrix[,-c(30,31)]
GPL6244_matrix <- as.matrix(GPL6244_matrix)
rownames(GPL6244_matrix) <- Names
GPL6244_matrix <- avereps(GPL6244_matrix)
GPL6244_matrix <- GPL6244_matrix[rowMeans(GPL6244_matrix)>0,]

# log2 transform
qx <- as.numeric(quantile(GPL6244_matrix, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { GPL6244_matrix[which(GPL6244_matrix <= 0)] <- NaN
  GPL6244_matrix <- log2(GPL6244_matrix) }

#取交集
data <- GPL6244_matrix
################################# 筛选出膜性肾病的样本
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE42955')
#meta <- read.csv('./GSE108109_meta.csv')
#final_data <- data[,meta$ID]

write.csv(data,file = "GSE42955-normalize_data.csv")
