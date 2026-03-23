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
library(mygene)
library(hugene10sttranscriptcluster.db)
library(org.Hs.eg.db)
library(tximport)
library(readr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

Sys.setenv ("VROOM_CONNECTION_SIZE"=99999999)

# load series and platform data from GEO

setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE203160')
################################# 
GSE203160_matrix <- fread('GSE203160-normalize_data.csv') %>% as.data.frame()

ensg_ids <- GSE203160_matrix$Gene
ensg_ids <- gsub("\\.\\d+$", "", ensg_ids)

gene_symbol <- mapIds(
  org.Hs.eg.db,          # 人类注释库
  keys = ensg_ids,        # 输入的ID
  keytype = "ENSEMBL",   # 输入ID的类型（Ensembl ID）
  column = "SYMBOL",     # 要转换的目标类型（基因名）
  multiVals = "first"    # 多个匹配时取第一个（避免返回列表）
)

gene_type <- mapIds(
  org.Hs.eg.db,          # 人类注释库
  keys = ensg_ids,        # 输入的ID
  keytype = "ENSEMBL",   # 输入ID的类型（Ensembl ID）
  column = "GENETYPE",     # 要转换的目标类型（基因名）
  multiVals = "first"    # 多个匹配时取第一个（避免返回列表）
)

Gene_list <- cbind(gene_symbol, names(gene_symbol),gene_type) %>% as.data.frame()
colnames(Gene_list) <- c('gene_symbol','ENSEMBL','genetype')
GSE203160_matrix$Gene <- gsub("\\.\\d+$", "", GSE203160_matrix$Gene)
GSE203160_matrix$MatchGene <- Gene_list$gene_symbol[match(GSE203160_matrix$Gene, Gene_list$ENSEMBL)]
GSE203160_matrix$GeneType <- Gene_list$genetype[match(GSE203160_matrix$Gene, Gene_list$ENSEMBL)]
GSE203160_matrix <- na.omit(GSE203160_matrix)
GSE203160_matrix <- GSE203160_matrix[GSE203160_matrix$GeneType == 'protein-coding',]
Names <- GSE203160_matrix$MatchGene
GSE203160_matrix <- GSE203160_matrix[,-c(1,17,18)]
GSE203160_matrix <- as.matrix(GSE203160_matrix)
rownames(GSE203160_matrix) <- Names

GSE203160_matrix <- avereps(GSE203160_matrix)
GSE203160_matrix <- GSE203160_matrix[rowMeans(GSE203160_matrix)>0,]
########################################
symbols <- rownames(GSE203160_matrix)

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = symbols,
  keytype = "SYMBOL",
  column = "ENTREZID"
)

gene_len <- data.frame(
  symbol = symbols,
  entrez = as.character(entrez_ids)
)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
exons_by_gene <- exonsBy(txdb, by="gene")
len <- sum(width(reduce(exons_by_gene)))

# 匹配 Symbol ↔ 长度
gene_len$length <- len[as.character(gene_len$entrez)]
gene_len <- gene_len[!is.na(gene_len$length) & gene_len$length > 0, ]

# ====================== 4. 对齐count ======================
common <- intersect(rownames(GSE203160_matrix), gene_len$symbol)
GSE203160_matrix_final <- GSE203160_matrix[common, ]
length_final <- gene_len$length[match(common, gene_len$symbol)]

# ====================== 5. 计算 TPM ======================
rate <- GSE203160_matrix_final / length_final
col_sums <- colSums(rate)
tpm <- t(t(rate) * 1e6 / col_sums)
tpm <- round(tpm,4)
tpm_log <- log2(tpm + 1)
# ====================== 6. 输出 ======================

#取交集
data <- tpm_log
################################# 筛选出膜性肾病的样本
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE203160')
#meta <- read.csv('./GSE108109_meta.csv')
#final_data <- data[,meta$ID]

write.csv(data,file = "GSE203160-normalize_tpm_data.csv")
