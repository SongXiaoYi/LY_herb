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
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggpubr)
############################################
Sys.setenv ("VROOM_CONNECTION_SIZE"=99999999)

# load series and platform data from GEO

setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE141910_norm_counts_TPM_GRCh38.p13_NCBI.tsv')
pbmc <- fread('./GSE141910_norm_counts_TPM_GRCh38.p13_NCBI.tsv') %>% as.data.frame()
############################ 获取注释文件
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = as.character(pbmc$GeneID),
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first" # 多映射取第一个
)

gene_type <- mapIds(
  org.Hs.eg.db,          # 人类注释库
  keys = as.character(pbmc$GeneID),
  keytype = "ENTREZID",
  column = "GENETYPE",    # 要转换的目标类型（基因名）
  multiVals = "first"    # 多个匹配时取第一个（避免返回列表）
)

Gene_list <- cbind(entrez_ids, names(entrez_ids), gene_type) %>% as.data.frame()
pbmc$GeneID <- Gene_list$entrez_ids
pbmc$GeneType <- Gene_list$gene_type
pbmc <- na.omit(pbmc)
pbmc <- pbmc[pbmc$GeneType == 'protein-coding',]
Names <- pbmc$GeneID
pbmc <- pbmc[,-c(1, 358)]
pbmc <- as.matrix(pbmc)
rownames(pbmc) <- Names
###############################
pbmc_meta <- fread('./phenoData.csv') %>% as.data.frame()
pbmc_meta <- pbmc_meta[pbmc_meta$etiology %in% c('NF','HCM'),]
pbmc <- pbmc[,colnames(pbmc) %in% pbmc_meta$GEO]
pbmc <- avereps(pbmc)
pbmc <- pbmc[rowMeans(pbmc)>0,]
#########################################
tpm_log <- log2(pbmc + 1)
########################################
library(homologene)
Genelist <- fread('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\WGCNA\\IMConn_rank2.csv') %>% as.data.frame()
Genelist <- Genelist$Gene
TurnGenelist <- mouse2human(Genelist)
##################################################
Data_matrix <- tpm_log[rownames(tpm_log) %in% TurnGenelist$humanGene,]
##############
pbmc <- t(Data_matrix)
Gene_matrix <- as.data.frame(pbmc)
Gene_matrix <- cbind(rownames(Gene_matrix), Gene_matrix$CCNA2) %>% as.data.frame()

colnames(Gene_matrix) <- c('ID','Gene')
Gene_matrix$Group <- pbmc_meta$etiology[match(Gene_matrix$ID, pbmc_meta$GEO)]
Gene_matrix$Gene <- as.numeric(Gene_matrix$Gene)
################################ Liner
newdata <- Gene_matrix
ko_color <- c('#3d4144', '#28aae2')
newdata$Group <- factor(newdata$Group, levels = c('NF','HCM'))
my_comparisons <- list(c('NF','HCM'))
a <- ggboxplot(newdata, x="Group", y="Gene", fill = "Group",
            ylab="Transcripts Per Million (TMP)",add = "jitter", shape = 16,xlab="",bxp.errorbar = TRUE,outlier.shape = NA)
a <- a + theme_bw() + rotate_x_text(51) + scale_fill_manual(values = ko_color)
a <- a + stat_compare_means(comparisons = my_comparisons, symnum.args=list(cutpoints = c(0,0.001,0.01,0.05, 1), symbols = c("***", "**","*", "ns")),
                            label = "p.signif")
a

setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE141910_norm_counts_TPM_GRCh38.p13_NCBI.tsv')
pdf(file = 'CCNA2_compare.pdf', width = 3.4, height = 3.05)
a
dev.off()























