#加载包
#####################################
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(harmony)
library(patchwork)
library(data.table)
library(magrittr)
library(ggplot2)
#####################################
setwd('I:\\LY\\geo\\GSE247468_RAW')
####################
###################################################  Control_1-case
Control_1.data <- Read10X(data.dir = "./Control_1_filtered_feature_bc_matrix")  #改目录
dense.size <- object.size(x = as.matrix(x = Control_1.data))
sparse.size <- object.size(x = Control_1.data)
Control_1 <- CreateAssayObject(counts = Control_1.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Control_1 <- CreateSeuratObject(Control_1, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Control_1@meta.data$Status<- "Control"
Control_1@meta.data$ID<- "Control_1"
##################################### 线粒体基因
Control_1[["percent.mt"]] <- PercentageFeatureSet(object = Control_1, pattern = "^MT-")     #注意MT
##################################### 红细胞基因
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Control_1@assays$RNA)) 
HB.genes <- rownames(Control_1@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Control_1[["percent.HB"]]<-PercentageFeatureSet(Control_1, features=HB.genes) 
##################################### 核糖体基因
#C<-GetAssayData(object = Control_1, slot = "counts")
#rb.genes <- rownames(Control_1)[grep("^RP[SL]",rownames(Control_1))]
#percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
#Control_1 <- AddMetaData(Control_1, percent.ribo, col.name = "percent.ribo")
#col.num <- length(levels(Control_1@active.ident))
##################################### 可视化
#violin <- VlnPlot(Control_1,
#                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","percent.ribo"), 
#                  cols =rainbow(col.num), 
#                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
#                  ncol = 5) + 
#  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
#####################################
Control_1 <- subset(x = Control_1, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 30 & percent.HB == 0)     #这里取了线粒体10%
Control_1 <- NormalizeData(object = Control_1, normalization.method = "LogNormalize", scale.factor = 10000)
Control_1 <- FindVariableFeatures(object = Control_1, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(Control_1)
Control_1 <- ScaleData(Control_1, features = all.genes)
Control_1 <- RunPCA(object= Control_1,npcs = 20,pc.genes=VariableFeatures(object = Control_1))  #这里取前20个PCs
pcSelect=20
Control_1 <- FindNeighbors(object = Control_1, dims = 1:pcSelect) 
Control_1 <- FindClusters(object = Control_1, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
Control_1 <- RunUMAP(Control_1, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(Control_1, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Control_1@meta.data$RNA_snn_res.0.4
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Control_1@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Control_1 <- doubletFinder_v3(Control_1, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Control_1 <- doubletFinder_v3(Control_1, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Control_1@meta.data$Doublet <- Control_1@meta.data$DF.classifications_0.25_0.005_330
table(Control_1@meta.data$DF.classifications_0.25_0.005_330)
###################################### 删除所有红细胞基因及线粒体基因
'%!in%' <- Negate('%in%')
MT.genes <- rownames(Control_1)[grep("^MT-",rownames(Control_1))]
Delete_genes <- c(MT.genes,HB.genes)
Control_1 <- Control_1[which(rownames(Control_1) %!in% Delete_genes),]
#保存
setwd('I:\\LY\\geo\\GSE247468_RAW')
saveRDS(Control_1, file = "./Control_1.rds")
#删除占内存的
rm(list = ls())
gc()


###################################################  Control_2-case
Control_2.data <- Read10X(data.dir = "./Control_2_filtered_feature_bc_matrix")  #改目录
dense.size <- object.size(x = as.matrix(x = Control_2.data))
sparse.size <- object.size(x = Control_2.data)
Control_2 <- CreateAssayObject(counts = Control_2.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Control_2 <- CreateSeuratObject(Control_2, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Control_2@meta.data$Status<- "Control"
Control_2@meta.data$ID<- "Control_2"
##################################### 线粒体基因
Control_2[["percent.mt"]] <- PercentageFeatureSet(object = Control_2, pattern = "^MT-")     #注意MT
##################################### 红细胞基因
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Control_2@assays$RNA)) 
HB.genes <- rownames(Control_2@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Control_2[["percent.HB"]]<-PercentageFeatureSet(Control_2, features=HB.genes) 
##################################### 核糖体基因
#C<-GetAssayData(object = Control_2, slot = "counts")
#rb.genes <- rownames(Control_2)[grep("^RP[SL]",rownames(Control_2))]
#percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
#Control_2 <- AddMetaData(Control_2, percent.ribo, col.name = "percent.ribo")
#col.num <- length(levels(Control_2@active.ident))
##################################### 可视化
#violin <- VlnPlot(Control_2,
#                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","percent.ribo"), 
#                  cols =rainbow(col.num), 
#                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
#                  ncol = 5) + 
#  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
#####################################
Control_2 <- subset(x = Control_2, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 30 & percent.HB == 0)     #这里取了线粒体10%
Control_2 <- NormalizeData(object = Control_2, normalization.method = "LogNormalize", scale.factor = 10000)
Control_2 <- FindVariableFeatures(object = Control_2, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(Control_2)
Control_2 <- ScaleData(Control_2, features = all.genes)
Control_2 <- RunPCA(object= Control_2,npcs = 20,pc.genes=VariableFeatures(object = Control_2))  #这里取前20个PCs
pcSelect=20
Control_2 <- FindNeighbors(object = Control_2, dims = 1:pcSelect) 
Control_2 <- FindClusters(object = Control_2, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
Control_2 <- RunUMAP(Control_2, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(Control_2, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Control_2@meta.data$RNA_snn_res.0.4
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Control_2@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Control_2 <- doubletFinder_v3(Control_2, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Control_2 <- doubletFinder_v3(Control_2, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Control_2@meta.data$Doublet <- Control_2@meta.data$DF.classifications_0.25_0.3_432
table(Control_2@meta.data$DF.classifications_0.25_0.3_432)
###################################### 删除所有红细胞基因及线粒体基因
'%!in%' <- Negate('%in%')
MT.genes <- rownames(Control_2)[grep("^MT-",rownames(Control_2))]
Delete_genes <- c(MT.genes,HB.genes)
Control_2 <- Control_2[which(rownames(Control_2) %!in% Delete_genes),]
#保存
setwd('I:\\LY\\geo\\GSE247468_RAW')
saveRDS(Control_2, file = "./Control_2.rds")
#删除占内存的
rm(list = ls())
gc()


###################################################  Heart_failure_1-case
Heart_failure_1.data <- Read10X(data.dir = "./Heart_failure_1_filtered_feature_bc_matrix")  #改目录
dense.size <- object.size(x = as.matrix(x = Heart_failure_1.data))
sparse.size <- object.size(x = Heart_failure_1.data)
Heart_failure_1 <- CreateAssayObject(counts = Heart_failure_1.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Heart_failure_1 <- CreateSeuratObject(Heart_failure_1, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Heart_failure_1@meta.data$Status<- "Case"
Heart_failure_1@meta.data$ID<- "Heart_failure_1"
##################################### 线粒体基因
Heart_failure_1[["percent.mt"]] <- PercentageFeatureSet(object = Heart_failure_1, pattern = "^MT-")     #注意MT
##################################### 红细胞基因
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Heart_failure_1@assays$RNA)) 
HB.genes <- rownames(Heart_failure_1@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Heart_failure_1[["percent.HB"]]<-PercentageFeatureSet(Heart_failure_1, features=HB.genes) 
##################################### 核糖体基因
#C<-GetAssayData(object = Heart_failure_1, slot = "counts")
#rb.genes <- rownames(Heart_failure_1)[grep("^RP[SL]",rownames(Heart_failure_1))]
#percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
#Heart_failure_1 <- AddMetaData(Heart_failure_1, percent.ribo, col.name = "percent.ribo")
#col.num <- length(levels(Heart_failure_1@active.ident))
##################################### 可视化
#violin <- VlnPlot(Heart_failure_1,
#                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","percent.ribo"), 
#                  cols =rainbow(col.num), 
#                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
#                  ncol = 5) + 
#  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
#####################################
Heart_failure_1 <- subset(x = Heart_failure_1, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 30 & percent.HB == 0)     #这里取了线粒体10%
Heart_failure_1 <- NormalizeData(object = Heart_failure_1, normalization.method = "LogNormalize", scale.factor = 10000)
Heart_failure_1 <- FindVariableFeatures(object = Heart_failure_1, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(Heart_failure_1)
Heart_failure_1 <- ScaleData(Heart_failure_1, features = all.genes)
Heart_failure_1 <- RunPCA(object= Heart_failure_1,npcs = 20,pc.genes=VariableFeatures(object = Heart_failure_1))  #这里取前20个PCs
pcSelect=20
Heart_failure_1 <- FindNeighbors(object = Heart_failure_1, dims = 1:pcSelect) 
Heart_failure_1 <- FindClusters(object = Heart_failure_1, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
Heart_failure_1 <- RunUMAP(Heart_failure_1, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(Heart_failure_1, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Heart_failure_1@meta.data$RNA_snn_res.0.4
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Heart_failure_1@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Heart_failure_1 <- doubletFinder_v3(Heart_failure_1, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Heart_failure_1 <- doubletFinder_v3(Heart_failure_1, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Heart_failure_1@meta.data$Doublet <- Heart_failure_1@meta.data$DF.classifications_0.25_0.28_373
table(Heart_failure_1@meta.data$DF.classifications_0.25_0.28_373)
###################################### 删除所有红细胞基因及线粒体基因
'%!in%' <- Negate('%in%')
MT.genes <- rownames(Heart_failure_1)[grep("^MT-",rownames(Heart_failure_1))]
Delete_genes <- c(MT.genes,HB.genes)
Heart_failure_1 <- Heart_failure_1[which(rownames(Heart_failure_1) %!in% Delete_genes),]
#保存
setwd('I:\\LY\\geo\\GSE247468_RAW')
saveRDS(Heart_failure_1, file = "./Heart_failure_1.rds")
#删除占内存的
rm(list = ls())
gc()

###################################################  Heart_failure_2-case
Heart_failure_2.data <- Read10X(data.dir = "./Heart_failure_2_filtered_feature_bc_matrix")  #改目录
dense.size <- object.size(x = as.matrix(x = Heart_failure_2.data))
sparse.size <- object.size(x = Heart_failure_2.data)
Heart_failure_2 <- CreateAssayObject(counts = Heart_failure_2.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Heart_failure_2 <- CreateSeuratObject(Heart_failure_2, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Heart_failure_2@meta.data$Status<- "Case"
Heart_failure_2@meta.data$ID<- "Heart_failure_2"
##################################### 线粒体基因
Heart_failure_2[["percent.mt"]] <- PercentageFeatureSet(object = Heart_failure_2, pattern = "^MT-")     #注意MT
##################################### 红细胞基因
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Heart_failure_2@assays$RNA)) 
HB.genes <- rownames(Heart_failure_2@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Heart_failure_2[["percent.HB"]]<-PercentageFeatureSet(Heart_failure_2, features=HB.genes) 
##################################### 核糖体基因
#C<-GetAssayData(object = Heart_failure_2, slot = "counts")
#rb.genes <- rownames(Heart_failure_2)[grep("^RP[SL]",rownames(Heart_failure_2))]
#percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
#Heart_failure_2 <- AddMetaData(Heart_failure_2, percent.ribo, col.name = "percent.ribo")
#col.num <- length(levels(Heart_failure_2@active.ident))
##################################### 可视化
#violin <- VlnPlot(Heart_failure_2,
#                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","percent.ribo"), 
#                  cols =rainbow(col.num), 
#                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
#                  ncol = 5) + 
#  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
#####################################
Heart_failure_2 <- subset(x = Heart_failure_2, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 30 & percent.HB == 0)     #这里取了线粒体10%
Heart_failure_2 <- NormalizeData(object = Heart_failure_2, normalization.method = "LogNormalize", scale.factor = 10000)
Heart_failure_2 <- FindVariableFeatures(object = Heart_failure_2, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(Heart_failure_2)
Heart_failure_2 <- ScaleData(Heart_failure_2, features = all.genes)
Heart_failure_2 <- RunPCA(object= Heart_failure_2,npcs = 20,pc.genes=VariableFeatures(object = Heart_failure_2))  #这里取前20个PCs
pcSelect=20
Heart_failure_2 <- FindNeighbors(object = Heart_failure_2, dims = 1:pcSelect) 
Heart_failure_2 <- FindClusters(object = Heart_failure_2, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
Heart_failure_2 <- RunUMAP(Heart_failure_2, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(Heart_failure_2, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Heart_failure_2@meta.data$RNA_snn_res.0.4
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Heart_failure_2@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Heart_failure_2 <- doubletFinder_v3(Heart_failure_2, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Heart_failure_2 <- doubletFinder_v3(Heart_failure_2, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Heart_failure_2@meta.data$Doublet <- Heart_failure_2@meta.data$DF.classifications_0.25_0.09_368
table(Heart_failure_2@meta.data$DF.classifications_0.25_0.09_368)
###################################### 删除所有红细胞基因及线粒体基因
'%!in%' <- Negate('%in%')
MT.genes <- rownames(Heart_failure_2)[grep("^MT-",rownames(Heart_failure_2))]
Delete_genes <- c(MT.genes,HB.genes)
Heart_failure_2 <- Heart_failure_2[which(rownames(Heart_failure_2) %!in% Delete_genes),]
#保存
setwd('I:\\LY\\geo\\GSE247468_RAW')
saveRDS(Heart_failure_2, file = "./Heart_failure_2.rds")
#删除占内存的
rm(list = ls())
gc()

###################################################  Heart_failure_3-case
Heart_failure_3.data <- Read10X(data.dir = "./Heart_failure_3_filtered_feature_bc_matrix")  #改目录
dense.size <- object.size(x = as.matrix(x = Heart_failure_3.data))
sparse.size <- object.size(x = Heart_failure_3.data)
Heart_failure_3 <- CreateAssayObject(counts = Heart_failure_3.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Heart_failure_3 <- CreateSeuratObject(Heart_failure_3, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Heart_failure_3@meta.data$Status<- "Case"
Heart_failure_3@meta.data$ID<- "Heart_failure_3"
##################################### 线粒体基因
Heart_failure_3[["percent.mt"]] <- PercentageFeatureSet(object = Heart_failure_3, pattern = "^MT-")     #注意MT
##################################### 红细胞基因
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Heart_failure_3@assays$RNA)) 
HB.genes <- rownames(Heart_failure_3@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Heart_failure_3[["percent.HB"]]<-PercentageFeatureSet(Heart_failure_3, features=HB.genes) 
##################################### 核糖体基因
#C<-GetAssayData(object = Heart_failure_3, slot = "counts")
#rb.genes <- rownames(Heart_failure_3)[grep("^RP[SL]",rownames(Heart_failure_3))]
#percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
#Heart_failure_3 <- AddMetaData(Heart_failure_3, percent.ribo, col.name = "percent.ribo")
#col.num <- length(levels(Heart_failure_3@active.ident))
##################################### 可视化
#violin <- VlnPlot(Heart_failure_3,
#                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","percent.ribo"), 
#                  cols =rainbow(col.num), 
#                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
#                  ncol = 5) + 
#  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
#####################################
Heart_failure_3 <- subset(x = Heart_failure_3, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 30 & percent.HB == 0)     #这里取了线粒体10%
Heart_failure_3 <- NormalizeData(object = Heart_failure_3, normalization.method = "LogNormalize", scale.factor = 10000)
Heart_failure_3 <- FindVariableFeatures(object = Heart_failure_3, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(Heart_failure_3)
Heart_failure_3 <- ScaleData(Heart_failure_3, features = all.genes)
Heart_failure_3 <- RunPCA(object= Heart_failure_3,npcs = 20,pc.genes=VariableFeatures(object = Heart_failure_3))  #这里取前20个PCs
pcSelect=20
Heart_failure_3 <- FindNeighbors(object = Heart_failure_3, dims = 1:pcSelect) 
Heart_failure_3 <- FindClusters(object = Heart_failure_3, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
Heart_failure_3 <- RunUMAP(Heart_failure_3, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(Heart_failure_3, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Heart_failure_3@meta.data$RNA_snn_res.0.4
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Heart_failure_3@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Heart_failure_3 <- doubletFinder_v3(Heart_failure_3, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Heart_failure_3 <- doubletFinder_v3(Heart_failure_3, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Heart_failure_3@meta.data$Doublet <- Heart_failure_3@meta.data$DF.classifications_0.25_0.25_432
table(Heart_failure_3@meta.data$DF.classifications_0.25_0.25_432)
###################################### 删除所有红细胞基因及线粒体基因
'%!in%' <- Negate('%in%')
MT.genes <- rownames(Heart_failure_3)[grep("^MT-",rownames(Heart_failure_3))]
Delete_genes <- c(MT.genes,HB.genes)
Heart_failure_3 <- Heart_failure_3[which(rownames(Heart_failure_3) %!in% Delete_genes),]
#保存
setwd('I:\\LY\\geo\\GSE247468_RAW')
saveRDS(Heart_failure_3, file = "./Heart_failure_3.rds")
#删除占内存的
rm(list = ls())
gc()


###################################################  Heart_failure_4-case
Heart_failure_4.data <- Read10X(data.dir = "./Heart_failure_4_filtered_feature_bc_matrix")  #改目录
dense.size <- object.size(x = as.matrix(x = Heart_failure_4.data))
sparse.size <- object.size(x = Heart_failure_4.data)
Heart_failure_4 <- CreateAssayObject(counts = Heart_failure_4.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Heart_failure_4 <- CreateSeuratObject(Heart_failure_4, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_ESCC")
Heart_failure_4@meta.data$Status<- "Case"
Heart_failure_4@meta.data$ID<- "Heart_failure_4"
##################################### 线粒体基因
Heart_failure_4[["percent.mt"]] <- PercentageFeatureSet(object = Heart_failure_4, pattern = "^MT-")     #注意MT
##################################### 红细胞基因
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Heart_failure_4@assays$RNA)) 
HB.genes <- rownames(Heart_failure_4@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Heart_failure_4[["percent.HB"]]<-PercentageFeatureSet(Heart_failure_4, features=HB.genes) 
##################################### 核糖体基因
#C<-GetAssayData(object = Heart_failure_4, slot = "counts")
#rb.genes <- rownames(Heart_failure_4)[grep("^RP[SL]",rownames(Heart_failure_4))]
#percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
#Heart_failure_4 <- AddMetaData(Heart_failure_4, percent.ribo, col.name = "percent.ribo")
#col.num <- length(levels(Heart_failure_4@active.ident))
##################################### 可视化
#violin <- VlnPlot(Heart_failure_4,
#                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","percent.ribo"), 
#                  cols =rainbow(col.num), 
#                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
#                  ncol = 5) + 
#  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
#####################################
Heart_failure_4 <- subset(x = Heart_failure_4, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 30 & percent.HB == 0)     #这里取了线粒体10%
Heart_failure_4 <- NormalizeData(object = Heart_failure_4, normalization.method = "LogNormalize", scale.factor = 10000)
Heart_failure_4 <- FindVariableFeatures(object = Heart_failure_4, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(Heart_failure_4)
Heart_failure_4 <- ScaleData(Heart_failure_4, features = all.genes)
Heart_failure_4 <- RunPCA(object= Heart_failure_4,npcs = 20,pc.genes=VariableFeatures(object = Heart_failure_4))  #这里取前20个PCs
pcSelect=20
Heart_failure_4 <- FindNeighbors(object = Heart_failure_4, dims = 1:pcSelect) 
Heart_failure_4 <- FindClusters(object = Heart_failure_4, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
Heart_failure_4 <- RunUMAP(Heart_failure_4, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(Heart_failure_4, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Heart_failure_4@meta.data$RNA_snn_res.0.4
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Heart_failure_4@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Heart_failure_4 <- doubletFinder_v3(Heart_failure_4, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Heart_failure_4 <- doubletFinder_v3(Heart_failure_4, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Heart_failure_4@meta.data$Doublet <- Heart_failure_4@meta.data$DF.classifications_0.25_0.21_307
table(Heart_failure_4@meta.data$DF.classifications_0.25_0.21_307)
###################################### 删除所有红细胞基因及线粒体基因
'%!in%' <- Negate('%in%')
MT.genes <- rownames(Heart_failure_4)[grep("^MT-",rownames(Heart_failure_4))]
Delete_genes <- c(MT.genes,HB.genes)
Heart_failure_4 <- Heart_failure_4[which(rownames(Heart_failure_4) %!in% Delete_genes),]
#保存
setwd('I:\\LY\\geo\\GSE247468_RAW')
saveRDS(Heart_failure_4, file = "./Heart_failure_4.rds")
#删除占内存的
rm(list = ls())
gc()

