#加载包
#####################################
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(harmony)
library(patchwork)
library(ggplot2)
#####################################
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE247468')
##################
Control_1 <- readRDS("./Control_1.rds")
Control_2 <- readRDS("./Control_2.rds")
Heart_failure_1 <- readRDS("./Heart_failure_1.rds")
Heart_failure_2 <- readRDS("./Heart_failure_2.rds")
Heart_failure_3 <- readRDS("./Heart_failure_3.rds")
Heart_failure_4 <- readRDS("./Heart_failure_4.rds")

Control_1@meta.data$Doublet <- Control_1@meta.data$DF.classifications_0.25_0.005_330
Control_1 <- subset(Control_1,subset = Doublet == "Singlet")

Control_2 <- subset(Control_2,subset = Doublet == "Singlet")
Heart_failure_1 <- subset(Heart_failure_1,subset = Doublet == "Singlet")
Heart_failure_2 <- subset(Heart_failure_2,subset = Doublet == "Singlet")
Heart_failure_3 <- subset(Heart_failure_3,subset = Doublet == "Singlet")
Heart_failure_4 <- subset(Heart_failure_4,subset = Doublet == "Singlet")

cellid  <- c('Control_1','Control_2','Heart_failure_1','Heart_failure_2','Heart_failure_3','Heart_failure_4')
sceList <- list(Control_2,Heart_failure_1,Heart_failure_2,Heart_failure_3,Heart_failure_4)

pbmc.big <- merge(Control_1, y = sceList, add.cell.ids = cellid, project = "10X_Heart")
saveRDS(pbmc.big, file = "sce.anchors.rds")

rm(list = ls())
gc()

#################################
sce.integrated <- readRDS('./sce.anchors.rds')
sce.integrated <- NormalizeData(object = sce.integrated, normalization.method = "LogNormalize", scale.factor = 10000)
sce.integrated <- FindVariableFeatures(object = sce.integrated, selection.method = "vst", nfeatures = 2500)
all.genes <- rownames(sce.integrated)
sce.integrated <- ScaleData(sce.integrated, features = all.genes)
sce.integrated <- RunPCA(object= sce.integrated,npcs = 20,pc.genes=VariableFeatures(object = sce.integrated))   #这里PCA设置为前100个

#开始去批次
sce.integrated$orig.ident <- sce.integrated$ID
sce.integrated <- RunHarmony(sce.integrated,"orig.ident",dims.use = 1:20,max.iter.harmony = 15,max.iter.cluster = 20)   #去批次这一部看结果
sce.integrated <- RunUMAP(sce.integrated, dims = 1:20, min.dist = 0.2,reduction = "harmony")   #尝试调整min.dist会得到更好的结果
sce.integrated <- FindNeighbors(sce.integrated, reduction = "harmony",dims = 1:20)

sce.integrated <- FindClusters(sce.integrated, resolution = 0.4)   

#sce.integrated <- RunTSNE(sce.integrated, dims = 1:20, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)    #尝试调整max_iter会得到更好的结果
sce.integrated <- RunUMAP(sce.integrated, dims = 1:20, min.dist = 0.75)   #尝试调整min.dist会得到更好的结果

#画图
DimPlot(sce.integrated,reduction = "umap", label = T) + NoLegend()
saveRDS(sce.integrated,'./pbmc.rds')
