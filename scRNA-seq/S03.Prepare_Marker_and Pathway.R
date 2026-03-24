#####################
library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(patchwork)
library(ggplot2)
library(ggbump)
library(knitr)
#library(plot1cell)
###########################
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE247468')
pbmc <- readRDS('./pbmc.rds')

############################## Identify Markers
logFCfilter=0.5
adjPvalFilter=0.05

pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]

sig.markers=pbmc.markers %>% group_by(cluster) %>% top_n(n = 50000, wt = avg_log2FC)

setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE247468')
write.csv(sig.markers,file="Seurat_cluster_markers.csv")

############################# Identify Pathway
#Pathway analysis
pbmc <- cerebroApp::getMarkerGenes(pbmc, groups = c('seurat_clusters'), assay = "RNA", organism = "hg")
gc()
pbmc <- cerebroApp::getEnrichedPathways(pbmc,databases = c("GO_Biological_Process_2021"),
                                        adj_p_cutoff = 0.05,
                                        max_terms = 100,
                                        URL_API = "http://maayanlab.cloud/Enrichr")

pathway <- pbmc@misc$enriched_pathways$cerebro_seurat_enrichr$seurat_clusters
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE247468')
write.csv(pathway,file="Seurat_cluster_pathway.csv")
