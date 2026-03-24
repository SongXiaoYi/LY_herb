#####################
library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(patchwork)
library(ggplot2)
library(ggbump)
library(knitr)
library(plot1cell)
###########################

#### load_data
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE247468')
pbmc <- readRDS('./data_preparation.rds')
DimPlot(pbmc, reduction = "umap", pt.size = 0.1,label = TRUE,label.size = 6) + NoLegend()

################### Annotation conduct
################ Cardiac Fibroblasts (PT)
Genes <- c('COL1A1','DCN')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ Loopp of Henle cells (LOH)
Genes <- c('CLDN16','UMOD')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ Intercalated cells (IC)
Genes <- c('DMRT2','SLC26A7')
Genes <- c('ATP6V1G3','SLC26A7')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ Mesangial cells (Mes)
Genes <- c('PDGFRB','MYL9')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ Macrophages (Mac)
Genes <- c('CD68','CD163')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ Mono
Genes <- c('FCN1','CD14')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ Principal cells (PC)
Genes <- c('AQP3','AQP2')
Genes <- c('AQP2','FXYD4')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ ECs
Genes <- c('PECAM1','FLT1')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ Epithelial 
Genes <- c('CLDN4','KRT8')
Genes <- c('KRT8','KRT18')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ DT 
Genes <- c('SLC12A3')
FeaturePlot(pbmc, features = Genes, label = TRUE,label.size = 6)

################ B cells
Genes <- c('IGKC','MS4A1')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ T cells
Genes <- c('CD3D','NKG7')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

Genes <- c('ENTPD1')
FeaturePlot(pbmc, features = Genes, label = TRUE,label.size = 6)

################ Podocytes
Genes <- c('NPHS2','WT1')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ Fibro
Genes <- c('COL1A1','DCN')
FeaturePlot(pbmc, features = Genes, blend = TRUE,label = TRUE,label.size = 6)

################ Pericyte
Genes <- c('RGS5')
FeaturePlot(pbmc, features = Genes, label = TRUE,label.size = 6)

####  First Annotation
## Cardiac Fibroblasts       0,6,7,12
## Endothelial cells         1,4,11,15,16
## Pericytes                 2
## Neuronal                  13
## T cell                    3,14,18
## Monocytes / Macrophages   9
## B cell                    19
## NK Cells                  10                  
## Neutrophils               5,17
## Cardiomyocyte             8

Idents(pbmc) <- pbmc$seurat_clusters
new.cluster.ids <- c("Fibroblast", "EC", "Pericyte", "T cell", "EC", "Neutrophil", "Fibroblast", "Fibroblast", "Cardiomyocyte", "Macrophage", "NK Cell", "EC",
                     "Fibroblast", "Neuronal", "T cell", "EC", "EC", "Neutrophil", "T cell", "B cell")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$celltype_first <- Idents(pbmc)

########################################## Annotation advance

#####################################################
pbmc$Cluster <- pbmc$celltype_first
Idents(pbmc) <- pbmc$Cluster
####################### Set Annotation

#################### 全部改大写
pbmc$Status <- pbmc$Status
pbmc$Sample <- pbmc$ID

####################################### circ_plot
#colnames(pbmc@meta.data)
circ_data <- prepare_circlize_data(pbmc, scale = 0.65)
set.seed(2023)

cluster_colors<- c("#9900CC", "#3C5488CC", "#0071C1", "#00B151", "#92D14F", "#ED7D33", "#997303", "#9D470c", "#C50001", "#660000")
Status_colors <- c('#000080', '#8B0000')
Sample_colors <- rand_color(length(names(table(pbmc$ID))))


################### 画图函数调整
plot_circlize <- function(
  data_plot,
  do.label = T,
  contour.levels = c(0.2,0.3),
  pt.size = 0.5,
  kde2d.n = 1000,
  contour.nlevels = 100,
  bg.color='#F9F2E4',
  col.use=NULL,
  label.cex = 0.5,
  repel=FALSE
  ) {
  data_plot %>%
    dplyr::group_by(Cluster) %>%
    summarise(x = median(x = x), y = median(x = y)) -> centers
  z <- MASS::kde2d(data_plot$x, data_plot$y, n=kde2d.n)
  celltypes<-names(table(data_plot$Cluster))
  cell_colors <- scales::hue_pal()(length(celltypes))
  if(!is.null(col.use)){
    cell_colors=col.use
    col_df<-data.frame(Cluster=celltypes, color2=col.use)
    cells_order<-rownames(data_plot)
    data_plot<-merge(data_plot, col_df, by="Cluster")
    rownames(data_plot)<-data_plot$cells
    data_plot<-data_plot[cells_order,]
    data_plot$Colors<-data_plot$color2
  }
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding=c(0,0,0,0), track.margin=c(0.01,0),"track.height" = 0.02, gap.degree =c(rep(2, (length(celltypes)-1)),12),points.overflow.warning=FALSE)
  circos.initialize(sectors =  data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y=data_plot$dim2, track.height = 0.022,bg.border=NA,panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter,
                CELL_META$cell.ylim[2]+ mm_y(4),
                CELL_META$sector.index,
                cex=0.8, col = 'black', facing = "bending.inside", niceFacing = T)
    circos.axis(labels.cex = 0.4, col = 'black', labels.col =  'black')
  })
  for(i in 1:length(celltypes)){
    dd<-data_plot[data_plot$Cluster==celltypes[i],]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), y1 = 0.5, col = cell_colors[i],  lwd=7, sector.index = celltypes[i])
  }
  text(x = 1, y=0.1, labels = "Cell type", cex = 0.7, col = 'black',srt=-87)
  points(data_plot$x,data_plot$y, pch = 19, col = alpha(data_plot$Colors,0.2), cex = pt.size);
  contour(z, drawlabels=F, nlevels= 100, levels = contour.levels,col = '#ae9c76', add=TRUE)
  if(do.label){
    if(repel){
      textplot(x=centers$x, y=centers$y, words =  centers$Cluster,cex = label.cex, new = F,show.lines=F)
    } else {
      text(centers$x,centers$y, labels=centers$Cluster, cex = label.cex, col = 'black')
    }
  } 
}

add_track <- function(
  data_plot, 
  group, 
  track_num, 
  colors = NULL
  ){
  if(track_num<2){
    stop("The first track is the cluster track. Please change the track_num to a value greater than 1")
  }
  circos.track(data_plot$Cluster, data_plot$x_polar2, y=data_plot$dim2, track.height = 0.022, bg.border=NA)
  celltypes<-names(table(data_plot$Cluster))
  group_names<-names(table(data_plot[,group]))
  if(is.null(colors)){
    col_group = scales::hue_pal()(length(group_names))
  } else {
    col_group = colors
  }
  for(i in 1:length(celltypes)) {
    data_plot_cl<-data_plot[data_plot$Cluster==celltypes[i],]
    dat_seg<-get_segment(data_plot_cl, group = group)
    dat_seg2<-c(dat_seg[-1]-1, nrow(data_plot_cl))
    scale_factor<-max(data_plot_cl$x_polar2)/nrow(data_plot_cl)
    dat_seg<-scale_factor*dat_seg
    dat_seg2<-scale_factor*dat_seg2
    circos.segments(x0 = dat_seg, y0 = 0, x1 = dat_seg2, y1 = 0, col = col_group, sector.index = celltypes[i], lwd=7)
  }
  text(x = (1-0.04*(track_num-1)), y=0.09, labels = group, cex = 0.7, col = 'black',srt=-87)
}


########################## 绘图-----主图
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE247468')
#png(filename =  'circlize_plot.pdf', width = 6, height = 6,units = 'in', res = 300)
pdf(file="./Tissue_specific_marker.pdf",width = 6, height = 6)
plot_circlize(circ_data,do.label = T, pt.size = 0.05, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = TRUE, label.cex = 1,contour.levels = c(0.1,0.2))
add_track(circ_data, group = "Status", colors = Status_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "Sample", colors = Sample_colors, track_num = 5) ## can change it to one of the columns in the meta data of your seurat object
dev.off()

pdf(file="./Umap_status.pdf",width = 3.75, height = 3.7)
DimPlot(pbmc, reduction = "umap", group.by = 'Status',pt.size = 0.1,label.size = 6) + NoLegend()
dev.off()

pdf(file="./Umap_sample.pdf",width = 3.75, height = 3.7)
DimPlot(pbmc, reduction = "umap", group.by = 'ID',pt.size = 0.1,label.size = 6) + NoLegend()
dev.off()

#ColourC <- c("#9900CC", "#4472C7", "#0071C1", "#00B151", "#92D14F", "#ED7D33",  "#C50001", "#9D470c", "#997303")
##ColourC <- c("#FED966", "#4472C7", "#0071C1", "#00B151", "#92D14F", "#ED7D33", "#636363", "#C50001", "#9D470c", "#997303")

#setwd('/GPUFS/sysu_hshan_9/atlas/Figure_for_article/Figure1')
#pdf(file = 'UMAP_celltype.pdf', width = 4.1, height = 3.6)
#DimPlot(pbmc, reduction = "umap", label = T, pt.size = 0.02) + theme_bw() + NoLegend() + 
#theme(panel.grid = element_blank(),panel.border = element_rect(fill=NA,color="black",size = 1,linetype="solid")) + scale_colour_manual(values = ColourC)
#dev.off()

########################## 绘图-----分图
#### 2. UMAP for all Status
Idents(pbmc) <- pbmc$Status
ColourStatus <- c('#000080','#8B0000')
ColourStatus <- rev(ColourStatus)
setwd('H:\\Mendilian\\RNA_experiment\\Annotation\\Annotation_figure')
pdf(file = 'UMAP_ColourStatus.pdf', width = 4.1, height = 3.6)
DimPlot(pbmc, reduction = "umap", pt.size = 0.02)  + theme_void() + scale_colour_manual(values = ColourStatus)
dev.off()

#pdf(file = 'UMAP_Month_legand.pdf', width = 4.1, height = 3.6)
#DimPlot(pbmc, reduction = "umap", pt.size = 0.1,) +  scale_colour_manual(values = ColourM)
#dev.off()

#### 3. UMAP for all Gender
Idents(pbmc) <- pbmc$Gender
ColourGender <- c('#FF1493', '#1E90FF')

setwd('H:\\Mendilian\\RNA_experiment\\Annotation\\Annotation_figure')
pdf(file = 'UMAP_ColourGender.pdf', width = 4.1, height = 3.6)
DimPlot(pbmc, reduction = "umap", pt.size = 0.02)  + theme_void() + scale_colour_manual(values = ColourGender)
dev.off()

#pdf(file = 'UMAP_tissue_legand.pdf', width = 4.1, height = 3.6)
#DimPlot(pbmc, reduction = "umap", pt.size = 0.02)  + scale_colour_manual(values = ColourT)
#dev.off()


#### 4. UMAP for Cohort
Idents(pbmc) <- pbmc$Cohort

ColourCohort <- c('#FF0000', '#008000', '#0000FF')

setwd('H:\\Mendilian\\RNA_experiment\\Annotation\\Annotation_figure')
pdf(file = 'UMAP_ColourCohort.pdf', width = 4.1, height = 3.6)
DimPlot(pbmc, reduction = "umap", pt.size = 0.02)  + theme_void() + scale_colour_manual(values = ColourCohort)
dev.off()

#pdf(file = 'UMAP_AV_legand.pdf', width = 4.1, height = 3.6)
#DimPlot(pbmc, reduction = "umap", pt.size = 0.1,)  + scale_colour_manual(values = ColourAV)
#dev.off()

#### 5. UMAP for Gender
Idents(pbmc) <- pbmc$Sample

ColourSample <- Sample_colors

setwd('H:\\Mendilian\\RNA_experiment\\Annotation\\Annotation_figure')
pdf(file = 'UMAP_ColourSample.pdf', width = 4.1, height = 3.6)
DimPlot(pbmc, reduction = "umap", pt.size = 0.02)  + theme_void() + scale_colour_manual(values = ColourSample)
dev.off()


################################################################ Expression

####  First Annotation
## Cardiac Fibroblasts       0,6,7,12
## Endothelial cells         1,4,11,15,16
## Pericytes                 2
## Neuronal                  13
## T cell                    3,14,18
## Monocytes / Macrophages   9
## B cell                    19
## NK Cells                  10                  
## Neutrophils               5,17
## Cardiomyocyte             8


Idents(pbmc) <- pbmc$celltype_first

feature_plot <- c('LUM','PECAM1','MYH11','MPZ','CD3D','CD68','MS4A1','NKG7','FCGR3B','MYH7')
cell_plot <- c('Fibroblast','EC','Pericyte','Neuronal','T cell','Macrophage','B cell','NK Cell','Neutrophil','Cardiomyocyte')

setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE247468')
#png(filename =  'vlnplot_multiple_genes.png', width = 6, height = 6,units = 'in', res = 300)
pdf(file =  'vlnplot_multiple_genes.pdf', width = 6.75, height = 5.95)
complex_vlnplot_multiple(pbmc, features = feature_plot, celltypes = cell_plot, group = "Status", add.dot=F, font.size = 7)
dev.off()

############################## 保存注释
# 创建一个3行2列的矩阵
pbmc@assays$RNA@scale.data <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
setwd('H:\\Mendilian\\RNA_experiment\\Annotation')
saveRDS(pbmc,'./Annotation_data.rds')
