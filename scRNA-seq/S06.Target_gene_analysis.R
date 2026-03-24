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
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(harmony)
library(patchwork)
library(ggplot2)
library(ggpubr)
###########################
library(Seurat)
library(dplyr)
library(harmony)
library(symphony)
library(ggplot2)
library(RColorBrewer)
#library(scRepertoire)
library(patchwork)

# stacked violin plot (adapted from https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/)
# remove the x-axis text and tick
# plot.margin to adjust the white space between each plot.
# ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0), 
          plot.title= element_blank(),
          axis.title.x = element_blank(), 
          plot.margin = plot.margin)
  return(p)
}

# extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

# main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = 9), 
          axis.ticks.x = element_line())
  
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}



feature_plot <- c('MKI67','CDK1','TOP2A','SPAG5','CCNA2')
VlnPlot(pbmc, features = 'MKI67')
#######################################################
pbmc_sub <- subset(pbmc,idents = 'Fibroblast')
Target_gene_matrix <- GetAssayData(pbmc_sub[["RNA"]], layer = "data")
Target_gene_matrix <- Target_gene_matrix[feature_plot,]
Target_gene_matrix <- t(Target_gene_matrix)
Target_gene_matrix <- as.matrix(Target_gene_matrix)
meta_datas <- pbmc_sub@meta.data
##################################
Target_gene_matrix <- as.data.frame(Target_gene_matrix)
Target_gene_matrix$group <- meta_datas$Status
########################################### Plot
MKI67_matrix <- Target_gene_matrix[Target_gene_matrix$MKI67 != 0 ,]
newdata <- MKI67_matrix
newdata$Group <- factor(newdata$group, level = c('Case','Control'))
my_comparisons <- list(c("Case", "Control"))
a=ggboxplot(newdata, x="group", y="MKI67", color = "group", add = "jitter", shape = 16,
            ylab="MKI67 expression", xlab="",outlier.shape = NA,bxp.errorbar = TRUE)
a=a+rotate_x_text(51) + theme(plot.title = element_text(size=12,hjust=0.5))
a + stat_compare_means(comparisons = my_comparisons) + NoLegend()
                         
pdf(file = 'MKI67.pdf', width = 2.4, height = 3.65)
a + stat_compare_means(comparisons = my_comparisons) + NoLegend()
dev.off()

                         
###################################
CDK1_matrix <- Target_gene_matrix[Target_gene_matrix$CDK1 != 0 ,]
newdata <- CDK1_matrix
newdata$Group <- factor(newdata$group, level = c('Case','Control'))
my_comparisons <- list(c("Case", "Control"))
a=ggboxplot(newdata, x="group", y="CDK1", color = "group", add = "jitter", shape = 16,
            ylab="CDK1 expression", xlab="",outlier.shape = NA,bxp.errorbar = TRUE)
a=a+rotate_x_text(51) + theme(plot.title = element_text(size=12,hjust=0.5))
a + stat_compare_means(comparisons = my_comparisons) + NoLegend()

pdf(file = 'CDK1.pdf', width = 2.4, height = 3.65)
a + stat_compare_means(comparisons = my_comparisons) + NoLegend()
dev.off()
###################################
TOP2A_matrix <- Target_gene_matrix[Target_gene_matrix$TOP2A != 0 ,]
newdata <- TOP2A_matrix
newdata$Group <- factor(newdata$group, level = c('Case','Control'))
my_comparisons <- list(c("Case", "Control"))
a=ggboxplot(newdata, x="group", y="TOP2A", color = "group", add = "jitter", shape = 16,
            ylab="TOP2A expression", xlab="",outlier.shape = NA,bxp.errorbar = TRUE)
a=a+rotate_x_text(51) + theme(plot.title = element_text(size=12,hjust=0.5))
a + stat_compare_means(comparisons = my_comparisons) + NoLegend()

pdf(file = 'TOP2A.pdf', width = 2.4, height = 3.65)
a + stat_compare_means(comparisons = my_comparisons) + NoLegend()
dev.off()                         
###################################
SPAG5_matrix <- Target_gene_matrix[Target_gene_matrix$SPAG5 != 0 ,]
newdata <- SPAG5_matrix
newdata$Group <- factor(newdata$group, level = c('Case','Control'))
my_comparisons <- list(c("Case", "Control"))
a=ggboxplot(newdata, x="group", y="SPAG5", color = "group", add = "jitter", shape = 16,
            ylab="SPAG5 expression", xlab="",outlier.shape = NA,bxp.errorbar = TRUE)
a=a+rotate_x_text(51) + theme(plot.title = element_text(size=12,hjust=0.5))
a + stat_compare_means(comparisons = my_comparisons) + NoLegend()

pdf(file = 'SPAG5.pdf', width = 2.4, height = 3.65)
a + stat_compare_means(comparisons = my_comparisons) + NoLegend()
dev.off()                         
###################################
CCNA2_matrix <- Target_gene_matrix[Target_gene_matrix$CCNA2 != 0 ,]
newdata <- CCNA2_matrix
newdata$Group <- factor(newdata$group, level = c('Case','Control'))
my_comparisons <- list(c("Case", "Control"))
a=ggboxplot(newdata, x="group", y="CCNA2", color = "group", add = "jitter", shape = 16,
            ylab="CCNA2 expression", xlab="",outlier.shape = NA,bxp.errorbar = TRUE)
a=a+rotate_x_text(51) + theme(plot.title = element_text(size=12,hjust=0.5))
a + stat_compare_means(comparisons = my_comparisons) + NoLegend()

pdf(file = 'CCNA2.pdf', width = 2.4, height = 3.65)
a + stat_compare_means(comparisons = my_comparisons) + NoLegend()
dev.off() 







                         

