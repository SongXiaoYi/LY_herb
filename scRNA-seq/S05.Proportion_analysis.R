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

#### load_data
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE247468')
Idents(pbmc) <- pbmc$celltype_first

metaData <- pbmc@meta.data

df <- metaData %>% group_by(Status,ID,celltype_first,.drop = FALSE) %>% summarise(n=n()) %>% mutate(freq = n/sum(n))
iOrd <-c('Control_1', 'Control_2', 'Heart_failure_1', 'Heart_failure_2', 'Heart_failure_3', 'Heart_failure_4')# in the order of first occurrence in df.
df$stim <- factor(df$ID,levels = iOrd)
#df[which(df$stim %in% c('CXN-T','LBZ-T','LFM_T_5','LJC-T','LSX-T','P1_T','P1_T2','P2_T','P3_T','P4_T','P5_T','P5_T2')), 'Tissue'] <- "T"
#df[which(df$stim %in% c('CXN-L','LBZ-L','LFM_LN_5','LJC-L','LSX-IIILN-L','LSX-IILN-L','P1_L','P2_L','P4_L','P5_L')), 'Tissue'] <- "L"
#df[which(df$stim %in% c('CXN-M','LBZ-M','LFM_M_5','LJC-M','P1_M','P2_M','P3_M','P4_M','P5_M')), 'Tissue'] <- "M"
df$Status <- factor(df$Status, levels = c('Case','Control'))
#df <- df[which(df$B_Subcluster %in% c('PC1','PC2','PC3','PC4','PC5','PC6')),]
my_comparisons <- list(c("Case", "Control"))
ggplot(df, aes(x = Status, y = freq, fill = celltype_first)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2,size=2) + 
    theme_classic()+ylab("Cell proportion (%)") + facet_grid(~celltype_first) + stat_compare_means(method = "wilcox.test")

my_comparisons <- list(c('Case','Control'))
p2 <- ggplot(df, aes(x = Status, y = freq, fill = Status)) + 
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15,size=2,alpha = 0.5) + 
    ylab("Cell proportion (%)") + 
    facet_wrap(~celltype_first, nrow = 1) + 
    #stat_compare_means(comparisons = my_comparisons,label.y = 0.55, method = "wilcox.test", label = "p.format") + 
    theme(legend.position = 'none', panel.border = element_rect(linewidth = 1, fill = NA), strip.background = element_rect(fill = "lightblue")) + xlab('')

setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE247468')
pdf(file = 'Cell_proportion.pdf', width = 9.4, height = 2.95)
p2
dev.off()


