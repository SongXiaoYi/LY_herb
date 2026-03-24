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
###################################
pbmc_sub <- subset(pbmc,idents = 'Fibroblast')
#########################################
Idents(pbmc_sub) <- pbmc_sub$Status
markers <- FindMarkers(pbmc_sub, ident.1 = 'Case', ident.2 = 'Control', min.pct = 0.1, logfc.threshold = 0)
###################################
###################################### 
####### WGCNA-PRAD
library(dplyr)
library(stringr)
library(patchwork)
library(knitr)
library(msigdbr) 
library(fgsea)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(limma)
library(clusterProfiler)
#library(DESeq2)
library(ggplotify)
library(cowplot)
library(edgeR)
library(tinyarray)
library(gridExtra)
library(Seurat)
library(WGCNA)
library(reshape2)
library('org.Mm.eg.db')
library('clusterProfiler')
library(homologene)
library(org.Hs.eg.db)
library(GSEAtopics)
#####################################
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\WGCNA')
TurnGenelist <- rownames(markers)
##################################
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\Pathway_enrichment')
#gmt <- getGmt('c2.cp.v2026.1.Hs.symbols.gmt')

read_gmt <- function(gmt_path) {
  gmt_list <- read.gmt(gmt_path)  # clusterProfiler内置函数
  # 转换为 通路名:基因列表 的格式
  gmt_dict <- split(gmt_list$gene, gmt_list$term)
  # 过滤空通路或基因数过少的通路
  gmt_dict <- gmt_dict[sapply(gmt_dict, length) >= 3]
  cat(sprintf("成功读取GMT文件，共解析到 %d 个通路\n", length(gmt_dict)))
  return(gmt_dict)
}

#' 执行基于GMT文件的富集分析
#' @param query_genes 待富集的基因列表（字符向量）
#' @param gmt_dict GMT文件解析后的列表
#' @param background_genes 背景基因集（可选）
#' @param pvalue_cutoff P值阈值
#' @param qvalue_cutoff Q值（FDR）阈值
#' @return 富集分析结果数据框
gmt_enrichment <- function(query_genes, 
                           gmt_dict, 
                           background_genes = NULL,
                           pvalue_cutoff = 0.05,
                           qvalue_cutoff = 0.05) {
  
  # 数据预处理：去重、过滤空值、统一大写
  query_genes <- unique(toupper(query_genes))
  query_genes <- query_genes[query_genes != ""]
  cat(sprintf("待富集基因总数: %d\n", length(query_genes)))
  
  # 设置背景基因集（默认使用GMT中所有基因）
  if (is.null(background_genes)) {
    background_genes <- unique(unlist(gmt_dict))
  } else {
    background_genes <- unique(toupper(background_genes))
    background_genes <- background_genes[background_genes != ""]
  }
  cat(sprintf("背景基因总数: %d\n", length(background_genes)))

  ABC <- stack(gmt_dict)
  ABC <- ABC[,c(2,1)]
  colnames(ABC) <- c('term','gene')
  # 执行富集分析（超几何检验）
  enrich_result <- enricher(
    gene = query_genes,
    TERM2GENE = ABC,  # 转换为clusterProfiler要求的格式
    universe = background_genes,
    pAdjustMethod = "BH",  # Benjamini-Hochberg校正（FDR）
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    minGSSize = 3,        # 通路最小基因数
    maxGSSize = 500       # 通路最大基因数（避免过宽通路）
  )
  print('完成')
  return(enrich_result)
}

#' 保存富集分析结果
#' @param result_df 富集结果数据框
#' @param output_path 输出文件路径（CSV格式）
save_enrichment <- function(result_df, output_path) {
  if (nrow(result_df) > 0) {
    write.csv(result_df, output_path, row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("结果已保存至: %s\n", output_path))
  } else {
    cat("无结果可保存\n")
  }
}

# ===================== 主运行流程 =====================
# 1. 设置参数
gmt_file <- "./c2.cp.v2026.1.Hs.symbols.gmt"  # 替换为你的GMT文件路径
gmt_file <- "./c2.cp.reactome.v2026.1.Hs.symbols.gmt"  # 替换为你的GMT文件路径
gmt_file <- "./c2.cp.wikipathways.v2026.1.Hs.symbols.gmt"  # 替换为你的GMT文件路径
gmt_file <- "./MAP3K6_hypertension_heart.gmt"
#output_file <- "enrichment_results_R.csv"  # 输出结果路径

# 2. 待富集的基因列表（替换为你的目标基因集）
# 示例：人类基因名（HGNC标准），实际使用时替换为你的基因列表
query_genes <- TurnGenelist
#query_genes <- TurnGenelist$humanGene

# 3. 读取GMT文件
gmt_dict <- read_gmt(gmt_file)
Gene_dictAAA = stack(gmt_dict)
                               
# 4. 执行富集分析
enrich_result <- gmt_enrichment(
  query_genes = query_genes,
  gmt_dict = gmt_dict,
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05
)
final_result <- as.data.frame(enrich_result@result)
# 5. 查看并保存结果
if (nrow(enrich_result) > 0) {
  # 查看前10个显著富集的通路
  cat("\n前10个显著富集的通路：\n")
  print(head(filter(enrich_result, Significant == "Yes"), 10))
  
  # 保存结果
  save_enrichment(enrich_result, output_file)
  
  # 可选：绘制富集气泡图（可视化结果）
  if (require("ggplot2")) {
    # 取前20个通路绘图
    plot_data <- head(enrich_result, 20)
    ggplot(plot_data, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "red", high = "blue") +
      labs(x = "Fold Enrichment", y = "Pathway", size = "Gene Count", color = "Adjusted P-value") +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8))
    ggsave("enrichment_bubble_plot.png", width = 10, height = 8, dpi = 300)
    cat("富集气泡图已保存至: enrichment_bubble_plot.png\n")
  }
}

setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\GEO\\Data\\GSE247468')
plot_data <- head(enrich_result, 20)
pdf(file = 'wikipathways_compare.pdf', width = 7.7, height = 4.1)
ggplot(plot_data, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Fold Enrichment", y = "Pathway", size = "Gene Count", color = "Adjusted P-value") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8)) + ylab('')
dev.off()






