library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(homologene)
###################
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\WGCNA')
Genelist <- read.csv('./IMConn_rank.csv')
Genelist <- Genelist$Gene

mouse2human_homologene <- function(mouse_genes) {
  # 核心转换函数：
  # inTax = 10090（小鼠TaxID），outTax = 9606（人类TaxID）
  result <- homologene(
    mouse_genes,
    inTax = 10090,
    outTax = 9606
  )
  
  # 整理结果列名（更直观）
  colnames(result) <- c("Mouse_Gene", "Human_Gene", "HID")
  # 去重+过滤空值
  result <- result[!is.na(result$Human_Gene) & result$Human_Gene != "", ]
  result <- unique(result)
  
  cat(sprintf("成功转换 %d / %d 个小鼠基因\n", nrow(result), length(mouse_genes)))
  return(result)
}

# 运行示例（你的目标基因）
mouse_genes <- Genelist
human_genes <- mouse2human_homologene(mouse_genes)
