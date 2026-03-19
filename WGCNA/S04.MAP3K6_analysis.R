# 加载包
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(homologene)
library(org.Hs.eg.db)
library(GSEAtopics)
#####################################
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\WGCNA')
Genelist <- read.csv('./IMConn_rank.csv')
Genelist <- Genelist$Gene
TurnGenelist <- mouse2human(Genelist)

# ===================== 核心函数 =====================
#' 读取MAP3K6专属GMT文件
read_map3k6_gmt <- function(gmt_path = "MAP3K6_hypertension_heart.gmt") {
  # 读取GMT文件
  gmt_data <- read.gmt(gmt_path)
  # 转换为通路-基因列表
  gmt_dict <- split(gmt_data$gene, gmt_data$term)
  # 过滤空通路
  gmt_dict <- gmt_dict[sapply(gmt_dict, length) >= 2]
  
  cat(sprintf("成功读取MAP3K6 GMT文件，共解析到 %d 个通路\n", length(gmt_dict)))
  cat("包含的通路：\n")
  print(names(gmt_dict))
  return(gmt_dict)
}

#' 针对MAP3K6相关基因做富集分析
map3k6_enrichment <- function(
  query_genes,  # 待富集基因（如你的目标基因列表）
  gmt_dict,
  pvalue_cutoff = 0.05
) {
  # 数据预处理：去重、大写、去空
  query_genes <- unique(toupper(query_genes))
  query_genes <- query_genes[query_genes != ""]
  cat(sprintf("\n待富集基因总数：%d\n", length(query_genes)))
  
  # 背景基因集（GMT中所有基因）
  background_genes <- unique(unlist(gmt_dict))
  cat(sprintf("背景基因总数：%d\n", length(background_genes)))
  
  # 检查交集
  overlap_genes <- intersect(query_genes, background_genes)
  cat(sprintf("与GMT通路的交集基因数：%d\n", length(overlap_genes)))
  if (length(overlap_genes) > 0) {
    cat("交集基因：", paste(overlap_genes, collapse = ", "), "\n")
  } else {
    stop("无交集基因！请检查输入基因名格式")
  }
  
  # 执行富集分析
  enrich_result <- enricher(
    gene = query_genes,
    TERM2GENE = stack(gmt_dict),
    universe = background_genes,
    pAdjustMethod = "BH",  # FDR校正
    pvalueCutoff = pvalue_cutoff,
    minGSSize = 2,         # 适配小基因集
    maxGSSize = 50
  )
  
  # 整理结果
  if (!is.null(enrich_result)) {
    result_df <- as.data.frame(enrich_result) %>%
      # 计算富集倍数
      mutate(
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
        BgRatio_num = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
        Fold_Enrichment = GeneRatio_num / BgRatio_num
      ) %>%
      select(-GeneRatio_num, -BgRatio_num) %>%
      arrange(p.adjust)
    
    cat(sprintf("富集分析完成，共得到 %d 个显著通路\n", nrow(result_df)))
    return(result_df)
  } else {
    cat("无显著富集结果\n")
    return(data.frame())
  }
}

#' 绘制MAP3K6富集分析气泡图
plot_map3k6_enrich <- function(result_df, output_path = "MAP3K6_enrich_bubble.png") {
  # 排序并绘图
  p <- ggplot(result_df, aes(x = Fold_Enrichment, y = reorder(Description, Fold_Enrichment))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    labs(
      x = "Fold Enrichment (富集倍数)",
      y = "Pathway (通路)",
      size = "Gene Count (基因数)",
      title = "MAP3K6 Related Pathway Enrichment (Hypertension-Heart Disease)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )
  
  # 保存图片
  ggsave(output_path, p, width = 10, height = 6, dpi = 300)
  cat("富集气泡图已保存至：", output_path, "\n")
  return(p)
}

# ===================== 主运行流程 =====================
# 1. 读取MAP3K6专属GMT文件
setwd('E:\\CQT2026012704-F001_20260313\\CQT2026012704-F001_20260313\\Experiment\\Pathway_enrichment')
gmt_dict <- read_map3k6_gmt("MAP3K6_hypertension_heart.gmt")

# 2. 输入你的目标基因列表（示例：MAP3K6+高血压心脏病相关基因）
# 替换为你实际的基因列表（如差异表达基因、候选基因）
query_genes <- TurnGenelist$humanGene
                               
# 3. 执行富集分析
enrich_result <- map3k6_enrichment(query_genes, gmt_dict)

# 4. 查看并保存结果
if (nrow(enrich_result) > 0) {
  # 打印结果
  cat("\nMAP3K6富集分析结果：\n")
  print(enrich_result[, c("Description", "Count", "Fold_Enrichment", "p.adjust")])
  
  # 保存结果到CSV
  write.csv(enrich_result, "MAP3K6_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
  cat("富集结果已保存至：MAP3K6_enrichment_results.csv\n")
  
  # 绘制气泡图
  plot_map3k6_enrich(enrich_result)
}
