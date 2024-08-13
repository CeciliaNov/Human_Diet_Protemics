# 加载所需的包
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(dplyr)
library(tibble)
library(enrichplot)
library(STRINGdb)
library(igraph)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(org.Mm.eg.db)

theme="E"
# 创建输出目录
dirname <- "20240811commonpathway"
dir.create(dirname)

# 人HFD差异分析读取数据
dep <- read.csv("D:/Cycy1/R_Proj/Diet/diet_R/diet_dep.csv")

# 计算hfhc diet中的差异蛋白质数据
sel <- dep[dep$diet == theme, ]

plasma_diff_proteins <- sel %>%
  rename_with(~ paste0("log2_", .), starts_with("fd")) %>%
  rowwise() %>%
  mutate(
    fd_values = list(exp(c_across(starts_with("log2_fd")))),
    fd_geo_mean = exp(mean(log(unlist(fd_values)), na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(
    logFC = log2(ifelse(fd_geo_mean > 0, fd_geo_mean, 1e-1)),
    pvalue = min(c_across(starts_with("pval")), na.rm = TRUE)
  ) %>%
  select(SYMBOL, logFC, pvalue) %>%
  filter(abs(logFC) > 0.26)

#动物HFD主动脉

# 计算不同饮食中的差异蛋白质数据
library(readr)
dep <- read.csv("D:/Cycy1/R_Proj/Diet/diet_urine_R/dep_diet.csv")
calculate_diff_proteins <- function(data, diet_label) {
  sel <- data[data$diet == diet_label, ]
  diff_proteins <- sel %>%
    rename_with(~ paste0("log2_", .), starts_with("fd")) %>%
    rowwise() %>%
    mutate(
      fd_values = list(exp(c_across(starts_with("log2_fd")))),
      fd_geo_mean = exp(mean(log(unlist(fd_values)), na.rm = TRUE))
    ) %>%
    ungroup() %>%
    mutate(
      logFC = log2(ifelse(fd_geo_mean > 0, fd_geo_mean, 1e-1)),
      pvalue = min(c_across(starts_with("pval")), na.rm = TRUE),
      UNIPROT=sel$protein
    ) %>%
    select(UNIPROT, logFC, pvalue) %>%
    filter(abs(logFC) > 0.26)
  return(diff_proteins)
}

urine_diff_proteins <- calculate_diff_proteins(dep, theme)



# ######################################################
library(clusterProfiler)
library(org.Hs.eg.db)

# 示例数据：人类差异蛋白的 SYMBOL 列表
plasma_proteins_list <- plasma_diff_proteins$SYMBOL

# 将 SYMBOL 转换为 ENTREZID
plasma_gene_symbols <- bitr(plasma_proteins_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 执行 GO 富集分析
plasma_go_enrichment <- enrichGO(gene = plasma_gene_symbols$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

# 执行 KEGG 富集分析
plasma_kegg_enrichment <- enrichKEGG(gene = plasma_gene_symbols$ENTREZID, organism = 'hsa', pAdjustMethod = "BH", qvalueCutoff = 0.05)

# 查看富集结果
head(plasma_go_enrichment)
head(plasma_kegg_enrichment)


#####################################
urine_proteins_list <- urine_diff_proteins$UNIPROT

# 将 SYMBOL 转换为 ENTREZID
urine_gene_symbols <- bitr(urine_proteins_list, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 执行 GO 富集分析
urine_go_enrichment <- enrichGO(gene = urine_gene_symbols$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

# 执行 KEGG 富集分析
urine_kegg_enrichment <- enrichKEGG(gene = urine_gene_symbols$ENTREZID, organism = 'hsa', pAdjustMethod = "BH", qvalueCutoff = 0.05)

# 查看富集结果
head(urine_go_enrichment)
head(urine_kegg_enrichment)

# 加载所需的库
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)


# 提取 GO 富集分析结果
plasma_go_results <- as.data.frame(plasma_go_enrichment)
urine_go_results <- as.data.frame(urine_go_enrichment)

# 提取 KEGG 富集分析结果
plasma_kegg_results <- as.data.frame(plasma_kegg_enrichment)
urine_kegg_results <- as.data.frame(urine_kegg_enrichment)

# 找到共同 GO 通路
common_go_pathways <- intersect(plasma_go_results$ID, urine_go_results$ID)
# 提取前5个共同通路
top_common_go <- head(common_go_pathways, 5)

# 找到不同 GO 通路
unique_plasma_go_pathways <- setdiff(plasma_go_results$ID, urine_go_results$ID)
unique_urine_go_pathways <- setdiff(urine_go_results$ID, plasma_go_results$ID)
# 提取前5个不同通路
top_unique_plasma_go <- head(unique_plasma_go_pathways, 5)
top_unique_urine_go <- head(unique_urine_go_pathways, 5)

# 标记通路类型
plasma_go_results$Type <- ifelse(plasma_go_results$ID %in% top_common_go, "Common", "Unique to plasma")
urine_go_results$Type <- ifelse(urine_go_results$ID %in% top_common_go, "Common", "Unique to urine")

# 合并数据框以便于可视化
all_go_pathways <- rbind(
  plasma_go_results %>% filter(ID %in% c(top_common_go, top_unique_plasma_go)),
  urine_go_results %>% filter(ID %in% c(top_common_go, top_unique_urine_go)) %>% mutate(Type = ifelse(ID %in% top_common_go, "Common", "Unique to urine"))
)

# 找到共同 KEGG 通路
common_kegg_pathways <- intersect(plasma_kegg_results$ID, urine_kegg_results$ID)
# 提取前5个共同通路
top_common_kegg <- head(common_kegg_pathways, 5)

# 找到不同 KEGG 通路
unique_plasma_kegg_pathways <- setdiff(plasma_kegg_results$ID, urine_kegg_results$ID)
unique_urine_kegg_pathways <- setdiff(urine_kegg_results$ID, plasma_kegg_results$ID)
# 提取前5个不同通路
top_unique_plasma_kegg <- head(unique_plasma_kegg_pathways, 5)
top_unique_urine_kegg <- head(unique_urine_kegg_pathways, 5)

# 标记通路类型
plasma_kegg_results$Type <- ifelse(plasma_kegg_results$ID %in% top_common_kegg, "Common", "Unique to plasma")
urine_kegg_results$Type <- ifelse(urine_kegg_results$ID %in% top_common_kegg, "Common", "Unique to urine")

# 合并数据框以便于可视化
all_kegg_pathways <- rbind(
  plasma_kegg_results %>% filter(ID %in% c(top_common_kegg, top_unique_plasma_kegg)),
  urine_kegg_results %>% filter(ID %in% c(top_common_kegg, top_unique_urine_kegg)) %>% mutate(Type = ifelse(ID %in% top_common_kegg, "Common", "Unique to urine"))
)
library(cowplot)
library(ggplot2)
# 创建 GO 通路的条形图
go_plot <- ggplot(all_go_pathways, aes(y = reorder(Description, p.adjust), x = -log10(p.adjust), fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 14), # 调整x轴文字大小
    axis.text.y = element_text(size = 14), # 调整y轴文字大小
    axis.title.x = element_text(size = 16), # 调整x轴标题大小
    axis.title.y = element_text(size = 16), # 调整y轴标题大小
    plot.title = element_text(size = 18), # 调整标题大小
    legend.text = element_text(size = 14), # 调整图例文字大小
    legend.title = element_text(size = 16) # 调整图例标题大小
  ) +
  labs(title = "Top 5 Common and Unique GO Pathways", x = "-log10(p.adjust)", y = "GO Pathway", fill = "Type") +
  scale_fill_manual(values = c("Common" = "darkblue", "Unique to plasma" = "darkred", "Unique to urine" = "darkgreen"))

# 创建 KEGG 通路的条形图
kegg_plot <- ggplot(all_kegg_pathways, aes(y = reorder(Description, p.adjust), x = -log10(p.adjust), fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 14), # 调整x轴文字大小
    axis.text.y = element_text(size = 14), # 调整y轴文字大小
    axis.title.x = element_text(size = 16), # 调整x轴标题大小
    axis.title.y = element_text(size = 16), # 调整y轴标题大小
    plot.title = element_text(size = 18), # 调整标题大小
    legend.text = element_text(size = 14), # 调整图例文字大小
    legend.title = element_text(size = 16) # 调整图例标题大小
  ) +
  labs(title = "Top 5 Common and Unique KEGG Pathways", x = "-log10(p.adjust)", y = "KEGG Pathway", fill = "Type") +
  scale_fill_manual(values = c("Common" = "darkblue", "Unique to plasma" = "darkred", "Unique to urine" = "darkgreen"))

# 拼接图像并保存
combined_plot <- plot_grid(go_plot, kegg_plot, nrow = 2)

# 保存图像，调整图形尺寸
ggsave(paste0(dirname,"/",theme,"_combined_pathways_plot.tiff"), combined_plot, dpi = 300, width = 12, height = 12)
