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
library(readxl)

# 创建输出目录
dirname <- "D:/Cycy1/R_Proj/Diet/diet_animal_R/202408HFDcommon"
dir.create(dirname)

# 人HFD差异分析读取数据
dep <- read.csv("D:/Cycy1/R_Proj/Diet/diet_R/diet_dep.csv")

# 计算hfhc diet中的差异蛋白质数据
sel_C <- dep[dep$diet == "C", ]

human_diff_proteins_hfhc <- sel_C %>%
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
mice <- read_excel("D:/Cycy1/R_Proj/Diet/diet_animal_R/差异分析结果/b31/2022-04-29_b17CD_HFD_17W_protein_ttest_diff.xlsx")

sel_P =mice

mice_diff_proteins_hfhc <- sel_P %>%
  mutate(
    UNIPROT = sel_P$...1,
    logFC=log2(sel_P$fd),
    pvalue=sel_P$p
  ) %>%
  select(UNIPROT, logFC,pvalue) %>%
  filter(abs(logFC) > 0.26)

s2e <- bitr(mice_diff_proteins_hfhc$UNIPROT, 
            fromType = "UNIPROT",
            toType = "SYMBOL",
            OrgDb = org.Mm.eg.db)#mouse

# 将 UNIPROT 转换为 SYMBOL
mice_diff_proteins_hfhc <- mice_diff_proteins_hfhc %>%
  left_join(s2e, by = "UNIPROT") %>%
  select(SYMBOL, logFC, pvalue)

# 将 SYMBOL 转换为大写
mice_diff_proteins_hfhc$UPSYMBOL <- toupper(mice_diff_proteins_hfhc$SYMBOL)


# 共同差异蛋白
common_proteins <- intersect(human_diff_proteins_hfhc$SYMBOL, mice_diff_proteins_hfhc$UPSYMBOL)

# 不同差异蛋白
unique_human_proteins <- setdiff(human_diff_proteins_hfhc$SYMBOL, mice_diff_proteins_hfhc$UPSYMBOL)
unique_mice_proteins <- setdiff(mice_diff_proteins_hfhc$SYMBOL, human_diff_proteins_hfhc$UPSYMBOL)

library(ggplot2)
library(dplyr)

# 创建数据框以用于可视化
proteins_comparison <- data.frame(
  Protein = c(common_proteins, unique_human_proteins, unique_mice_proteins),
  Type = c(rep("Common", length(common_proteins)), 
           rep("Unique to human", length(unique_human_proteins)), 
           rep("Unique to mice", length(unique_mice_proteins)))
)

# 绘制条形图
ggplot(proteins_comparison, aes(x = Protein, fill = Type)) +
  geom_bar(stat = "count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Comparison of Differential Proteins", x = "Protein", y = "Count")

# 加载所需的库
library(clusterProfiler)
library(org.Hs.eg.db)

# 示例数据：人类差异蛋白的 SYMBOL 列表
human_proteins_list <- human_diff_proteins_hfhc$SYMBOL

# 将 SYMBOL 转换为 ENTREZID
human_gene_symbols <- bitr(human_proteins_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 执行 GO 富集分析
human_go_enrichment <- enrichGO(gene = human_gene_symbols$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

# 执行 KEGG 富集分析
human_kegg_enrichment <- enrichKEGG(gene = human_gene_symbols$ENTREZID, organism = 'hsa', pAdjustMethod = "BH", qvalueCutoff = 0.05)

# 查看富集结果
head(human_go_enrichment)
head(human_kegg_enrichment)
# 加载所需的库
library(clusterProfiler)
library(org.Mm.eg.db)

# 示例数据：小鼠差异蛋白的 SYMBOL 列表
mice_proteins_list <- mice_diff_proteins_hfhc$SYMBOL

# 将 SYMBOL 转换为 ENTREZID
mice_gene_symbols <- bitr(mice_proteins_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# 执行 GO 富集分析
mice_go_enrichment <- enrichGO(gene = mice_gene_symbols$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

# 执行 KEGG 富集分析
mice_kegg_enrichment <- enrichKEGG(gene = mice_gene_symbols$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)

# 查看富集结果
head(mice_go_enrichment)
head(mice_kegg_enrichment)

# 加载所需的库
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


# 假设人类和小鼠的差异蛋白富集分析结果已存储在human_go_enrichment和mice_go_enrichment中

# 提取 GO 富集分析结果
human_go_results <- as.data.frame(human_go_enrichment)
mice_go_results <- as.data.frame(mice_go_enrichment)

# 提取 KEGG 富集分析结果
human_kegg_results <- as.data.frame(human_kegg_enrichment)
mice_kegg_results <- as.data.frame(mice_kegg_enrichment)

# 找到共同 GO 通路
common_go_pathways <- intersect(human_go_results$ID, mice_go_results$ID)
# 提取前5个共同通路
top_common_go <- head(common_go_pathways, 5)

# 找到不同 GO 通路
unique_human_go_pathways <- setdiff(human_go_results$ID, mice_go_results$ID)
unique_mice_go_pathways <- setdiff(mice_go_results$ID, human_go_results$ID)
# 提取前5个不同通路
top_unique_human_go <- head(unique_human_go_pathways, 5)
top_unique_mice_go <- head(unique_mice_go_pathways, 5)

# 标记通路类型
human_go_results$Type <- ifelse(human_go_results$ID %in% top_common_go, "Common", "Unique to human")
mice_go_results$Type <- ifelse(mice_go_results$ID %in% top_common_go, "Common", "Unique to mice")

# 合并数据框以便于可视化
all_go_pathways <- rbind(
  human_go_results %>% filter(ID %in% c(top_common_go, top_unique_human_go)),
  mice_go_results %>% filter(ID %in% c(top_common_go, top_unique_mice_go)) %>% mutate(Type = ifelse(ID %in% top_common_go, "Common", "Unique to mice"))
)

# 找到共同 KEGG 通路
common_kegg_pathways <- intersect(human_kegg_results$ID, mice_kegg_results$ID)
# 提取前5个共同通路
top_common_kegg <- head(common_kegg_pathways, 5)

# 找到不同 KEGG 通路
unique_human_kegg_pathways <- setdiff(human_kegg_results$ID, mice_kegg_results$ID)
unique_mice_kegg_pathways <- setdiff(mice_kegg_results$ID, human_kegg_results$ID)
# 提取前5个不同通路
top_unique_human_kegg <- head(unique_human_kegg_pathways, 5)
top_unique_mice_kegg <- head(unique_mice_kegg_pathways, 5)

# 标记通路类型
human_kegg_results$Type <- ifelse(human_kegg_results$ID %in% top_common_kegg, "Common", "Unique to human")
mice_kegg_results$Type <- ifelse(mice_kegg_results$ID %in% top_common_kegg, "Common", "Unique to mice")

# 合并数据框以便于可视化
all_kegg_pathways <- rbind(
  human_kegg_results %>% filter(ID %in% c(top_common_kegg, top_unique_human_kegg)),
  mice_kegg_results %>% filter(ID %in% c(top_common_kegg, top_unique_mice_kegg)) %>% mutate(Type = ifelse(ID %in% top_common_kegg, "Common", "Unique to mice"))
)
library(cowplot)
library(ggplot2)

# 创建 GO 通路的条形图
go_plot <- ggplot(all_go_pathways, aes(x = -log10(p.adjust), y = reorder(Description, p.adjust), fill = Type)) +
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
  labs(title = "Top 5 Common and Unique GO Pathways", x = "GO Pathway", y = "-log10(p.adjust)", fill = "Type") +
  scale_fill_manual(values = c("Common" = "darkblue", "Unique to human" = "darkred", "Unique to mice" = "darkgreen"))

# 创建 KEGG 通路的条形图
kegg_plot <- ggplot(all_kegg_pathways, aes(x = -log10(p.adjust), y = reorder(Description, p.adjust), fill = Type)) +
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
  labs(title = "Top 5 Common and Unique KEGG Pathways", x = "KEGG Pathway", y = "-log10(p.adjust)", fill = "Type") +
  scale_fill_manual(values = c("Common" = "darkblue", "Unique to human" = "darkred", "Unique to mice" = "darkgreen"))

# 拼接图像并保存
combined_plot <- plot_grid(go_plot, kegg_plot, nrow = 2, rel_heights = c(1, 1.2))  # 调整图形的相对高度比例

# 保存图像，调整图形尺寸
ggsave(paste0(dirname, "/ARTO_combined_pathways_plot.tiff"), combined_plot, dpi = 300, width = 12, height = 12)  # 减小图形尺寸
