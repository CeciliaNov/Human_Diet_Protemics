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

# 创建输出目录
dirname <- "20240728enrichplot"
dir.create(dirname)

# 读取数据
dep <- read.csv("diet_dep.csv")

# 计算不同饮食中的差异蛋白质数据
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
      pvalue = min(c_across(starts_with("pval")), na.rm = TRUE)
    ) %>%
    select(SYMBOL, logFC, pvalue) %>%
    filter(abs(logFC) > 0.58)
  return(diff_proteins)
}

diff_proteins_high_glucose <- calculate_diff_proteins(dep, "A")
diff_proteins_ketogenic <- calculate_diff_proteins(dep, "B")
diff_proteins_hfhc <- calculate_diff_proteins(dep, "C")
diff_proteins_LFLC <- calculate_diff_proteins(dep, "E")
diff_proteins_normal <- calculate_diff_proteins(dep, "D")

# 进行富集分析
perform_enrichment <- function(genes) {
  enrich_result <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  return(enrich_result)
}

enrich_high_glucose <- perform_enrichment(diff_proteins_high_glucose$SYMBOL)
enrich_ketogenic <- perform_enrichment(diff_proteins_ketogenic$SYMBOL)
enrich_hfhc <- perform_enrichment(diff_proteins_hfhc$SYMBOL)
enrich_LFLC <- perform_enrichment(diff_proteins_LFLC$SYMBOL)
enrich_normal <- perform_enrichment(diff_proteins_normal$SYMBOL)

# 选择显著通路
extract_significant_pathways <- function(enrich_result, diet_name) {
  significant_pathways <- enrich_result@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(Diet = diet_name)
  return(significant_pathways)
}

significant_high_glucose <- extract_significant_pathways(enrich_high_glucose, "High-glucose")
significant_ketogenic <- extract_significant_pathways(enrich_ketogenic, "Ketogenic")
significant_hfhc <- extract_significant_pathways(enrich_hfhc, "HFHC")
significant_LFLC <- extract_significant_pathways(enrich_LFLC, "LFLC")
significant_normal <- extract_significant_pathways(enrich_normal, "Normal")

# 合并结果
combined_results <- bind_rows(significant_high_glucose, significant_ketogenic, significant_hfhc, significant_LFLC, significant_normal)

# 共同变化的通路
common_pathways <- Reduce(intersect, list(significant_high_glucose$Description, significant_ketogenic$Description, significant_hfhc$Description, significant_LFLC$Description, significant_normal$Description))

# 不同变化的通路
unique_high_glucose <- setdiff(significant_high_glucose$Description, c(significant_ketogenic$Description, significant_hfhc$Description, significant_LFLC$Description, significant_normal$Description))
unique_ketogenic <- setdiff(significant_ketogenic$Description, c(significant_high_glucose$Description, significant_hfhc$Description, significant_LFLC$Description, significant_normal$Description))
unique_hfhc <- setdiff(significant_hfhc$Description, c(significant_high_glucose$Description, significant_ketogenic$Description, significant_LFLC$Description, significant_normal$Description))
unique_LFLC <- setdiff(significant_LFLC$Description, c(significant_high_glucose$Description, significant_ketogenic$Description, significant_hfhc$Description, significant_normal$Description))
unique_normal <- setdiff(significant_normal$Description, c(significant_high_glucose$Description, significant_ketogenic$Description, significant_hfhc$Description, significant_LFLC$Description))

# 提取前5个共同变化的通路
top_common_pathways <- combined_results %>%
  filter(Description %in% common_pathways) %>%
  arrange(pvalue) %>%
  distinct(Description) %>%
  head(5) %>%
  pull(Description)

# 提取共同变化的通路数据并按 pvalue 排序选择前5种
common_results <- combined_results %>%
  filter(Description %in% top_common_pathways) %>%
  mutate(Category = "Common")

# 提取不同变化的通路数据并按 pvalue 排序选择前5种
unique_high_glucose_results <- combined_results %>%
  filter(Description %in% unique_high_glucose) %>%
  arrange(pvalue) %>%
  head(5) %>%
  mutate(Category = "Unique to High-glucose")

unique_ketogenic_results <- combined_results %>%
  filter(Description %in% unique_ketogenic) %>%
  arrange(pvalue) %>%
  head(5) %>%
  mutate(Category = "Unique to Ketogenic")

unique_hfhc_results <- combined_results %>%
  filter(Description %in% unique_hfhc) %>%
  arrange(pvalue) %>%
  head(5) %>%
  mutate(Category = "Unique to HFHC")

unique_LFLC_results <- combined_results %>%
  filter(Description %in% unique_LFLC) %>%
  arrange(pvalue) %>%
  head(5) %>%
  mutate(Category = "Unique to LFLC")

unique_normal_results <- combined_results %>%
  filter(Description %in% unique_normal) %>%
  arrange(pvalue) %>%
  head(5) %>%
  mutate(Category = "Unique to Normal")

# 合并共同变化和不同变化的结果
plot_data <- bind_rows(common_results, unique_high_glucose_results, unique_ketogenic_results, unique_hfhc_results, unique_LFLC_results, unique_normal_results)

# 创建一个包含所有数据的ggplot
plot <- ggplot(plot_data, aes(x = Diet, y = reorder(Description, pvalue), size = Count, color = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  labs(
    title = "Enriched GO Terms for Different Diets",
    x = "Diet",
    y = "GO Term",
    color = "-log10(p.adjust)",
    size = "Gene Count"
  ) +
  facet_wrap(~Category, scales = "free_y", ncol = 1)

# 保存图像
ggsave(paste0(dirname,"/diet_comparison.tiff"), plot, width = 10, height = 15)

# 绘制共同通路的 cnetplot
common_result <- enrich_high_glucose
common_result@result <- enrich_high_glucose@result %>%
  filter(Description %in% common_pathways) %>%
  arrange(pvalue) %>%
  head(5)

cnet_common <- cnetplot(common_result, showCategory = 5) + 
  ggtitle("Common Pathways") + 
  theme_minimal()

# 绘制 High-glucose 特有通路的 cnetplot
unique_high_glucose_result <- enrich_high_glucose
unique_high_glucose_result@result <- enrich_high_glucose@result %>% 
  filter(Description %in% unique_high_glucose) %>% 
  arrange(pvalue) %>% 
  head(5)

cnet_high_glucose <- cnetplot(unique_high_glucose_result, showCategory = 5) + 
  ggtitle("Unique to High-glucose") + 
  theme_minimal()

# 绘制 Ketogenic 特有通路的 cnetplot
unique_ketogenic_result <- enrich_ketogenic
unique_ketogenic_result@result <- enrich_ketogenic@result %>% 
  filter(Description %in% unique_ketogenic) %>% 
  arrange(pvalue) %>% 
  head(5)

cnet_ketogenic <- cnetplot(unique_ketogenic_result, showCategory = 5) + 
  ggtitle("Unique to Ketogenic") + 
  theme_minimal()

# 绘制 HFHC 特有通路的 cnetplot
unique_hfhc_result <- enrich_hfhc
unique_hfhc_result@result <- enrich_hfhc@result %>% 
  filter(Description %in% unique_hfhc) %>% 
  arrange(pvalue) %>% 
  head(5)

cnet_hfhc <- cnetplot(unique_hfhc_result, showCategory = 5) + 
  ggtitle("Unique to HFHC") + 
  theme_minimal()

# 绘制 LFLC 特有通路的 cnetplot
unique_LFLC_result <- enrich_LFLC
unique_LFLC_result@result <- enrich_LFLC@result %>% 
  filter(Description %in% unique_LFLC) %>% 
  arrange(pvalue) %>% 
  head(5)

cnet_LFLC <- cnetplot(unique_LFLC_result, showCategory = 5) + 
  ggtitle("Unique to LFLC") + 
  theme_minimal()

# 合并图形
# Example of setting explicit dimensions in ggplot2
combined_cnetplot <- grid.arrange(cnet_common, cnet_high_glucose, cnet_ketogenic, cnet_hfhc, cnet_LFLC, nrow=1)


# 保存图像
ggsave(paste0(dirname,"/diet_cnetplot.tiff"), plot = combined_cnetplot, width = 35, height = 8, dpi = 300)
