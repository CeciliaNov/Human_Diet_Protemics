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

# 加载必要的包
library(dplyr)
library(readr)  # 如果使用read_csv函数读取CSV文件

# 定义文件路径及对应的饮食标签
file_paths <- list(
  "Result_gee_lgFC/A_gee_lgFC0.26_sel.csv" = "A",
  "Result_gee_lgFC/B_gee_lgFC0.26_sel.csv" = "B",
  "Result_gee_lgFC/C_gee_lgFC0.26_sel.csv" = "C",
  "Result_gee_lgFC/D_gee_lgFC0.26_sel.csv" = "D",
  "Result_gee_lgFC/E_gee_lgFC0.26_sel.csv" = "E"
)

# 读取每个文件并添加饮食列
data_list <- lapply(names(file_paths), function(file) {
  diet_label <- file_paths[[file]]
  df <- read_csv(file) %>%  # 或使用 read.csv(file) 取决于你的文件格式
    mutate(diet = diet_label)
  return(df)
})

# 合并所有数据
dep <- bind_rows(data_list)

# 保存合并结果到一个新文件
write_csv(dep, "dep_diet.csv")  # 或使用 write.csv(combined_data, "combined_data.csv", row.names = FALSE)




# 计算不同饮食中的差异蛋白质数据
library(readr)
dep <- read_csv("dep_diet.csv")
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

diff_proteins_high_glucose <- calculate_diff_proteins(dep, "A")
diff_proteins_ketogenic <- calculate_diff_proteins(dep, "B")
diff_proteins_hfhc <- calculate_diff_proteins(dep, "C")
diff_proteins_LFLC <- calculate_diff_proteins(dep, "E")
diff_proteins_normal <- calculate_diff_proteins(dep, "D")




library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 定义一个函数来转换 UNIPROT 到 SYMBOL
convert_uniprot_to_symbol <- function(data) {
  uniprot_ids <- unique(data$UNIPROT)
  uniprot_to_symbol <- bitr(uniprot_ids, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  
  # 合并基因符号到数据框
  data <- data %>%
    left_join(uniprot_to_symbol, by = c("UNIPROT" = "UNIPROT")) %>%
    rename(SYMBOL = SYMBOL)
  
  return(data)
}

# 批量转换 UNIPROT 到 SYMBOL
diff_proteins_high_glucose <- convert_uniprot_to_symbol(diff_proteins_high_glucose)
diff_proteins_ketogenic <- convert_uniprot_to_symbol(diff_proteins_ketogenic)
diff_proteins_hfhc <- convert_uniprot_to_symbol(diff_proteins_hfhc)
diff_proteins_LFLC <- convert_uniprot_to_symbol(diff_proteins_LFLC)
diff_proteins_normal <- convert_uniprot_to_symbol(diff_proteins_normal)



# 合并所有饮食的差异蛋白质数据
all_diff_proteins <- bind_rows(
  mutate(diff_proteins_high_glucose, diet = "High Glucose"),
  mutate(diff_proteins_ketogenic, diet = "Ketogenic"),
  mutate(diff_proteins_hfhc, diet = "HFHC"),
  mutate(diff_proteins_LFLC, diet = "LFLC"),
  mutate(diff_proteins_normal, diet = "Normal")
)

# 可视化
library(ggplot2)

# Define custom colors
custom_colors <- c("High Glucose" = "darkblue", "HFHC" = "#ff7f0e", "LFLC" = "#2ca02c", "Ketogenic" = "#d62728")

# Create the plot with custom colors
p <- ggplot(all_diff_proteins, aes(x = SYMBOL, y = logFC, color = diet)) +
  geom_point(aes(size = -log10(pvalue))) +
  theme_minimal(base_size = 14) +  # Increase base font size
  labs(
    title = "Differential Proteins log2FC Across Different Diets",
    x = "Symbol",
    y = "log2FC",
    color = "Diet",
    size = "-log10(pvalue)"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase x-axis text size
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    plot.title = element_text(size = 16, face = "bold"),  # Increase title size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12)  # Increase legend text size
  )

print(p)


# Save the image with high resolution
ggsave(paste0(dirname, "/urine.tiff"), plot = p, width = 8, height = 6, dpi = 300)


###人数据####################################################
# 读取数据
dep <- read.csv("D:/Cycy1/R_Proj/Diet/diet_R/diet_dep.csv")

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


# 合并所有饮食的差异蛋白质数据
plasma_data <- bind_rows(
  mutate(diff_proteins_high_glucose, diet = "High Glucose"),
  mutate(diff_proteins_ketogenic, diet = "Ketogenic"),
  mutate(diff_proteins_hfhc, diet = "HFHC"),
  mutate(diff_proteins_LFLC, diet = "LFLC"),
  mutate(diff_proteins_normal, diet = "Normal")
)


# 过滤出基因符号一致的数据
filtered_urine_data <- all_diff_proteins %>%
  filter(SYMBOL %in% plasma_data$SYMBOL)

filtered_plasma_data <- plasma_data %>%
  filter(SYMBOL %in% filtered_urine_data$SYMBOL)

# 查看过滤后的数据
print(filtered_urine_data)
print(filtered_plasma_data)

# 为了可视化，将过滤后的尿液和等离子数据合并
comparison_data <- filtered_urine_data %>%
  select(SYMBOL, logFC.urine = logFC, pvalue.urine = pvalue, diet.urine = diet) %>%
  inner_join(
    filtered_plasma_data %>% select(SYMBOL, logFC.plasma = logFC, pvalue.plasma = pvalue, diet.plasma = diet),
    by = "SYMBOL"
  )

# 将尿液和等离子的数据分别作为单独的行添加到一个新的data frame中
urine_plot_data <- comparison_data %>%
  select(SYMBOL, logFC = logFC.urine, pvalue = pvalue.urine, diet = diet.urine) %>%
  mutate(type = "urine")

plasma_plot_data <- comparison_data %>%
  select(SYMBOL, logFC = logFC.plasma, pvalue = pvalue.plasma, diet = diet.plasma) %>%
  mutate(type = "plasma")

plot_data <- bind_rows(urine_plot_data, plasma_plot_data)

# 可视化尿液和等离子之间的logFC相关性，使用不同的点形状来表示不同的类型
# Define custom colors
custom_colors <- c("High Glucose" = "darkblue", "HFHC" = "#ff7f0e", "LFLC" = "#2ca02c", "Ketogenic" = "#d62728","Normal"="lightblue")
p2 <- ggplot(plot_data, aes(x = SYMBOL, y = logFC, color = diet, shape = type)) +
  geom_point(aes(size = -log10(pvalue))) +
  scale_shape_manual(values = c(16, 1)) +# Use solid and hollow points
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme_minimal(base_size = 14) +  # Increase base font size
  labs(
    title = "Comparison Between Urine and Plasma Data",
    x = "Gene Symbol",
    y = "log2FC",
    color = "Diet",
    shape = "Type",
    size = "-log10(pvalue)"
  ) +
  guides(label = "none") +  # Remove labels from legend
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase x-axis text size
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    plot.title = element_text(size = 16, face = "bold"),  # Increase title size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12)  # Increase legend text size
  )

print(p2)

# Save the image with high resolution
ggsave(paste0(dirname, "/urine_plasma1.tiff"), plot = p2, width = 6, height = 6, dpi = 300)




# 进行富集分析



# 计算不同饮食中的差异蛋白质数据
library(readr)
dep <- read_csv("dep_diet.csv")
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

diff_proteins_high_glucose <- calculate_diff_proteins(dep, "A")
diff_proteins_ketogenic <- calculate_diff_proteins(dep, "B")
diff_proteins_hfhc <- calculate_diff_proteins(dep, "C")
diff_proteins_LFLC <- calculate_diff_proteins(dep, "E")
diff_proteins_normal <- calculate_diff_proteins(dep, "D")


library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 定义一个函数来转换 UNIPROT 到 SYMBOL
convert_uniprot_to_symbol <- function(data) {
  uniprot_ids <- unique(data$UNIPROT)
  uniprot_to_symbol <- bitr(uniprot_ids, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  
  # 合并基因符号到数据框
  data <- data %>%
    left_join(uniprot_to_symbol, by = c("UNIPROT" = "UNIPROT")) %>%
    rename(SYMBOL = SYMBOL)
  
  return(data)
}

# 批量转换 UNIPROT 到 SYMBOL
diff_proteins_high_glucose <- convert_uniprot_to_symbol(diff_proteins_high_glucose)
diff_proteins_ketogenic <- convert_uniprot_to_symbol(diff_proteins_ketogenic)
diff_proteins_hfhc <- convert_uniprot_to_symbol(diff_proteins_hfhc)
diff_proteins_LFLC <- convert_uniprot_to_symbol(diff_proteins_LFLC)
diff_proteins_normal <- convert_uniprot_to_symbol(diff_proteins_normal)

#################################################
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


# 合并共同变化和不同变化的结果
plot_data <- bind_rows(common_results, unique_high_glucose_results, unique_ketogenic_results, unique_hfhc_results, unique_LFLC_results)

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
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# 绘制 Ketogenic 特有通路的 cnetplot
unique_ketogenic_result <- enrich_ketogenic
unique_ketogenic_result@result <- enrich_ketogenic@result %>% 
  filter(Description %in% unique_ketogenic) %>% 
  arrange(pvalue) %>% 
  head(5)

cnet_ketogenic <- cnetplot(unique_ketogenic_result, showCategory = 5) + 
  ggtitle("Unique to Ketogenic") + 
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# 绘制 HFHC 特有通路的 cnetplot
unique_hfhc_result <- enrich_hfhc
unique_hfhc_result@result <- enrich_hfhc@result %>% 
  filter(Description %in% unique_hfhc) %>% 
  arrange(pvalue) %>% 
  head(5)

cnet_hfhc <- cnetplot(unique_hfhc_result, showCategory = 5) + 
  ggtitle("Unique to HFHC") + 
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# 绘制 LFLC 特有通路的 cnetplot
unique_LFLC_result <- enrich_LFLC
unique_LFLC_result@result <- enrich_LFLC@result %>% 
  filter(Description %in% unique_LFLC) %>% 
  arrange(pvalue) %>% 
  head(5)

cnet_LFLC <- cnetplot(unique_LFLC_result, showCategory = 5) + 
  ggtitle("Unique to LFLC") + 
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Combine the plots
library(cowplot)
combined_cnetplot <- plot_grid(cnet_high_glucose, cnet_ketogenic, cnet_hfhc, cnet_LFLC, nrow = 1)

# Save the combined plot with appropriate dimensions and high resolution
ggsave(paste0(dirname, "/diet_cnetplot.tiff"), plot = combined_cnetplot, width = 30, height = 10, dpi = 300)
c