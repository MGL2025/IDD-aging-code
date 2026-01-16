library(ggvenn)
library(RColorBrewer)
library(viridis)
library(ComplexUpset)
library(ggplot2)
# 1. 设置工作目录
setwd("I:\\2025教程新n\\14.Q1区7.8分，MR+单细胞\\09.差异分析和研究方向基因取交集")

# 2. 创建输出文件夹
output_folder <- "output_folder"
if (!dir.exists(output_folder)) dir.create(output_folder)

# 3. 获取所有TXT文件并存储在列表中
txt_files <- list.files(pattern = "\\.txt$", full.names = TRUE)

# 4. 读取所有文件并存入列表
gene_list <- list()

for (file in txt_files) {
  # 获取文件名（不带后缀）
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # 读取文件内容
  rt <- read.table(file, header = FALSE, sep = "\t", check.names = FALSE)
  
  # 将基因数据添加到gene_list中
  gene_list[[file_name]] <- as.vector(rt[, 1])
  
  # 打印日志
  cat(file_name, "基因数：", length(gene_list[[file_name]]), "\n")
}

# 5. 利用 viridis 调色板生成颜色
nColors <- length(gene_list)
myColors <- viridis(nColors)

# 6. 绘制Venn图并保存为PDF
pdf(file = file.path(output_folder, "venn.pdf"), width = 10, height = 10)  # 增大图表尺寸

ggvenn(gene_list, show_percentage = TRUE,
       stroke_color = "white", stroke_size = 1.5,
       fill_color = myColors,       # 使用viridis调色板
       set_name_color = "black",      # 设置集合名称为黑色
       set_name_size = 8,             # 调整集合名称的字体大小
       text_size = 6,                 # 调整文本大小
       text_color = "black")          # 设置文本颜色为黑色

dev.off()

cat("Venn图保存至：", file.path(output_folder, "venn.pdf"), "\n")

# 7. 获取所有文件的全局交集基因并保存（CSV格式）
intersect_genes <- Reduce(intersect, gene_list)
cat("全局交集基因数：", length(intersect_genes), "\n")

# 将全局交集转为数据框保存
global_intersect_df <- data.frame(Gene = intersect_genes, stringsAsFactors = FALSE)
write.csv(global_intersect_df,
          file = file.path(output_folder, "IntersectionGenes.csv"),
          row.names = FALSE,
          quote = FALSE)

cat("全局交集基因已保存至CSV文件：", file.path(output_folder, "IntersectionGenes.csv"), "\n")

# 8. 生成每两个TXT文件之间的交集基因列表，并保存为CSV表格
pairwise_intersections <- data.frame(
  File1 = character(),
  File2 = character(),
  Intersection_Count = integer(),
  Intersection_Genes = character(),
  stringsAsFactors = FALSE
)

files <- names(gene_list)
for (i in 1:(length(files) - 1)) {
  for (j in (i + 1):length(files)) {
    file1 <- files[i]
    file2 <- files[j]
    
    common_genes <- intersect(gene_list[[file1]], gene_list[[file2]])
    
    pairwise_intersections <- rbind(pairwise_intersections,
                                    data.frame(
                                      File1 = file1,
                                      File2 = file2,
                                      Intersection_Count = length(common_genes),
                                      Intersection_Genes = paste(common_genes, collapse = ";"),
                                      stringsAsFactors = FALSE
                                    ))
  }
}

write.csv(pairwise_intersections,
          file = file.path(output_folder, "PairwiseIntersectionGenes.csv"),
          row.names = FALSE,
          quote = FALSE)

cat("每两个TXT文件间的交集基因表格已保存至CSV文件：", file.path(output_folder, "PairwiseIntersectionGenes.csv"), "\n")

# 假设 gene_list 是一个命名列表，每个元素包含一个基因向量
all_genes <- unique(unlist(gene_list))
membership_df <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)

# 对于 gene_list 中的每个集合，增加一列指示该基因是否存在（1 表示存在，0 表示不存在）
for (set_name in names(gene_list)) {
  membership_df[[set_name]] <- as.integer(membership_df$Gene %in% gene_list[[set_name]])
}

# 查看前几行数据检查是否正确构建
head(membership_df)
library(ComplexUpset)
library(ggplot2)

upset_plot <- ComplexUpset::upset(
  membership_df,
  intersect = names(gene_list),  # 使用 gene_list 的列作为集合
  name = "Gene",                 # 主列名称
  sort_intersections_by = "degree",  # 按交集数量排序
  base_annotations = list(
    "Intersection Size" = ComplexUpset::intersection_size(
      text = element_text(size = 12)  # 设置文字大小
    )
  )
)

# 显示图形
print(upset_plot)
# 假设 output_folder 已经定义为保存文件的目标文件夹
ggsave(
  filename = file.path(output_folder, "ComplexUpSetPlot.pdf"),
  plot = upset_plot,
  width = 8,
  height = 6
)
