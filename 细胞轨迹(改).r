# 设置工作目录
workDir <- "D:/2.椎间盘衰老(1)/Q1区纯生信7.8分/33.单细胞细胞轨迹分析"
setwd(workDir)
message("工作目录设置为：", getwd())

# 感兴趣的基因列表
showGenes <- c("SF3A3", "GSTZ1")

# 加载数据
load(file = "filtered_seurat_2.rda")

# 检查数据
print(filtered_seurat)
table(filtered_seurat$CellType)

# 安装并加载必要的包
if (!requireNamespace("monocle3", quietly = TRUE)) {
  install.packages("monocle3")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)

# 提取Chondrocytes细胞
chondrocytes_seurat <- subset(filtered_seurat, subset = cell_type == "Chondrocytes")
print(paste("提取到", ncol(chondrocytes_seurat), "个Chondrocytes细胞"))

# 将Seurat对象转换为Monocle3对象
# 使用GetAssayData获取counts数据，确保使用正确的slot
expression_matrix <- GetAssayData(chondrocytes_seurat, assay = "RNA", slot = "counts")

# 创建细胞metadata
cell_metadata <- chondrocytes_seurat@meta.data

# 创建基因metadata
gene_annotation <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

# 创建CellDataSet对象
cds <- new_cell_data_set(
  expression_data = expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

# 检查Type列的内容
print("Type列的内容:")
print(table(colData(cds)$Type))

# 使用Type列作为分组信息
colData(cds)$group <- colData(cds)$Type
print("使用Type列作为分组信息:")
print(table(colData(cds)$group))

# 现在继续轨迹分析的其余部分
# 预处理
cds <- preprocess_cds(cds, 
                      num_dim = 14,  # 使用与Seurat相同的PC数量
                      method = "PCA",
                      norm_method = "log")

print("预处理完成")

# 降维和聚类
cds <- reduce_dimension(cds, 
                        preprocess_method = "PCA",
                        reduction_method = "UMAP")

print("降维完成")

cds <- cluster_cells(cds, 
                     resolution = 1e-5,
                     reduction_method = "UMAP")

print("聚类完成")

# 轨迹分析
# cds <- learn_graph(cds, 
#                    use_partition = TRUE,
#                    close_loop = TRUE)

cds <- learn_graph(cds, 
                   use_partition = FALSE,#强制一个分区（整个数据一个轨迹）
                   close_loop = TRUE,  # 关闭闭环，简化轨迹
                   learn_graph_control = list(
                     minimal_branch_len = 8,  # 增加最小分支长度，减少小分支
                     ncenter = 50  # 减少中心点数量
                   ))
print("轨迹学习完成")
################################################
library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
load(file = "chondrocytes_monocle_cds.rda")
# 感兴趣的基因列表
showGenes <- c("SF3A3", "GSTZ1")
# 绘制基础轨迹图
p1 <- plot_cells(cds,
                 color_cells_by = "cluster",
                 label_groups_by_cluster = TRUE,
                 label_leaves = TRUE,
                 label_branch_points = TRUE,
                 graph_label_size = 3,
                 cell_size = 0.5) +
  ggtitle("Chondrocytes Trajectory Analysis") +
  theme(plot.title = element_text(hjust = 0.5))

print(p1)

# 保存基础轨迹图
ggsave("chondrocytes_trajectory_basic.pdf", p1, width = 5, height = 4)
ggsave("chondrocytes_trajectory_basic.png", p1, width = 5, height = 4, dpi = 300)

# 按分组（CON vs IDD）着色的轨迹图
p2 <- plot_cells(cds,
                 color_cells_by = "group",
                 label_cell_groups = FALSE,  # 不显示分组标签
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 graph_label_size = 3,
                 cell_size = 0.5) +
  scale_color_manual(
    values = c("CON" = "#4DBBD5", "IDD" = "#F39B7F"),
    name = "Group"  # 设置图例标题
  ) +
  ggtitle("Chondrocytes Trajectory by Group (CON vs IDD)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"  # 确保图例显示在右侧
  )

print(p2)

# 保存分组轨迹图
ggsave("chondrocytes_trajectory_by_group.pdf", p2, width = 5, height = 4)
ggsave("chondrocytes_trajectory_by_group.png", p2, width = 5, height = 4, dpi = 300)

# 3. 计算伪时间
# 首先需要选择轨迹的根节点
# 我们可以根据先验知识选择根节点，或者让Monocle自动选择
# 这里我们让Monocle自动选择根节点
cds <- order_cells(cds)

print("伪时间计算完成")

# 4. 伪时间轨迹图（按伪时间颜色渐变）
p3 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 graph_label_size = 3,
                 cell_size = 0.5,
                 show_trajectory_graph = TRUE) +
  scale_color_viridis_c() +
  ggtitle("Chondrocytes Trajectory by Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5))

print(p3)
ggsave("chondrocytes_trajectory_pseudotime.pdf", p3, width = 5, height = 4)
ggsave("chondrocytes_trajectory_pseudotime.png", p3, width = 5, height = 4, dpi = 300)

# 5. State轨迹图
###第一张图用来看state
p4 <- plot_cells(cds,
                 color_cells_by = "partition",
                 label_cell_groups = TRUE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 graph_label_size = 3,
                 cell_size = 0.5,
                 group_label_size = 4) +
  ggtitle("Chondrocytes Trajectory by State") +
  theme(plot.title = element_text(hjust = 0.5))

print(p4)
ggsave("chondrocytes_trajectory_state.pdf", p4, width = 10, height = 8)
ggsave("chondrocytes_trajectory_state.png", p4, width = 10, height = 8, dpi = 300)

# 方法2：使用自定义颜色和主题
p4_simple <- plot_cells(cds,
                        color_cells_by = "partition",
                        label_cell_groups = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = FALSE,
                        graph_label_size = 3,
                        cell_size = 0.5) +
  scale_color_manual(values = c("#F8766D", "#B79F00", "#00BFC4")) +
  labs(color = "State") +  # 只修改图例标题
  ggtitle("Chondrocytes Trajectory by State") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "right")

print(p4_simple)
ggsave("chondrocytes_trajectory_state_no_labels.pdf", p4_simple, width = 5, height = 4)
ggsave("chondrocytes_trajectory_state_no_labels.png", p4_simple, width = 5, height = 4, dpi = 300)

# 检查靶基因是否在数据中
available_genes <- showGenes[showGenes %in% rownames(cds)]
print(paste("可用的靶基因:", paste(available_genes, collapse = ", ")))

# 6. 靶基因在轨迹上的表达图（按表达量渐变颜色）
for (gene in available_genes) {
  p_gene <- plot_cells(cds,
                       genes = gene,
                       label_cell_groups = FALSE,
                       label_leaves = FALSE,
                       label_branch_points = FALSE,
                       cell_size = 0.5,
                       show_trajectory_graph = TRUE) +
    scale_color_viridis_c() +
    ggtitle(paste("Expression of", gene, "on Chondrocytes Trajectory")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p_gene)
  ggsave(paste0("chondrocytes_trajectory_", gene, ".pdf"), p_gene, width = 5, height = 4)
  ggsave(paste0("chondrocytes_trajectory_", gene, ".png"), p_gene, width = 5, height = 4, dpi = 300)
}

# 7. 靶基因随Pseudotime变化表达量的点线图
for (gene in available_genes) {
  # 提取基因表达数据
  gene_expr <- exprs(cds)[which(rownames(cds) == gene), ]
  
  # 提取伪时间
  pseudotime <- pseudotime(cds)
  
  # 创建数据框
  plot_data <- data.frame(
    pseudotime = pseudotime,
    expression = as.numeric(gene_expr),
    group = colData(cds)$group
  )
  
  # 移除NA值
  plot_data <- plot_data[!is.na(plot_data$pseudotime), ]
  
  # 绘制点线图
  p_line <- ggplot(plot_data, aes(x = pseudotime, y = expression)) +
    geom_point(aes(color = group), alpha = 0.6, size = 0.8) +
    geom_smooth(method = "loess", color = "black", se = TRUE) +
    scale_color_manual(values = c("CON" = "#4DBBD5", "IDD" = "#F39B7F")) +
    labs(title = paste("Expression of", gene, "along Pseudotime"),
         x = "Pseudotime",
         y = "Expression Level") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p_line)
  ggsave(paste0("chondrocytes_", gene, "_pseudotime_expression.pdf"), p_line, width = 6, height = 3)
  ggsave(paste0("chondrocytes_", gene, "_pseudotime_expression.png"), p_line, width = 6, height = 3, dpi = 300)
  
  # 可选：按分组分别绘制平滑曲线
  p_line_by_group <- ggplot(plot_data, aes(x = pseudotime, y = expression, color = group)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "loess", se = TRUE) +
    scale_color_manual(values = c("CON" = "#4195C1", "IDD" = "#FFBC90")) +
    labs(title = paste("Expression of", gene, "along Pseudotime by Group"),
         x = "Pseudotime",
         y = "Expression Level") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p_line_by_group)
  ggsave(paste0("chondrocytes_", gene, "_pseudotime_expression_by_group.pdf"), p_line_by_group, width = 6, height = 3)
  ggsave(paste0("chondrocytes_", gene, "_pseudotime_expression_by_group.png"), p_line_by_group, width = 6, height = 3, dpi = 300)
}

# 8. 组合基因图
if(length(available_genes) > 1) {
  p_combined <- plot_cells(cds,
                           genes = available_genes,
                           label_cell_groups = FALSE,
                           label_leaves = FALSE,
                           label_branch_points = FALSE,
                           cell_size = 0.5) +
    ggtitle("Combined Expression of Target Genes on Chondrocytes Trajectory") +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p_combined)
  ggsave("chondrocytes_trajectory_combined_genes.pdf", p_combined, width = 12, height = 8)
  ggsave("chondrocytes_trajectory_combined_genes.png", p_combined, width = 12, height = 8, dpi = 300)
}

# 保存Monocle对象
save(cds, file = "chondrocytes_monocle_cds.rda")
load(file = "chondrocytes_monocle_cds.rda")
