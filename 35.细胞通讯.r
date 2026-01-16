#                        单细胞RNA测序数据分析流程                              #
#                    Single Cell RNA-seq Analysis Pipeline                    #
################################################################################

# ==============================================================================
# 01. R包加载
# ==============================================================================

# 核心单细胞分析包
library(Seurat)        # 单细胞数据分析核心包
library(limma)         # 差异表达分析
library(SingleR)       # 单细胞类型注释
library(celldex)       # 细胞类型注释数据集
library(monocle)       # 单细胞轨迹分析
library(monocle3)      # Monocle3轨迹分析
# 数据处理包
library(dplyr)         # 数据处理和管道操作
library(magrittr)      # 管道操作符
library(tidyr)         # 数据整理
library(stringr)       # 字符串处理

# 可视化包
library(ggplot2)       # 基础绘图
library(ggpubr)        # 统计绘图和发表级图形
library(RColorBrewer)  # 调色板
library(viridis)       # 现代配色方案
library(scales)        # 图形标度
library(patchwork)     # 图形组合
library(gridExtra)     # 网格布局

# 统计分析包
library(rstatix)       # 统计检验

# 交互和辅助包
library(DT)            # 交互式表格
library(progress)      # 进度条显示

# 轨迹分析相关包
library(SeuratWrappers) # Seurat封装器
library(reshape2)       # 数据重塑

# 可选的文件输出包
library(openxlsx)      # Excel文件读写

message("所有必需的R包已成功加载完成！")

# ==============================================================================
# 02. 参数设置和工作目录
# ==============================================================================

# 设置工作目录 - 请根据实际路径修改
workDir <- "D:\2.椎间盘衰老(1)\Q1区纯生信7.8分\35.细胞通讯"

# 检查并设置工作目录
if (!file.exists(workDir)) {
  stop("工作目录不存在，请检查路径：", workDir)
} else {
  setwd(workDir)
  message("工作目录设置为：", getwd())
}

# 分析参数设置
analysis_params <- list(
  pcSelect = 14,                      # PCA主成分数选择
  logFCfilter = 1,                    # logFC筛选阈值
  adjPvalFilter = 0.05,               # 调整后p值筛选阈值
  min_cells = 10,                     # 基因至少在多少个细胞中表达
  min_features = 40,                  # 细胞至少表达多少个基因
  cluster_resolution = 0.5,           # 聚类分辨率
  n_variable_features = 2000          # 高变基因数量
)

# 感兴趣的基因列表
showGenes <- c("SF3A3", "GSTZ1")

# ==============================================================================
# 03. 数据读取和Seurat对象创建
# ==============================================================================
# 加载数据
load(file = "filtered_seurat_2.rda")


#基站
# ==============================================================================
# 10. 细胞通讯分析 (Cell Communication Analysis) - 简化版本
# ==============================================================================
# 加载细胞通讯分析所需的包
library(CellChat)
library(patchwork)
library(igraph)

message("开始细胞通讯分析...")

# 细胞通讯分析函数 - 使用clusterAnn.txt注释
perform_cellchat_analysis <- function(seurat_obj) {
  
  # 读取细胞注释文件
  if (!file.exists("clusterAnn.txt")) {
    stop("未找到clusterAnn.txt文件，请确保文件存在于当前工作目录")
  }
  
  cluster_ann <- read.table("clusterAnn.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  message("使用clusterAnn.txt中的细胞注释")
  message(paste("clusterAnn.txt包含", nrow(cluster_ann), "个cluster注释"))
  
  # 获取表达数据
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  
  # 获取元数据
  meta <- seurat_obj@meta.data
  
  # 获取当前Seurat对象的cluster信息
  current_clusters <- as.character(Idents(seurat_obj))
  
  # 检查cluster注释文件的列名
  if ("labels" %in% colnames(cluster_ann)) {
    celltype_col <- "labels"
  } else if ("celltype" %in% colnames(cluster_ann)) {
    celltype_col <- "celltype"
  } else {
    celltype_col <- colnames(cluster_ann)[2]
  }
  
  # 创建cluster到celltype的映射
  cluster_to_celltype <- setNames(cluster_ann[[celltype_col]], cluster_ann$Cluster)
  
  # 为每个细胞分配细胞类型
  cell_types <- cluster_to_celltype[current_clusters]
  
  # 检查是否有未匹配的cluster
  unmatched_clusters <- unique(current_clusters[is.na(cell_types)])
  if (length(unmatched_clusters) > 0) {
    message(paste("未匹配的cluster:", paste(unmatched_clusters, collapse = ", ")))
    cell_types[is.na(cell_types)] <- paste0("Cluster_", current_clusters[is.na(cell_types)])
  }
  
  # 将细胞类型添加到元数据中
  meta$annotated_cell_type <- cell_types
  
  # 根据样本名称创建分组
  meta$group <- ifelse(grepl("Control|Normal|Ctrl|control|normal|ctrl|CON", rownames(meta)), 
                       "CON", "IDD")
  
  # 验证细胞类型分配
  message("细胞类型分配验证:")
  cell_type_summary <- table(meta$annotated_cell_type, meta$group)
  print(cell_type_summary)
  
  # 分别为对照组和实验组创建CellChat对象
  groups <- unique(meta$group)
  cellchat_list <- list()
  
  for (group in groups) {
    message(paste("\n处理", group, "组..."))
    
    # 筛选当前组的细胞
    cells_subset <- rownames(meta)[meta$group == group]
    data_subset <- data.input[, cells_subset]
    meta_subset <- meta[cells_subset, ]
    
    # 创建CellChat对象
    cellchat <- createCellChat(object = data_subset, meta = meta_subset, group.by = "annotated_cell_type")
    
    # 设置配体受体数据库
    CellChatDB <- CellChatDB.human
    cellchat@DB <- CellChatDB
    
    # 预处理表达数据
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # 计算通讯概率
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    # 推断细胞通讯网络
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    message(paste(group, "组分析完成"))
    message(paste("  细胞类型数量:", length(levels(cellchat@idents))))
    message(paste("  细胞总数:", length(cellchat@idents)))
    
    cellchat_list[[group]] <- cellchat
  }
  
  return(cellchat_list)
}

# 执行细胞通讯分析
cellchat_results <- perform_cellchat_analysis(filtered_seurat)
######################################
save(analysis_params, annotation_results, cellchat_results, filtered_seurat, marker_results,
     stats_summary, file = "cellchat_result.rda")
load(file = "cellchat_result.rda")
file = "cellchat_result.rda"
# ==============================================================================
# 11. 生成四个核心图形和统计表格（支持颜色自定义）
# ==============================================================================

# 简化的可视化函数 - 支持颜色自定义
visualize_cellchat_core <- function(cellchat_list, 
                                    output_dir = "CellChat_Results",
                                    color_palette = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", 
                                                      "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"),
                                    vertex_size_max = 15,  # 新增：节点大小控制
                                    edge_width_max = 8) {  # 新增：边宽度控制
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 初始化统计数据框
  communication_stats <- data.frame(
    Group = character(),
    Cell_Types = numeric(),
    Total_Cells = numeric(),
    Active_Interactions = numeric(),
    Total_Strength = numeric(),
    Average_Strength = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 为每个组生成网络图
  for (group_name in names(cellchat_list)) {
    cellchat <- cellchat_list[[group_name]]
    
    message(paste("生成", group_name, "组的细胞通讯网络图..."))
    
    # 检查是否有有效的通讯网络
    if (is.null(cellchat@net$count) || sum(cellchat@net$count) == 0) {
      message(paste("警告:", group_name, "组没有找到有效的细胞通讯"))
      
      # 记录空数据
      communication_stats <- rbind(communication_stats, data.frame(
        Group = group_name,
        Cell_Types = length(levels(cellchat@idents)),
        Total_Cells = length(cellchat@idents),
        Active_Interactions = 0,
        Total_Strength = 0,
        Average_Strength = 0
      ))
      next
    }
    
    # 设置网络布局参数
    groupSize <- as.numeric(table(cellchat@idents))
    n_celltypes <- length(levels(cellchat@idents))
    
    # 颜色设置逻辑
    if (!is.null(color_palette)) {
      # 使用用户提供的颜色
      if (length(color_palette) >= n_celltypes) {
        cell_colors <- color_palette[1:n_celltypes]
        message(paste("✓ 使用自定义颜色:", paste(cell_colors, collapse = ", ")))
      } else {
        warning(paste("自定义颜色数量不足，使用默认颜色。需要", n_celltypes, "种颜色，但只提供了", length(color_palette)))
        cell_colors <- rainbow(n_celltypes)
      }
    } else {
      # 使用默认颜色
      cell_colors <- rainbow(n_celltypes)
      message(paste("使用默认彩虹颜色，包含", n_celltypes, "种颜色"))
    }
    
    # 图1: 细胞通讯数量网络图
    tryCatch({
      pdf(file = file.path(output_dir, paste0(group_name, "_Communication_Count_Network.pdf")), 
          width = 5, height = 5)
      
      par(mfrow = c(1,1), xpd = TRUE, mar = c(2,2,4,2))
      netVisual_circle(cellchat@net$count, 
                       vertex.weight = groupSize, 
                       weight.scale = TRUE, 
                       label.edge = FALSE,
                       title.name = paste("Number of Interactions -", group_name),
                       vertex.label.cex = 1.0,
                       edge.width.max = edge_width_max,
                       vertex.size.max = vertex_size_max,
                       color.use = cell_colors,  # 应用自定义颜色
                       sources.use = NULL, 
                       targets.use = NULL)
      
      dev.off()
      message(paste("✓ 完成", group_name, "通讯数量网络图"))
    }, error = function(e) {
      message(paste("✗ 生成", group_name, "通讯数量网络图失败:", e$message))
      if (dev.cur() != 1) dev.off()
    })
    
    # 图2: 细胞通讯强度网络图  
    tryCatch({
      pdf(file = file.path(output_dir, paste0(group_name, "_Communication_Strength_Network.pdf")), 
          width = 5, height = 5)
      
      par(mfrow = c(1,1), xpd = TRUE, mar = c(2,2,4,2))
      netVisual_circle(cellchat@net$weight, 
                       vertex.weight = groupSize, 
                       weight.scale = TRUE, 
                       label.edge = FALSE,
                       title.name = paste("Interaction Strength -", group_name),
                       vertex.label.cex = 1.0,
                       edge.width.max = edge_width_max,
                       vertex.size.max = vertex_size_max,
                       color.use = cell_colors,  # 应用自定义颜色
                       sources.use = NULL, 
                       targets.use = NULL)
      
      dev.off()
      message(paste("✓ 完成", group_name, "通讯强度网络图"))
    }, error = function(e) {
      message(paste("✗ 生成", group_name, "通讯强度网络图失败:", e$message))
      if (dev.cur() != 1) dev.off()
    })
    
    # 统计信息收集
    total_interactions <- sum(cellchat@net$count > 0)
    total_strength <- sum(cellchat@net$weight)
    average_strength <- ifelse(total_interactions > 0, total_strength / total_interactions, 0)
    
    # 添加到统计数据框
    communication_stats <- rbind(communication_stats, data.frame(
      Group = group_name,
      Cell_Types = length(levels(cellchat@idents)),
      Total_Cells = length(cellchat@idents),
      Active_Interactions = total_interactions,
      Total_Strength = round(total_strength, 3),
      Average_Strength = round(average_strength, 3)
    ))
  }
  
  # 保存统计表格
  write.csv(communication_stats, file = file.path(output_dir, "Communication_Statistics.csv"), 
            row.names = FALSE)
  
  # 创建详细的细胞类型统计表
  create_detailed_stats_table(cellchat_list, output_dir)
  
  # 保存CellChat对象
  saveRDS(cellchat_list, file = file.path(output_dir, "CellChat_Objects.rds"))
  
  message("✓ 完成所有统计表格生成")
  
  return(communication_stats)
}

# 创建详细统计表格的函数
create_detailed_stats_table <- function(cellchat_list, output_dir) {
  
  # 细胞类型详细统计
  celltype_stats_list <- list()
  
  for (group_name in names(cellchat_list)) {
    cellchat <- cellchat_list[[group_name]]
    
    # 获取细胞类型信息
    cell_types <- levels(cellchat@idents)
    celltype_counts <- as.numeric(table(cellchat@idents))
    
    # 创建细胞类型统计表
    celltype_stats <- data.frame(
      Group = group_name,
      Cell_Type = cell_types,
      Cell_Count = celltype_counts,
      Percentage = round(celltype_counts / sum(celltype_counts) * 100, 2),
      stringsAsFactors = FALSE
    )
    
    celltype_stats_list[[group_name]] <- celltype_stats
  }
  
  # 合并所有组的细胞类型统计
  all_celltype_stats <- do.call(rbind, celltype_stats_list)
  
  # 保存细胞类型统计表
  write.csv(all_celltype_stats, file = file.path(output_dir, "Cell_Type_Statistics.csv"), 
            row.names = FALSE)
  
  # 创建通讯矩阵统计表
  if (length(cellchat_list) >= 1) {
    for (group_name in names(cellchat_list)) {
      cellchat <- cellchat_list[[group_name]]
      
      if (!is.null(cellchat@net$count) && sum(cellchat@net$count) > 0) {
        # 通讯数量矩阵
        count_matrix <- cellchat@net$count
        rownames(count_matrix) <- levels(cellchat@idents)
        colnames(count_matrix) <- levels(cellchat@idents)
        
        write.csv(count_matrix, 
                  file = file.path(output_dir, paste0(group_name, "_Communication_Count_Matrix.csv")))
        
        # 通讯强度矩阵
        weight_matrix <- cellchat@net$weight
        rownames(weight_matrix) <- levels(cellchat@idents)
        colnames(weight_matrix) <- levels(cellchat@idents)
        
        write.csv(weight_matrix, 
                  file = file.path(output_dir, paste0(group_name, "_Communication_Strength_Matrix.csv")))
      }
    }
  }
}

# ==============================================================================
# 执行简化的细胞通讯分析
# ==============================================================================

# 执行核心可视化分析
message("=== 开始执行细胞通讯可视化 ===")
stats_summary <- visualize_cellchat_core(cellchat_results)


