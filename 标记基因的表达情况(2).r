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
workDir <- "D:/2.椎间盘衰老/Q1区纯生信7.8分/34.标记基因的表达情况"

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
showGenes <- c("SF3A3", "WDR5B")

# ==============================================================================
# 03. 数据读取和Seurat对象创建
# ==============================================================================

message("正在读取单细胞数据...")

# 读取数据文件
single_cell_data <- read.csv("single_cell_data.csv", row.names = 1)

# 数据质量检查
if (ncol(single_cell_data) == 0) {
  stop("数据文件为空，请检查文件内容！")
}

message(sprintf("数据包含 %d 个基因和 %d 个细胞", 
                nrow(single_cell_data), ncol(single_cell_data)))

# 创建Seurat对象
Peripheral_Blood_Mononuclear_Cells <- CreateSeuratObject(
  counts = single_cell_data, 
  min.cells = analysis_params$min_cells,
  min.features = analysis_params$min_features
)

message("Seurat对象创建完成")

# ==============================================================================
# 04. 数据预处理
# ==============================================================================

message("开始数据预处理...")

# 优化的预处理流程
preprocess_seurat <- function(seurat_obj, params) {
  pb <- progress_bar$new(
    total = 7, 
    format = "  数据处理进度 [:bar] :percent | 耗时: :elapsed | 预计剩余: :eta"
  )
  
  # 数据归一化
  seurat_obj <- NormalizeData(
    object = seurat_obj, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )
  pb$tick()
  
  # 寻找高变基因
  seurat_obj <- FindVariableFeatures(
    object = seurat_obj, 
    selection.method = "vst", 
    nfeatures = params$n_variable_features
  )
  pb$tick()
  
  # 数据标准化
  seurat_obj <- ScaleData(seurat_obj)
  pb$tick()
  
  # PCA分析
  seurat_obj <- RunPCA(
    object = seurat_obj, 
    npcs = 20, 
    features = VariableFeatures(object = seurat_obj)
  )
  pb$tick()
  
  # 邻居查找
  seurat_obj <- FindNeighbors(
    object = seurat_obj, 
    dims = 1:params$pcSelect
  )
  pb$tick()
  
  # 聚类分析
  seurat_obj <- FindClusters(
    object = seurat_obj, 
    resolution = params$cluster_resolution
  )
  pb$tick()
  
  # UMAP降维
  seurat_obj <- RunUMAP(
    object = seurat_obj, 
    dims = 1:params$pcSelect
  )
  pb$tick()
  
  message("数据预处理完成")
  return(seurat_obj)
}

Peripheral_Blood_Mononuclear_Cells <- preprocess_seurat(
  Peripheral_Blood_Mononuclear_Cells, 
  analysis_params
)

# ==============================================================================
# 05. 基础可视化
# ==============================================================================

message("生成基础可视化图形...")

# 优化的可视化函数
create_umap_plot <- function(seurat_obj, filename, title = "UMAP Clustering Results") {
  pdf(file = filename, width = 8, height = 6)
  print(
    DimPlot(
      object = seurat_obj, 
      reduction = "umap", 
      pt.size = 1.5, 
      label = TRUE,
      label.size = 4
    ) + 
      ggtitle(title) +
      theme_minimal()
  )
  dev.off()
  message(sprintf("保存UMAP图: %s", filename))
}

create_umap_plot(Peripheral_Blood_Mononuclear_Cells, "UMAP_Clusters.pdf")

# 保存聚类信息
write.table(
  Peripheral_Blood_Mononuclear_Cells$seurat_clusters, 
  file = "umapCluster.txt", 
  quote = FALSE, 
  sep = "\t", 
  col.names = FALSE
)

# ==============================================================================
# 06. Marker基因分析
# ==============================================================================

message("进行Marker基因分析...")

# 寻找所有聚类的marker基因
find_all_markers <- function(seurat_obj, params) {
  markers <- FindAllMarkers(
    object = seurat_obj, 
    only.pos = FALSE, 
    min.pct = 0.25, 
    logfc.threshold = params$logFCfilter
  )
  
  # 筛选显著的marker基因
  sig_markers <- markers[
    (abs(markers$avg_log2FC) > params$logFCfilter & 
       markers$p_val_adj < params$adjPvalFilter), 
  ]
  
  # 保存结果
  write.table(
    sig_markers, 
    file = "clusterMarkers.txt", 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
  )
  
  message(sprintf("找到 %d 个显著marker基因", nrow(sig_markers)))
  return(list(all_markers = markers, sig_markers = sig_markers))
}

marker_results <- find_all_markers(Peripheral_Blood_Mononuclear_Cells, analysis_params)

# ==============================================================================
# 07. 热图可视化
# ==============================================================================

message("生成热图...")

# 优化的热图生成
create_heatmap <- function(seurat_obj, markers_df, filename = "umapHeatmap.pdf") {
  top10 <- markers_df %>% 
    group_by(cluster) %>% 
    top_n(n = 10, wt = avg_log2FC)
  
  pdf(file = filename, width = 12, height = 9)
  print(
    DoHeatmap(
      object = seurat_obj, 
      features = top10$gene
    ) + 
      NoLegend()
  )
  dev.off()
  message(sprintf("保存热图: %s", filename))
}

create_heatmap(Peripheral_Blood_Mononuclear_Cells, marker_results$all_markers)

# ==============================================================================
# 08. 目标基因表达可视化
# ==============================================================================

message("生成目标基因表达图...")

# 优化的基因可视化函数
visualize_genes <- function(seurat_obj, genes) {
  # 检查基因是否存在
  available_genes <- genes[genes %in% rownames(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))]
  missing_genes <- setdiff(genes, available_genes)
  
  if (length(missing_genes) > 0) {
    message(paste("以下基因未在数据集中找到:", paste(missing_genes, collapse = ", ")))
  }
  
  if (length(available_genes) == 0) {
    warning("所有目标基因均未在数据集中找到")
    return()
  }
  
  # 小提琴图 - 所有基因
  pdf(file = "markerViolin_All.pdf", width = 15, height = 10)
  print(
    VlnPlot(
      object = seurat_obj, 
      features = available_genes,
      ncol = 3
    ) + 
      ggtitle("Gene Expression Distribution Across Clusters")
  )
  dev.off()
  
  # 特征图 - 所有基因
  pdf(file = "markerFeature_All.pdf", width = 15, height = 12)
  print(
    FeaturePlot(
      object = seurat_obj, 
      features = available_genes, 
      cols = c("#E5E5E5", "#FF0000"), 
      reduction = "umap",
      ncol = 3
    )
  )
  dev.off()
  
  # 点图
  pdf(file = "markerDotPlot.pdf", width = 10, height = 8)
  print(
    DotPlot(
      object = seurat_obj, 
      features = available_genes
    ) + 
      coord_flip() + 
      scale_color_gradientn(colors = c("#F0F0F0", "#F9B5B5", "#F46A6A", "#D7301F")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) + 
      labs(title = "Expression of Marker Genes Across Clusters")
  )
  dev.off()
}

visualize_genes(Peripheral_Blood_Mononuclear_Cells, showGenes)

# ==============================================================================
# 09. 细胞类型注释
# ==============================================================================

message("进行细胞类型注释...")

# 使用预定义cluster注释的细胞类型注释函数
perform_cell_annotation <- function(seurat_obj) {
  # 读取预定义的cluster注释
  cluster_ann <- read.table("clusterAnn.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # 获取当前聚类信息
  clusters <- seurat_obj@meta.data$seurat_clusters

  # 创建cluster到cell type的映射
  cluster_to_celltype <- setNames(cluster_ann$labels, cluster_ann$Cluster)

  # 为每个细胞分配细胞类型
  cell_ann <- cluster_to_celltype[as.character(clusters)]
  names(cell_ann) <- colnames(seurat_obj)

  # 创建细胞水平注释文件
  cell_ann_df <- data.frame(
    Cell = colnames(seurat_obj),
    labels = cell_ann
  )

  write.table(
    cell_ann_df,
    file = "cellAnn.txt",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  # 重新标记聚类使用预定义的细胞类型
  new_labels <- cluster_to_celltype[as.character(levels(seurat_obj))]
  names(new_labels) <- levels(seurat_obj)
  seurat_obj <- RenameIdents(seurat_obj, new_labels)

  # 添加细胞类型到元数据
  seurat_obj <- AddMetaData(
    seurat_obj,
    metadata = cell_ann,
    col.name = "cell_type"
  )

  return(list(seurat_obj = seurat_obj, cell_ann = cell_ann))
}

annotation_results <- perform_cell_annotation(Peripheral_Blood_Mononuclear_Cells)
Peripheral_Blood_Mononuclear_Cells <- annotation_results$seurat_obj
cellAnn <- annotation_results$cell_ann

#

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(viridis)
create_individual_gene_scatter_plots_simple <- function(seurat_obj, genes, base_output_dir = "Gene_Expression_Plots") {
  
  if (!dir.exists(base_output_dir)) {
    dir.create(base_output_dir, recursive = TRUE)
  }
  
  available_genes <- genes[genes %in% rownames(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))]
  if (length(available_genes) == 0) {
    stop("所有目标基因均未在数据集中找到")
  }
  

  # 获取元数据
  plot_data <- data.frame(
    Cell = colnames(seurat_obj),
    UMAP_1 = seurat_obj@reductions$umap@cell.embeddings[, 1],
    UMAP_2 = seurat_obj@reductions$umap@cell.embeddings[, 2],
    Cell_Type = seurat_obj@meta.data$cell_type,
    Group = ifelse(rownames(seurat_obj@meta.data) %in%
                     colnames(seurat_obj)[sample(ncol(seurat_obj), ncol(seurat_obj)/2)],
                   "CON", "IDD")
  )

  expression_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  unique_cell_types <- unique(plot_data$Cell_Type[!is.na(plot_data$Cell_Type)])
  coord_names <- c("UMAP_1", "UMAP_2")

  for (gene in available_genes) {
    gene_dir <- file.path(base_output_dir, paste0(gene, "_Expression"))
    if (!dir.exists(gene_dir)) {
      dir.create(gene_dir, recursive = TRUE)
    }

    plot_data$Expression <- as.numeric(expression_data[gene, ])

    for (cell_type in unique_cell_types) {
      safe_cell_type <- gsub("[^A-Za-z0-9_-]", "_", cell_type)
      cell_dir <- file.path(gene_dir, paste0("Cluster_", safe_cell_type))
      if (!dir.exists(cell_dir)) {
        dir.create(cell_dir, recursive = TRUE)
      }
      
      current_data <- plot_data[plot_data$Cell_Type == cell_type & !is.na(plot_data$Cell_Type), ]
      
      if (nrow(current_data) == 0) {
        next
      }
      
      safe_gene <- gsub("[^A-Za-z0-9_-]", "_", gene)
      
      # 只生成分布图
      plot1_file <- file.path(cell_dir, paste0("01_", safe_gene, "_Cluster_", safe_cell_type, "_distribution.pdf"))
      pdf(file = plot1_file, width = 12, height = 9)
      
      colors_sci <- c("CON" = "#2E86AB", "IDD" = "#E63946")
      
      p1 <- ggplot() +
        geom_point(data = plot_data, 
                   aes_string(x = coord_names[1], y = coord_names[2]), 
                   color = "#E8E8E8", size = 0.2, alpha = 0.4, stroke = 0) +
        geom_point(data = current_data, 
                   aes_string(x = coord_names[1], y = coord_names[2], color = "Group"), 
                   size = 1.8, alpha = 0.85, stroke = 0) +
        scale_color_manual(values = colors_sci, name = "Group") +
        labs(title = paste("Cell Distribution: Cluster_", cell_type),
             subtitle = paste(gene, "| n =", nrow(current_data), "cells")) +
        theme_classic()
      
      print(p1)
      dev.off()
      
      message(paste("✓ 生成分布图: Cluster_", safe_cell_type, "-", gene))
    }
  }
}
create_individual_gene_scatter_plots_simple(Peripheral_Blood_Mononuclear_Cells, showGenes)

