#                        单细胞RNA测序数据分析流程                              #
#                    Single Cell RNA-seq Analysis Pipeline                    #
################################################################################

# ==============================================================================
# 01. 包加载
# ==============================================================================
# ==============================================================================
# R包加载汇总
# ==============================================================================

# 核心单细胞分析包
library(Seurat)        # 单细胞数据分析核心包
library(limma)         # 差异表达分析
library(SingleR)       # 单细胞类型注释
library(celldex)       # 细胞类型注释数据集
library(monocle)       # 单细胞轨迹分析

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

# 统计分析包
library(rstatix)       # 统计检验

# 交互和辅助包
library(DT)            # 交互式表格
library(progress)      # 进度条显示

# 可选的文件输出包
library(openxlsx)      # Excel文件读写（如果需要Excel输出）

message("所有必需的R包已成功加载完成！")

# ==============================================================================
# 02. 参数设置和工作目录
# ==============================================================================

# 设置工作目录 - 请根据实际路径修改
workDir <- "D:/2.椎间盘衰老(1)/Q1区纯生信7.8分/32.单细胞分析-聚类注释、基因表达和分布"

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

###先去31文件夹把seurat.rda复制到32文件夹
load(file = "seurat.rda")

# 感兴趣的基因列表
showGenes <- c("SF3A3", "GSTZ1")

# ==============================================================================
# 03. 数据读取和Seurat对象创建
# ==============================================================================
# 细胞聚类
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:15) # dims参数根据ElbowPlot结果选择
filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.5) # resolution值影响聚类数量，可以尝试不同的分辨率

# 运行UMAP进行非线性降维可视化
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:15)

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

create_umap_plot(filtered_seurat, "UMAP_Clusters.pdf")

# 保存聚类信息
write.table(
  filtered_seurat$seurat_clusters, 
  file = "umapCluster.txt", 
  quote = FALSE, 
  sep = "\t", 
  col.names = FALSE
)

# ==============================================================================
# 06. Marker基因分析####这一步很久
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

marker_results <- find_all_markers(filtered_seurat, analysis_params)

# ==============================================================================
# 07. 热图可视化
# ==============================================================================

message("生成热图...")

# 优化的热图生成
create_heatmap <- function(seurat_obj, markers_df, filename = "umapHeatmap.pdf") {
  top10 <- markers_df %>% 
    group_by(cluster) %>% 
    top_n(n = 10, wt = avg_log2FC)
  
  pdf(file = filename, width = 20, height = 20)
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

create_heatmap(filtered_seurat, marker_results$all_markers)

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
  pdf(file = "markerViolin_All.pdf", width = 12, height = 10)
  print(
    VlnPlot(
      object = seurat_obj, 
      features = available_genes,
      ncol = 2
    ) + 
      ggtitle("Gene Expression Distribution Across Clusters")
  )
  dev.off()
  
  # 特征图 - 所有基因
  pdf(file = "markerFeature_All.pdf", width = 15, height = 10)
  print(
    FeaturePlot(
      object = seurat_obj, 
      features = available_genes, 
      cols = c("#E5E5E5", "#FF0000"), 
      reduction = "umap",
      ncol = 2
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

visualize_genes(filtered_seurat, showGenes)

# ==============================================================================
# 09. 细胞类型注释
# ==============================================================================

message("进行细胞类型注释...")

# 优化的细胞类型注释函数
perform_cell_annotation <- function(seurat_obj) {
  # 获取表达数据和聚类信息
  counts <- GetAssayData(object = seurat_obj, slot = "data")
  clusters <- seurat_obj@meta.data$seurat_clusters
  
  # 使用SingleR进行细胞类型注释
  ref <- celldex::HumanPrimaryCellAtlasData()
  
  # 聚类水平注释
  singler_clusters <- SingleR(
    test = counts, 
    ref = ref, 
    labels = ref$label.main, 
    clusters = clusters
  )
  
  cluster_ann <- data.frame(
    Cluster = rownames(singler_clusters), 
    labels = singler_clusters$labels
  )
  
  write.table(
    cluster_ann, 
    file = "clusterAnn.txt", 
    quote = FALSE, 
    sep = "\t", 
    row.names = FALSE
  )
  
  # 细胞水平注释
  singler_cells <- SingleR(
    test = counts, 
    ref = ref, 
    labels = ref$label.main
  )
  
  cell_ann <- singler_cells$labels
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
  
  # 重新标记聚类
  new_labels <- singler_clusters$labels
  new_labels <- gsub("_", " ", new_labels)
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

annotation_results <- perform_cell_annotation(filtered_seurat)
filtered_seurat <- annotation_results$seurat_obj
cellAnn <- annotation_results$cell_ann

# ==============================================================================
# 10. 注释后可视化
# ==============================================================================

message("生成注释后的可视化图...")

# 细胞类型注释UMAP图
pdf(file = "cellTypeAnnotation.pdf", width = 10, height = 8)
print(
  DimPlot(
    filtered_seurat, 
    reduction = "umap", 
    pt.size = 1.5, 
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle("Cell Type Annotations") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)
dev.off()

# 添加分组信息
Type <- gsub("(.*?)\\..*", "\\1", colnames(filtered_seurat))
names(Type) <- colnames(filtered_seurat)
filtered_seurat <- AddMetaData(
  object = filtered_seurat, 
  metadata = Type, 
  col.name = "Type"
)

# 分组可视化
pdf(file = "groupComparison.pdf", width = 14, height = 6)
print(
  DimPlot(
    filtered_seurat, 
    reduction = "umap", 
    pt.size = 1, 
    label = TRUE, 
    split.by = "Type"
  ) +
    ggtitle("Cell Types: CON vs IDD") +
    theme_minimal()
)
dev.off()

# ==============================================================================
# 11. 差异表达分析###这个也很久
# ==============================================================================

message("进行组间差异表达分析...")

# 优化的差异表达分析函数
perform_differential_analysis <- function(seurat_obj, cell_ann, type_info, params) {
  # 添加组合标签
  groups <- paste0(type_info, "_", cell_ann)
  names(groups) <- colnames(seurat_obj)
  seurat_obj <- AddMetaData(
    object = seurat_obj, 
    metadata = groups, 
    col.name = "group"
  )
  
  # 为每种细胞类型进行差异分析
  unique_cell_types <- unique(cell_ann)
  
  for (cellName in unique_cell_types) {
    con_name <- paste0("CON_", cellName)
    treat_name <- paste0("IDD_", cellName)
    
    # 检查标签是否存在
    if (!(con_name %in% seurat_obj$group) | !(treat_name %in% seurat_obj$group)) {
      message(sprintf("跳过 %s - 缺少CON组或IDD组", cellName))
      next
    }
    
    # 获取细胞
    con_cells <- WhichCells(seurat_obj, expression = group == con_name)
    treat_cells <- WhichCells(seurat_obj, expression = group == treat_name)
    
    if (length(con_cells) >= 3 & length(treat_cells) >= 3) {
      message(sprintf("分析 %s - CON组: %d 细胞, IDD组: %d 细胞", 
                      cellName, length(con_cells), length(treat_cells)))
      
      tryCatch({
        pbmc_markers <- FindMarkers(
          seurat_obj, 
          ident.1 = treat_cells, 
          ident.2 = con_cells, 
          group.by = 'group', 
          logfc.threshold = 0.1
        )
        
        sig_markers_group <- pbmc_markers[
          (abs(pbmc_markers$avg_log2FC) > params$logFCfilter & 
             pbmc_markers$p_val_adj < params$adjPvalFilter), 
        ]
        
        if (nrow(sig_markers_group) > 0) {
          sig_markers_group <- cbind(Gene = rownames(sig_markers_group), sig_markers_group)
          
          # 清理文件名中的特殊字符
          safe_cell_name <- gsub("[^A-Za-z0-9_]", "_", cellName)
          
          write.table(
            sig_markers_group, 
            file = paste0(safe_cell_name, "_diffGenes.txt"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE
          )
          
          message(sprintf("保存差异基因结果: %s_diffGenes.txt", safe_cell_name))
        } else {
          message(sprintf("未找到 %s 的显著差异基因", cellName))
        }
      }, error = function(e) {
        message(sprintf("分析 %s 时出错: %s", cellName, e$message))
      })
    } else {
      message(sprintf("跳过 %s - 细胞数量不足", cellName))
    }
  }
  
  return(seurat_obj)
}

filtered_seurat <- perform_differential_analysis(
  filtered_seurat, 
  cellAnn, 
  Type, 
  analysis_params
)

# ==============================================================================
# 12. 综合统计和报告
# ==============================================================================

message("生成分析统计报告...")

# 优化的统计报告函数
generate_statistics_report <- function(seurat_obj) {
  # 细胞类型统计
  cell_type_counts <- table(Idents(seurat_obj), 
                            seurat_obj$Type)
  
  cell_type_stats <- data.frame(
    CellType = rownames(cell_type_counts),
    Control_Count = cell_type_counts[, "CON"],
    Treatment_Count = cell_type_counts[, "IDD"],
    Total_Count = rowSums(cell_type_counts)
  )
  
  cell_type_stats$Control_Percentage <- round(
    cell_type_stats$Control_Count / sum(cell_type_stats$Control_Count) * 100, 2
  )
  cell_type_stats$Treatment_Percentage <- round(
    cell_type_stats$Treatment_Count / sum(cell_type_stats$Treatment_Count) * 100, 2
  )
  
  write.csv(cell_type_stats, "CellType_Statistics.csv", row.names = FALSE)
  
  # 打印分析摘要
  message("\n=== 分析完成摘要 ===")
  message(sprintf("总细胞数: %d", ncol(seurat_obj)))
  message(sprintf("总基因数: %d", nrow(seurat_obj)))
  message(sprintf("聚类数: %d", length(levels(seurat_obj$seurat_clusters))))
  message(sprintf("细胞类型数: %d", length(unique(seurat_obj$cell_type))))
  
  return(cell_type_stats)
}

stats_report <- generate_statistics_report(filtered_seurat)

message("\n所有分析已完成！结果文件已保存到工作目录。")
message("主要输出文件包括：")
message("- UMAP聚类图和细胞类型注释图")
message("- 基因表达分布图和点图")
message("- 差异表达分析结果")
message("- 统计报告和Seurat对象")
###########################################################
##########################################################
# 保存Seurat对象
save(filtered_seurat, cellAnn,  analysis_params, annotation_results, marker_results, 
     file = "filtered_seurat_1.rda")
load(file = "filtered_seurat_1.rda")

####################################################################
# ==============================================================================
# 13. 注释细胞与Top4 Marker基因相关性分析
# ==============================================================================

message("开始注释细胞与Top4 Marker基因相关性分析...")

# 创建输出目录
output_dirs <- c("CellType_Top4_Markers_Analysis", 
                 "CellType_Top4_DotPlots", 
                 "CellType_Top4_Statistics")

for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(sprintf("创建目录: %s", dir))
  }
}

# 优化的Top4 Marker基因分析函数
analyze_celltype_top4_markers <- function(seurat_obj, logfc_threshold = 0.5, 
                                          adj_pval_threshold = 0.05) {
  
  message("正在为每个细胞类型寻找Top4 marker基因...")
  
  # 使用当前的细胞类型注释作为identity
  current_ident <- Idents(seurat_obj)
  
  # 寻找每个细胞类型的marker基因
  celltype_markers <- FindAllMarkers(
    object = seurat_obj, 
    only.pos = TRUE,  # 只保留正向表达的基因
    min.pct = 0.25, 
    logfc.threshold = logfc_threshold,
    test.use = "wilcox",
    verbose = TRUE
  )
  
  # 筛选显著的marker基因
  sig_celltype_markers <- celltype_markers[
    (celltype_markers$p_val_adj < adj_pval_threshold & 
       celltype_markers$avg_log2FC > logfc_threshold), 
  ]
  
  message(sprintf("找到 %d 个显著的细胞类型marker基因", nrow(sig_celltype_markers)))
  
  # 为每个细胞类型选择Top4 marker基因
  top4_markers_summary <- sig_celltype_markers %>% 
    group_by(cluster) %>% 
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 4) %>%  # 选择前4个，如果不足4个则全部选择
    ungroup()
  
  # 保存详细的marker基因结果
  write.csv(
    sig_celltype_markers, 
    "CellType_Top4_Markers_Analysis/All_CellType_Markers.csv", 
    row.names = FALSE
  )
  
  write.csv(
    top4_markers_summary, 
    "CellType_Top4_Markers_Analysis/Top4_Markers_Summary.csv", 
    row.names = FALSE
  )
  
  return(list(
    all_markers = sig_celltype_markers,
    top4_markers = top4_markers_summary
  ))
}

# 执行Top4 marker基因分析
top4_results <- analyze_celltype_top4_markers(
  filtered_seurat,
  logfc_threshold = 0.5,
  adj_pval_threshold = 0.05
)

# 生成每个细胞类型的详细统计报告
generate_celltype_statistics <- function(seurat_obj, top4_markers_df) {
  
  message("生成每个细胞类型的详细统计...")
  
  # 获取所有细胞类型
  unique_celltypes <- levels(Idents(seurat_obj))
  
  # 创建详细统计数据框
  detailed_stats <- data.frame()
  
  for (celltype in unique_celltypes) {
    # 获取该细胞类型的top4基因
    celltype_genes <- top4_markers_df %>% 
      filter(cluster == celltype) %>%
      pull(gene)
    
    # 如果没有找到基因，记录并跳过
    if (length(celltype_genes) == 0) {
      warning(sprintf("细胞类型 %s 没有找到显著的marker基因", celltype))
      
      # 添加空记录
      empty_record <- data.frame(
        CellType = celltype,
        Gene = "No significant markers found",
        avg_log2FC = NA,
        pct.1 = NA,
        pct.2 = NA,
        p_val_adj = NA,
        Rank = NA,
        Gene_Count = 0
      )
      detailed_stats <- rbind(detailed_stats, empty_record)
      next
    }
    
    # 获取该细胞类型的统计信息
    celltype_stats <- top4_markers_df %>% 
      filter(cluster == celltype) %>%
      arrange(desc(avg_log2FC)) %>%
      mutate(Rank = row_number()) %>%
      select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj, Rank)
    
    # 添加基因数量信息
    celltype_stats$Gene_Count <- nrow(celltype_stats)
    
    # 重命名列
    names(celltype_stats)[names(celltype_stats) == "cluster"] <- "CellType"
    names(celltype_stats)[names(celltype_stats) == "gene"] <- "Gene"
    
    detailed_stats <- rbind(detailed_stats, celltype_stats)
  }
  
  # 保存详细统计
  write.csv(
    detailed_stats, 
    "CellType_Top4_Statistics/CellType_Top4_Detailed_Statistics.csv", 
    row.names = FALSE
  )
  
  return(detailed_stats)
}

# 生成统计报告
detailed_statistics <- generate_celltype_statistics(
  filtered_seurat, 
  top4_results$top4_markers
)

# 生成综合Top4气泡图
create_comprehensive_dotplot <- function(seurat_obj, top4_markers_df, 
                                         filename = "CellType_Top4_Comprehensive_DotPlot.pdf") {
  
  # 获取所有top4基因
  all_top4_genes <- unique(top4_markers_df$gene)
  
  if (length(all_top4_genes) == 0) {
    warning("没有找到任何Top4基因，跳过综合气泡图生成")
    return()
  }
  
  message(sprintf("生成包含 %d 个Top4基因的综合气泡图", length(all_top4_genes)))
  
  # 创建综合气泡图
  pdf(file = filename, width =12, height = 15)
  
  tryCatch({
    plot <- DotPlot(
      object = seurat_obj, 
      features = all_top4_genes
    ) + 
      coord_flip() + 
      scale_color_gradientn(
        colors = c("#F7F7F7", "#D9F0A3", "#ADDD8E", "#78C679", 
                   "#41AB5D", "#238443", "#005A32"),
        name = "Average\nExpression"
      ) +
      scale_size_continuous(
        name = "Percent\nExpressed",
        range = c(1, 12),
        breaks = c(0, 25, 50, 75, 100),
        labels = c("0%", "25%", "50%", "75%", "100%")
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20)),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        plot.margin = margin(20, 20, 20, 20)
      ) + 
      labs(
        title = "Top 4 Marker Genes Expression Across Cell Types",
        x =  "Cell Types", 
        y = "Marker Genes"
      )
    
    print(plot)
    
  }, error = function(e) {
    message(sprintf("生成综合气泡图时出错: %s", e$message))
  })
  
  dev.off()
  message(sprintf("综合气泡图已保存: %s", filename))
}

# 生成综合气泡图
create_comprehensive_dotplot(
  filtered_seurat, 
  top4_results$top4_markers,
  "CellType_Top4_DotPlots/CellType_Top4_Comprehensive_DotPlot.pdf"
)

# 为每个细胞类型生成单独的气泡图
create_individual_celltype_dotplots <- function(seurat_obj, top4_markers_df) {
  
  message("为每个细胞类型生成单独的Top4基因气泡图...")
  
  # 获取所有细胞类型
  unique_celltypes <- levels(Idents(seurat_obj))
  
  for (celltype in unique_celltypes) {
    # 获取当前细胞类型的top4基因
    current_top4 <- top4_markers_df %>% 
      filter(cluster == celltype) %>% 
      arrange(desc(avg_log2FC)) %>%
      pull(gene)
    
    # 如果没有找到基因，跳过
    if (length(current_top4) == 0) {
      message(sprintf("细胞类型 %s 没有找到显著的marker基因，跳过气泡图生成", celltype))
      next
    }
    
    # 清理细胞类型名称用于文件名
    clean_celltype <- gsub("[^A-Za-z0-9_]", "_", celltype)
    pdf_file <- file.path("CellType_Top4_DotPlots", 
                          paste0(clean_celltype, "_Top4_Markers_DotPlot.pdf"))
    
    message(sprintf("生成 %s 的气泡图，包含 %d 个基因", celltype, length(current_top4)))
    
    tryCatch({
      pdf(file = pdf_file, width = max(8, length(current_top4) * 1.5), height = 6)
      
      plot <- DotPlot(
        object = seurat_obj, 
        features = current_top4
      ) + 
        coord_flip() + 
        scale_color_gradientn(
          colors = c("#F0F8FF", "#B0E0E6", "#87CEEB", "#4682B4", "#1E90FF", "#0000CD"),
          name = "Average\nExpression"
        ) +
        scale_size_continuous(
          name = "Percent\nExpressed",
          range = c(2, 15)
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 13, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                                    margin = margin(b = 20)),
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 10),
          legend.position = "right",
          panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.3),
          panel.grid.minor = element_blank()
        ) + 
        labs(
          title = sprintf("Top %d Marker Genes for %s", length(current_top4), celltype),
          x = "Marker Genes", 
          y = "Cell Types"
        )
      
      print(plot)
      dev.off()
      
      message(sprintf("已生成细胞类型 %s 的气泡图: %s", celltype, pdf_file))
      
      # 同时保存该细胞类型的基因列表
      celltype_genes_df <- top4_markers_df %>% 
        filter(cluster == celltype) %>%
        arrange(desc(avg_log2FC)) %>%
        select(gene, avg_log2FC, pct.1, pct.2, p_val_adj)
      
      csv_file <- file.path("CellType_Top4_Statistics", 
                            paste0(clean_celltype, "_Top4_Genes.csv"))
      
      write.csv(celltype_genes_df, csv_file, row.names = FALSE)
      message(sprintf("已保存细胞类型 %s 的基因列表: %s", celltype, csv_file))
      
    }, error = function(e) {
      message(sprintf("生成细胞类型 %s 的气泡图时出错: %s", celltype, e$message))
    })
  }
}

# 生成每个细胞类型的单独气泡图
create_individual_celltype_dotplots(
  filtered_seurat, 
  top4_results$top4_markers
)

# 生成Top4基因的表达热图
create_top4_heatmap <- function(seurat_obj, top4_markers_df, 
                                filename = "CellType_Top4_Expression_Heatmap.pdf") {
  
  all_top4_genes <- unique(top4_markers_df$gene)
  
  if (length(all_top4_genes) == 0) {
    warning("没有找到任何Top4基因，跳过热图生成")
    return()
  }
  
  message(sprintf("生成包含 %d 个Top4基因的表达热图", length(all_top4_genes)))
  
  tryCatch({
    # 计算平均表达
    avg_exp <- AverageExpression(
      seurat_obj, 
      features = all_top4_genes,
      return.seurat = FALSE
    )
    
    exp_matrix <- as.matrix(avg_exp$RNA)
    
    # 标准化表达矩阵
    exp_scaled <- t(scale(t(exp_matrix)))
    
    pdf(file = filename, width = 12, height = max(8, length(all_top4_genes) * 0.3))
    
    # 使用基础heatmap函数
    heatmap(
      exp_scaled,
      col = colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", 
                               "#E0F3F8", "#FFFFCC", "#FEE090", "#FDAE61", 
                               "#F46D43", "#D73027", "#A50026"))(100),
      main = "Top 4 Marker Genes Expression Heatmap",
      cexRow = 0.8,
      cexCol = 0.9,
      margins = c(10, 12)
    )
    
    dev.off()
    message(sprintf("表达热图已保存: %s", filename))
    
  }, error = function(e) {
    message(sprintf("生成热图时出错: %s", e$message))
  })
}

# 生成Top4基因表达热图
create_top4_heatmap(
  filtered_seurat, 
  top4_results$top4_markers,
  "CellType_Top4_DotPlots/CellType_Top4_Expression_Heatmap.pdf"
)

# 生成最终的分析摘要报告
generate_top4_analysis_summary <- function(top4_markers_df, detailed_stats_df) {
  
  message("生成Top4 Marker基因分析摘要报告...")
  
  summary_stats <- detailed_stats_df %>%
    group_by(CellType) %>%
    summarise(
      Marker_Gene_Count = sum(!is.na(avg_log2FC)),
      Max_LogFC = max(avg_log2FC, na.rm = TRUE),
      Min_Pval_adj = min(p_val_adj, na.rm = TRUE),
      Avg_Pct1 = mean(pct.1, na.rm = TRUE),
      Avg_Pct2 = mean(pct.2, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(Marker_Gene_Count))
  
  # 保存摘要统计
  write.csv(
    summary_stats, 
    "CellType_Top4_Statistics/CellType_Analysis_Summary.csv", 
    row.names = FALSE
  )
  
  # 打印摘要信息到控制台
  message("\n=== Top4 Marker基因分析完成摘要 ===")
  message(sprintf("分析的细胞类型数: %d", nrow(summary_stats)))
  message(sprintf("总计Top4基因数: %d", sum(summary_stats$Marker_Gene_Count)))
  
  print(summary_stats)
  
  return(summary_stats)
}

# 生成分析摘要
analysis_summary <- generate_top4_analysis_summary(
  top4_results$top4_markers, 
  detailed_statistics
)

message("\n=== Marker基因相关性分析完成！ ===")

# ==============================================================================
# 14. 对照组与实验组细胞比例分析
# ==============================================================================

message("开始CON组与IDD组细胞比例分析...")

# 创建输出目录
proportion_dirs <- c("CellType_Proportion_Analysis", 
                     "Group_Proportion_Plots", 
                     "Proportion_Statistics")

for (dir in proportion_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(sprintf("创建目录: %s", dir))
  }
}

# 细胞比例分析函数
analyze_celltype_proportions <- function(seurat_obj) {
  
  message("计算各组细胞类型比例...")
  
  # 获取细胞类型和分组信息
  cell_types <- Idents(seurat_obj)
  groups <- seurat_obj@meta.data$Type
  
  # 创建交叉表
  cross_table <- table(cell_types, groups)
  
  # 计算各种比例
  # 1. 每组内各细胞类型的比例
  prop_within_group <- prop.table(cross_table, margin = 2) * 100
  
  # 2. 每个细胞类型在各组中的分布比例
  prop_within_celltype <- prop.table(cross_table, margin = 1) * 100
  
  # 3. 总体比例
  prop_total <- prop.table(cross_table) * 100
  
  # 转换为数据框格式
  # 组内比例数据框
  prop_within_group_df <- as.data.frame(prop_within_group)
  colnames(prop_within_group_df) <- c("CellType", "Group", "Percentage_Within_Group")
  
  # 细胞类型内比例数据框
  prop_within_celltype_df <- as.data.frame(prop_within_celltype)
  colnames(prop_within_celltype_df) <- c("CellType", "Group", "Percentage_Within_CellType")
  
  # 总体比例数据框
  prop_total_df <- as.data.frame(prop_total)
  colnames(prop_total_df) <- c("CellType", "Group", "Percentage_Total")
  
  # 细胞计数数据框
  counts_df <- as.data.frame(cross_table)
  colnames(counts_df) <- c("CellType", "Group", "Cell_Count")
  
  # 合并所有信息
  comprehensive_stats <- merge(counts_df, prop_within_group_df, by = c("CellType", "Group"))
  comprehensive_stats <- merge(comprehensive_stats, prop_within_celltype_df, by = c("CellType", "Group"))
  comprehensive_stats <- merge(comprehensive_stats, prop_total_df, by = c("CellType", "Group"))
  
  # 添加总计信息
  total_control <- sum(groups == "CON")
  total_treat <- sum(groups == "IDD")
  total_cells <- length(groups)
  
  comprehensive_stats$Total_Cells_In_Group <- ifelse(
    comprehensive_stats$Group == "CON", 
    total_control, 
    total_treat
  )
  
  comprehensive_stats$Total_Cells_Overall <- total_cells
  
  # 重新排列列顺序
  comprehensive_stats <- comprehensive_stats[, c(
    "CellType", "Group", "Cell_Count", "Total_Cells_In_Group",
    "Percentage_Within_Group", "Percentage_Within_CellType", 
    "Percentage_Total", "Total_Cells_Overall"
  )]
  
  # 按细胞类型和组排序
  comprehensive_stats <- comprehensive_stats[order(comprehensive_stats$CellType, comprehensive_stats$Group), ]
  
  return(list(
    cross_table = cross_table,
    prop_within_group = prop_within_group,
    prop_within_celltype = prop_within_celltype,
    comprehensive_stats = comprehensive_stats
  ))
}

# 执行比例分析
proportion_results <- analyze_celltype_proportions(filtered_seurat)

# 保存详细统计结果
write.csv(
  proportion_results$comprehensive_stats,
  "Proportion_Statistics/Comprehensive_CellType_Proportions.csv",
  row.names = FALSE
)

# 生成组内比例柱状图（Control vs Treat中各细胞类型的占比）
create_within_group_proportion_plot <- function(stats_df, filename) {
  
  message("生成组内细胞类型比例图...")
  
  pdf(file = filename, width = 7, height = 8)
  
  p <- ggplot(stats_df, aes(x = Group, y = Percentage_Within_Group, fill = CellType)) +
    geom_col(position = "stack", alpha = 0.8, color = "white", linewidth = 0.3) +
    scale_fill_viridis_d(name = "Cell Type", option = "turbo") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 20)),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = "Cell Type Proportions Within Each Group",
      x = "Group",
      y = "Percentage (%)",
      subtitle = " "
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, 10),
      labels = paste0(seq(0, 100, 10), "%")
    ) +
    geom_text(
      aes(label = ifelse(Percentage_Within_Group > 3, 
                         paste0(round(Percentage_Within_Group, 1), "%"), "")),
      position = position_stack(vjust = 0.5),
      size = 3.5,
      color = "white",
      fontface = "bold"
    )
  
  print(p)
  dev.off()
  
  message(sprintf("组内比例图已保存: %s", filename))
}

# 生成细胞类型间比例对比图（每个细胞类型在Control vs Treat中的分布）
create_between_group_proportion_plot <- function(stats_df, filename) {
  
  message("生成细胞类型间组别比例图...")
  
  pdf(file = filename, width = 16, height = 10)
  
  p <- ggplot(stats_df, aes(x = CellType, y = Percentage_Within_CellType, fill = Group)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    scale_fill_manual(
      values = c("Control" = "#2E86AB", "Treat" = "#A23B72"),
      name = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 20)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      title = "Group Distribution Within Each Cell Type",
      x = "Cell Type",
      y = "Percentage (%)",
      subtitle = "Bar chart showing Control vs Treatment distribution for each cell type"
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, 10),
      labels = paste0(seq(0, 100, 10), "%")
    ) +
    geom_text(
      aes(label = paste0(round(Percentage_Within_CellType, 1), "%")),
      position = position_dodge(width = 0.7),
      vjust = -0.3,
      size = 3.5,
      fontface = "bold"
    )
  
  print(p)
  dev.off()
  
  message(sprintf("细胞类型间比例图已保存: %s", filename))
}

# 生成绝对数量对比图
create_absolute_count_plot <- function(stats_df, filename) {
  
  message("生成细胞绝对数量对比图...")
  
  pdf(file = filename, width = 16, height = 10)
  
  p <- ggplot(stats_df, aes(x = CellType, y = Cell_Count, fill = Group)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    scale_fill_manual(
      values = c("CON" = "#1B4F72", "IDD" = "#C0392B"),
      name = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 20)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      title = "Absolute Cell Counts by Cell Type and Group",
      x = "Cell Type",
      y = "Cell Count",
      subtitle = "Bar chart showing actual cell numbers for each cell type in Control vs Treatment"
    ) +
    geom_text(
      aes(label = Cell_Count),
      position = position_dodge(width = 0.7),
      vjust = -0.3,
      size = 3.5,
      fontface = "bold"
    )
  
  print(p)
  dev.off()
  
  message(sprintf("绝对数量图已保存: %s", filename))
}

# 生成热图样式的比例图
create_proportion_heatmap <- function(prop_matrix, filename, title) {
  
  message(sprintf("生成比例热图: %s", title))
  
  pdf(file = filename, width = 12, height = max(8, nrow(prop_matrix) * 0.5))
  
  # 转换为矩阵格式用于热图
  heatmap_data <- as.matrix(prop_matrix)
  
  # 创建颜色映射
  colors <- colorRampPalette(c("#FFFFFF", "#FFF2CC", "#FFE599", "#FFD966", 
                               "#F1C232", "#E69138", "#CC6600", "#B45F06"))(100)
  
  heatmap(
    heatmap_data,
    col = colors,
    main = title,
    xlab = "Group",
    ylab = "Cell Type",
    cexRow = 0.9,
    cexCol = 1.2,
    margins = c(8, 12),
    scale = "none",
    Rowv = NA,
    Colv = NA
  )
  
  dev.off()
  
  message(sprintf("热图已保存: %s", filename))
}

# 生成所有比例图
create_within_group_proportion_plot(
  proportion_results$comprehensive_stats,
  "Group_Proportion_Plots/Within_Group_CellType_Proportions.pdf"
)

create_between_group_proportion_plot(
  proportion_results$comprehensive_stats,
  "Group_Proportion_Plots/Between_Group_CellType_Proportions.pdf"
)

create_absolute_count_plot(
  proportion_results$comprehensive_stats,
  "Group_Proportion_Plots/Absolute_CellType_Counts.pdf"
)

# 生成热图
create_proportion_heatmap(
  proportion_results$prop_within_group,
  "Group_Proportion_Plots/Within_Group_Proportion_Heatmap.pdf",
  "Cell Type Proportions Within Each Group (%)"
)

create_proportion_heatmap(
  proportion_results$prop_within_celltype,
  "Group_Proportion_Plots/Group_Distribution_Within_CellType_Heatmap.pdf",
  "Group Distribution Within Each Cell Type (%)"
)

# 生成统计摘要报告
generate_proportion_summary <- function(stats_df) {
  
  message("生成比例分析摘要报告...")
  
  # 按组统计
  group_summary <- stats_df %>%
    group_by(Group) %>%
    summarise(
      Total_Cells = sum(Cell_Count),
      CellType_Count = n_distinct(CellType),
      Most_Abundant_CellType = CellType[which.max(Percentage_Within_Group)],
      Max_Percentage = max(Percentage_Within_Group),
      .groups = 'drop'
    )
  
  # 按细胞类型统计  
  celltype_summary <- stats_df %>%
    group_by(CellType) %>%
    summarise(
      Total_Cells = sum(Cell_Count),
      Control_Cells = sum(Cell_Count[Group == "CON"]),
      Treat_Cells = sum(Cell_Count[Group == "IDD"]),
      Control_Percentage = round(Control_Cells / Total_Cells * 100, 2),
      Treat_Percentage = round(Treat_Cells / Total_Cells * 100, 2),
      Fold_Change = ifelse(Control_Cells > 0, Treat_Cells / Control_Cells, Inf),
      .groups = 'drop'
    ) %>%
    arrange(desc(Total_Cells))
  
  # 计算显著性差异（简单的比例检验）
  celltype_summary$Significant_Difference <- sapply(1:nrow(celltype_summary), function(i) {
    if (celltype_summary$Control_Cells[i] > 0 && celltype_summary$Treat_Cells[i] > 0) {
      test_result <- prop.test(
        c(celltype_summary$Control_Cells[i], celltype_summary$Treat_Cells[i]),
        c(sum(stats_df$Cell_Count[stats_df$Group == "CON"]),
          sum(stats_df$Cell_Count[stats_df$Group == "IDD"]))
      )
      return(test_result$p.value < 0.05)
    } else {
      return(TRUE)  # 如果某组为0，认为有显著差异
    }
  })
  
  # 保存摘要报告
  write.csv(group_summary, "Proportion_Statistics/Group_Summary.csv", row.names = FALSE)
  write.csv(celltype_summary, "Proportion_Statistics/CellType_Summary.csv", row.names = FALSE)
  
  # 打印摘要信息
  message("\n=== 细胞比例分析摘要 ===")
  message(sprintf("总细胞数: %d", sum(stats_df$Cell_Count)))
  message(sprintf("对照组细胞数: %d", sum(stats_df$Cell_Count[stats_df$Group == "Control"])))
  message(sprintf("实验组细胞数: %d", sum(stats_df$Cell_Count[stats_df$Group == "Treat"])))
  message(sprintf("细胞类型数: %d", length(unique(stats_df$CellType))))
  
  print("=== 各组细胞类型分布 ===")
  print(group_summary)
  
  print("=== 细胞类型在组间的分布情况 ===")
  print(celltype_summary)
  
  # 找出比例差异最大的细胞类型
  max_diff_celltype <- celltype_summary[which.max(abs(celltype_summary$Control_Percentage - celltype_summary$Treat_Percentage)), ]
  message(sprintf("\n比例差异最大的细胞类型: %s", max_diff_celltype$CellType))
  message(sprintf("对照组占比: %.2f%%, 实验组占比: %.2f%%", 
                  max_diff_celltype$Control_Percentage, max_diff_celltype$Treat_Percentage))
  
  return(list(group_summary = group_summary, celltype_summary = celltype_summary))
}

# 生成摘要报告
summary_results <- generate_proportion_summary(proportion_results$comprehensive_stats)

# 生成Excel格式的综合报告（如果有openxlsx包）
if (require(openxlsx, quietly = TRUE)) {
  wb <- createWorkbook()
  
  # 添加工作表
  addWorksheet(wb, "Comprehensive_Stats")
  addWorksheet(wb, "Group_Summary") 
  addWorksheet(wb, "CellType_Summary")
  addWorksheet(wb, "Cross_Table")
  
  # 写入数据
  writeData(wb, "Comprehensive_Stats", proportion_results$comprehensive_stats)
  writeData(wb, "Group_Summary", summary_results$group_summary)
  writeData(wb, "CellType_Summary", summary_results$celltype_summary)
  writeData(wb, "Cross_Table", as.data.frame(proportion_results$cross_table))
  
  # 保存Excel文件
  saveWorkbook(wb, "CellType_Proportion_Analysis/Complete_Proportion_Analysis.xlsx", overwrite = TRUE)
  message("Excel综合报告已保存: Complete_Proportion_Analysis.xlsx")
}

message("\n=== 细胞比例分析完成！ ===")

# ==============================================================================
# 15. showGenes在正常组和实验组的箱线图分析（基于groupComparison.pdf的细胞注释）
# ==============================================================================

message("开始绘制showGenes在正常组和实验组的箱线图（与groupComparison.pdf风格一致）...")

# 创建输出目录
boxplot_dirs <- c("ShowGenes_BoxPlots", 
                  "ShowGenes_Statistics",
                  "ShowGenes_Individual_Plots")

for (dir in boxplot_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(sprintf("创建目录: %s", dir))
  }
}


# 使用Seurat对象中的细胞类型注释（与groupComparison.pdf一致）
get_consistent_cell_annotations <- function(seurat_obj) {
  # 使用当前的Idents，这与groupComparison.pdf中显示的细胞类型一致
  cell_types <- as.character(Idents(seurat_obj))
  names(cell_types) <- colnames(seurat_obj)
  
  message(sprintf("使用的细胞类型: %s", paste(unique(cell_types), collapse = ", ")))
  message(sprintf("细胞类型数量: %d", length(unique(cell_types))))
  
  return(cell_types)
}

# 获取与groupComparison.pdf一致的细胞注释
consistent_cell_ann <- get_consistent_cell_annotations(filtered_seurat)

# 设置与groupComparison.pdf一致的颜色方案
setup_celltype_colors <- function(cell_types) {
  unique_celltypes <- unique(cell_types)
  nColors <- length(unique_celltypes)
  
  # 利用 colorRampPalette 从 RColorBrewer 的 Set3 调色板中扩展出足够数量的颜色
  myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)
  names(myColors) <- unique_celltypes
  
  message(sprintf("为 %d 个细胞类型分配了颜色", nColors))
  return(myColors)
}

# 生成细胞类型颜色
celltype_colors <- setup_celltype_colors(consistent_cell_ann)

# 箱线图生成函数（突出对照组和实验组差异）
create_groupcomparison_style_boxplots <- function(seurat_obj, genes, cell_ann, type_info) {
  
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
  
  message(sprintf("开始分析 %d 个可用基因", length(available_genes)))
  
  # 获取表达数据
  expr_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  
  # 创建数据框
  plot_data <- data.frame(
    Cell = colnames(seurat_obj),
    CellType = cell_ann,
    Group = type_info,
    stringsAsFactors = FALSE
  )
  
  # 添加基因表达数据
  for (gene in available_genes) {
    plot_data[[gene]] <- as.numeric(expr_data[gene, ])
  }
  
  # 统计结果存储
  all_stats <- data.frame()
  
  # 为每个基因生成箱线图
  for (gene in available_genes) {
    message(sprintf("正在处理基因: %s", gene))
    
    # 准备当前基因的数据
    current_data <- plot_data[, c("Cell", "CellType", "Group", gene)]
    colnames(current_data)[4] <- "Expression"
    
    # 计算表达范围，用于设置y轴范围
    expr_values <- current_data$Expression
    y_max <- quantile(expr_values, 0.95, na.rm = TRUE)
    y_min <- 0
    
    # 为每个细胞类型计算统计检验
    cell_types <- unique(current_data$CellType)
    gene_stats <- data.frame()
    
    for (ct in cell_types) {
      ct_data <- current_data[current_data$CellType == ct, ]
      
      # 检查是否有足够的细胞进行统计
      control_cells <- sum(ct_data$Group == "CON")
      treat_cells <- sum(ct_data$Group == "IDD")
      
      if (control_cells >= 3 && treat_cells >= 3) {
        # 执行Wilcoxon秩和检验
        tryCatch({
          test_result <- wilcox.test(
            Expression ~ Group, 
            data = ct_data,
            alternative = "two.sided"
          )
          
          p_value <- test_result$p.value
          
          # 计算平均表达和中位数
          control_mean <- mean(ct_data$Expression[ct_data$Group == "CON"], na.rm = TRUE)
          treat_mean <- mean(ct_data$Expression[ct_data$Group == "IDD"], na.rm = TRUE)
          control_median <- median(ct_data$Expression[ct_data$Group == "CON"], na.rm = TRUE)
          treat_median <- median(ct_data$Expression[ct_data$Group == "IDD"], na.rm = TRUE)
          
          # 计算fold change
          fold_change <- ifelse(control_mean > 0, treat_mean / control_mean, Inf)
          log2_fc <- ifelse(control_mean > 0, log2(treat_mean / control_mean), NA)
          
          # 添加到统计结果
          stat_row <- data.frame(
            Gene = gene,
            CellType = ct,
            Control_Mean = control_mean,
            Treatment_Mean = treat_mean,
            Control_Median = control_median,
            Treatment_Median = treat_median,
            Fold_Change = fold_change,
            Log2_FC = log2_fc,
            P_value = p_value,
            P_adj = p.adjust(p_value, method = "BH"),
            Significant = p_value < 0.05,
            P_symbol = case_when(
              p_value < 0.001 ~ "***",
              p_value < 0.01 ~ "**", 
              p_value < 0.05 ~ "*",
              TRUE ~ ""
            ),
            Control_Cells = control_cells,
            Treatment_Cells = treat_cells
          )
          
          gene_stats <- rbind(gene_stats, stat_row)
          
        }, error = function(e) {
          message(sprintf("统计检验失败 - 基因: %s, 细胞类型: %s, 错误: %s", 
                          gene, ct, e$message))
        })
      } else {
        # 细胞数不足的情况
        stat_row <- data.frame(
          Gene = gene,
          CellType = ct,
          Control_Mean = ifelse(control_cells > 0, 
                                mean(ct_data$Expression[ct_data$Group == "CON"], na.rm = TRUE), 
                                NA),
          Treatment_Mean = ifelse(treat_cells > 0,
                                  mean(ct_data$Expression[ct_data$Group == "IDD"], na.rm = TRUE),
                                  NA),
          Control_Median = ifelse(control_cells > 0,
                                  median(ct_data$Expression[ct_data$Group == "CON"], na.rm = TRUE),
                                  NA),
          Treatment_Median = ifelse(treat_cells > 0,
                                    median(ct_data$Expression[ct_data$Group == "IDD"], na.rm = TRUE),
                                    NA),
          Fold_Change = NA,
          Log2_FC = NA,
          P_value = NA,
          P_adj = NA,
          Significant = FALSE,
          P_symbol = "",
          Control_Cells = control_cells,
          Treatment_Cells = treat_cells
        )
        
        gene_stats <- rbind(gene_stats, stat_row)
      }
    }
    
    # 保存当前基因的统计结果
    write.csv(gene_stats, 
              file.path("ShowGenes_Statistics", paste0(gene, "_boxplot_statistics.csv")),
              row.names = FALSE)
    
    # 合并到总统计结果
    all_stats <- rbind(all_stats, gene_stats)
    
    # 创建箱线图（突出CON vs IDD差异）
    tryCatch({
      
      # 计算细胞类型数量用于调整图形大小
      n_celltypes <- length(unique(current_data$CellType))
      plot_width <- max(12, n_celltypes * 1.5)
      plot_height <- 8
      
      # 创建基础图形，按Group分组显示
      p <- ggplot(current_data, aes(x = CellType, y = Expression, fill = Group)) +
        
        # 箱线图 - 并排显示Control和Treatment
        geom_boxplot(alpha = 0.8, 
                     position = position_dodge(width = 0.8),
                     outlier.size = 0.8,
                     outlier.alpha = 0.6,
                     color = "black",
                     size = 0.5) +
        
        # 添加散点表示单个细胞（可选）
        geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
                    size = 0.3, 
                    alpha = 0.4) +
        
        # 为Control/Treat分组使用对比色（更明显的颜色差异）
        scale_fill_manual(values = c("CON" = "#4CAF50", "IDD" = "#F44336"), 
                          name = "Group") +
        
        # 主题设置
        theme_classic(base_size = 12) +
        theme(
          # 标题样式
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold", 
                                    margin = margin(b = 20)),
          plot.subtitle = element_text(hjust = 0.5, size = 12, 
                                       margin = margin(b = 15), color = "gray50"),
          
          # 坐标轴样式
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, face = "bold", 
                                      margin = margin(t = 10)),
          axis.title.y = element_text(size = 13, face = "bold", 
                                      margin = margin(r = 10)),
          
          # 坐标轴线条
          axis.line = element_line(color = "black", size = 0.8),
          axis.ticks = element_line(color = "black", size = 0.6),
          axis.ticks.length = unit(0.2, "cm"),
          
          # 图例样式
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          legend.position = "top",
          legend.margin = margin(b = 10),
          legend.box.background = element_rect(color = "black", size = 0.5),
          
          # 面板设置
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.y = element_line(color = "gray90", size = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          
          # 边距设置
          plot.margin = margin(20, 20, 20, 20)
        ) +
        
        # 标签设置
        labs(
          title = paste0(gene, " Expression: CON vs IDD"),
          subtitle = paste0("Comparison across cell types (n = ", nrow(current_data), " cells)"),
          x = "Cell Type",
          y = paste0(gene, " Expression Level")
        ) +
        
        # Y轴设置
        scale_y_continuous(
          labels = scales::number_format(accuracy = 0.1),
          expand = expansion(mult = c(0.05, 0.1))
        )
      
      # 添加显著性标记和连接线
      if (nrow(gene_stats) > 0) {
        sig_data <- gene_stats[gene_stats$P_symbol != "", ]
        
        if (nrow(sig_data) > 0) {
          # 计算注释位置
          max_expr_by_celltype <- current_data %>%
            group_by(CellType) %>%
            summarise(max_expr = max(Expression, na.rm = TRUE), .groups = 'drop')
          
          sig_data <- merge(sig_data, max_expr_by_celltype, by = "CellType")
          sig_data$y_pos <- sig_data$max_expr * 1.1
          
          # 添加显著性连接线
          p <- p + 
            geom_segment(data = sig_data,
                         aes(x = as.numeric(as.factor(CellType)) - 0.2,
                             xend = as.numeric(as.factor(CellType)) + 0.2,
                             y = y_pos,
                             yend = y_pos),
                         inherit.aes = FALSE,
                         color = "black", size = 0.8) +
            
            # 添加显著性标记
            geom_text(data = sig_data,
                      aes(x = CellType, y = y_pos + y_max * 0.02, label = P_symbol),
                      inherit.aes = FALSE,
                      size = 5, color = "#D32F2F", fontface = "bold",
                      vjust = 0)
        }
      }
      
      # 保存高质量PDF
      pdf_file <- file.path("ShowGenes_Individual_Plots", paste0(gene, "_BoxPlot_CONVsIDD.pdf"))
      
      pdf(file = pdf_file, width = plot_width, height = plot_height)
      print(p)
      dev.off()
      
      message(sprintf("已保存基因 %s 的箱线图: %s", gene, pdf_file))
      
      # 同时保存高分辨率PNG
      png_file <- file.path("ShowGenes_Individual_Plots", paste0(gene, "_BoxPlot_CONVsIDD.png"))
      ggsave(png_file, plot = p, width = plot_width, height = plot_height, 
             dpi = 300, bg = "white")
      
    }, error = function(e) {
      message(sprintf("生成基因 %s 的箱线图时出错: %s", gene, e$message))
    })
  }
  
  # 保存所有基因的统计结果
  write.csv(all_stats, "ShowGenes_Statistics/All_Genes_BoxPlot_Statistics.csv", row.names = FALSE)
  
  return(all_stats)
}

# 执行分析
showgenes_stats_boxplot <- create_groupcomparison_style_boxplots(
  filtered_seurat, 
  showGenes, 
  consistent_cell_ann, 
  Type
)

# 创建综合箱线图（所有基因在一个图中）
create_comprehensive_boxplots <- function(seurat_obj, genes, cell_ann, type_info, stats_df) {
  
  message("创建综合箱线图...")
  
  # 检查基因是否存在
  available_genes <- genes[genes %in% rownames(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))]
  
  if (length(available_genes) == 0) {
    warning("没有可用基因用于综合图")
    return()
  }
  
  # 获取表达数据
  expr_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  
  # 创建长格式数据
  plot_data_long <- data.frame()
  
  for (gene in available_genes) {
    gene_data <- data.frame(
      Cell = colnames(seurat_obj),
      CellType = cell_ann,
      Group = type_info,
      Gene = gene,
      Expression = as.numeric(expr_data[gene, ]),
      stringsAsFactors = FALSE
    )
    
    plot_data_long <- rbind(plot_data_long, gene_data)
  }
  
  # 创建分面箱线图
  tryCatch({
    p <- ggplot(plot_data_long, aes(x = CellType, y = Expression, fill = Group)) +
      
      # 箱线图
      geom_boxplot(alpha = 0.8, 
                   position = position_dodge(width = 0.8),
                   outlier.size = 0.5,
                   outlier.alpha = 0.5,
                   color = "black",
                   size = 0.4) +
      
      # 分面设置
      facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
      
      # 颜色设置（突出Control vs Treatment差异）
      scale_fill_manual(values = c("CON" = "#4CAF50", "IDD" = "#F44336"), 
                        name = "Group") +
      
      # 主题
      theme_classic(base_size = 10) +
      theme(
        # 主标题
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", 
                                  margin = margin(b = 25)),
        plot.subtitle = element_text(hjust = 0.5, size = 12, 
                                     margin = margin(b = 20), color = "gray50"),
        
        # 坐标轴
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, face = "bold",
                                    margin = margin(t = 10)),
        axis.title.y = element_text(size = 11, face = "bold",
                                    margin = margin(r = 10)),
        
        # 坐标轴线条
        axis.line = element_line(color = "black", size = 0.6),
        axis.ticks = element_line(color = "black", size = 0.4),
        
        # 图例
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "top",
        legend.margin = margin(b = 15),
        legend.box.background = element_rect(color = "black", size = 0.5),
        
        # 分面标签
        strip.background = element_rect(fill = "gray90", color = "black", size = 0.5),
        strip.text = element_text(size = 10, face = "bold"),
        
        # 面板
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(color = "gray90", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        
        # 边距
        plot.margin = margin(25, 25, 25, 25)
      ) +
      
      # 标签
      labs(
        title = "ShowGenes Expression: CON vs IDD Comparison",
        subtitle = "Box plots showing expression differences across cell types",
        x = "Cell Type",
        y = "Gene Expression Level"
      ) +
      
      # Y轴刻度
      scale_y_continuous(
        labels = scales::number_format(accuracy = 0.1),
        expand = expansion(mult = c(0.05, 0.1))
      )
    
    # 添加显著性标记
    if (!is.null(stats_df) && nrow(stats_df) > 0) {
      sig_data <- stats_df[stats_df$P_symbol != "", ]
      
      if (nrow(sig_data) > 0) {
        # 计算每个基因-细胞类型组合的最大表达值
        max_expr_data <- plot_data_long %>%
          group_by(Gene, CellType) %>%
          summarise(max_expr = max(Expression, na.rm = TRUE), .groups = 'drop')
        
        sig_data <- merge(sig_data, max_expr_data, by = c("Gene", "CellType"))
        sig_data$y_pos <- sig_data$max_expr * 1.1
        
        p <- p + 
          geom_text(data = sig_data,
                    aes(x = CellType, y = y_pos, label = P_symbol),
                    inherit.aes = FALSE,
                    size = 3.5, color = "#D32F2F", fontface = "bold",
                    vjust = 0)
      }
    }
    
    # 计算图形大小
    n_genes <- length(available_genes)
    plot_width <- 20
    plot_height <- max(15, ceiling(n_genes / 3) * 6)
    
    # 保存综合图
    pdf_file <- "ShowGenes_BoxPlots/Comprehensive_ShowGenes_BoxPlots.pdf"
    pdf(file = pdf_file, width = plot_width, height = plot_height)
    print(p)
    dev.off()
    
    # 保存PNG版本
    png_file <- "ShowGenes_BoxPlots/Comprehensive_ShowGenes_BoxPlots.png"
    ggsave(png_file, plot = p, width = plot_width, height = plot_height, 
           dpi = 300, bg = "white")
    
    message(sprintf("已保存综合箱线图: %s", pdf_file))
    
  }, error = function(e) {
    message(sprintf("生成综合箱线图时出错: %s", e$message))
  })
}

# 生成综合箱线图
create_comprehensive_boxplots(
  filtered_seurat, 
  showGenes, 
  consistent_cell_ann, 
  Type,
  showgenes_stats_boxplot
)

# 生成统计报告
generate_boxplot_report <- function(stats_df) {
  
  message("生成箱线图统计报告...")
  
  if (is.null(stats_df) || nrow(stats_df) == 0) {
    message("没有统计结果可供分析")
    return()
  }
  
  # 总体摘要
  total_tests <- nrow(stats_df)
  significant_tests <- sum(stats_df$Significant, na.rm = TRUE)
  
  # 按基因汇总
  gene_summary <- stats_df %>%
    filter(!is.na(P_value)) %>%
    group_by(Gene) %>%
    summarise(
      Total_CellTypes_Tested = n(),
      Significant_CellTypes = sum(Significant, na.rm = TRUE),
      Sig_Percentage = round(Significant_CellTypes / Total_CellTypes_Tested * 100, 2),
      Min_Pvalue = min(P_value, na.rm = TRUE),
      Mean_Log2FC = mean(Log2_FC, na.rm = TRUE),
      Max_AbsLog2FC = max(abs(Log2_FC), na.rm = TRUE),
      Upregulated_CellTypes = sum(Log2_FC > 0, na.rm = TRUE),
      Downregulated_CellTypes = sum(Log2_FC < 0, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(Significant_CellTypes), Min_Pvalue)
  
  # 按细胞类型汇总
  celltype_summary <- stats_df %>%
    filter(!is.na(P_value)) %>%
    group_by(CellType) %>%
    summarise(
      Total_Genes_Tested = n(),
      Significant_Genes = sum(Significant, na.rm = TRUE),
      Sig_Percentage = round(Significant_Genes / Total_Genes_Tested * 100, 2),
      Min_Pvalue = min(P_value, na.rm = TRUE),
      Mean_Log2FC = mean(Log2_FC, na.rm = TRUE),
      Upregulated_Genes = sum(Log2_FC > 0, na.rm = TRUE),
      Downregulated_Genes = sum(Log2_FC < 0, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(Significant_Genes))
  
  # 显著差异的基因-细胞类型组合
  significant_combinations <- stats_df %>%
    filter(Significant == TRUE) %>%
    select(Gene, CellType, Log2_FC, P_value, P_adj, Control_Mean, Treatment_Mean, 
           Control_Median, Treatment_Median, Fold_Change, P_symbol) %>%
    arrange(P_value)
  
  # 保存详细报告
  write.csv(gene_summary, "ShowGenes_Statistics/Gene_Summary_BoxPlot.csv", row.names = FALSE)
  write.csv(celltype_summary, "ShowGenes_Statistics/CellType_Summary_BoxPlot.csv", row.names = FALSE)
  write.csv(significant_combinations, "ShowGenes_Statistics/Significant_Combinations_BoxPlot.csv", row.names = FALSE)
  
  # 控制台输出摘要
  message("\n=== ShowGenes 箱线图分析摘要 ===")
  message(sprintf("总检验次数: %d", total_tests))
  message(sprintf("显著检验次数: %d (%.2f%%)", significant_tests, significant_tests/total_tests*100))
  message(sprintf("分析的基因数: %d", nrow(gene_summary)))
  message(sprintf("分析的细胞类型数: %d", nrow(celltype_summary)))
  
  # 显示上调和下调情况
  if (nrow(significant_combinations) > 0) {
    upregulated_count <- sum(significant_combinations$Log2_FC > 0, na.rm = TRUE)
    downregulated_count <- sum(significant_combinations$Log2_FC < 0, na.rm = TRUE)
    message(sprintf("上调的基因-细胞类型组合: %d", upregulated_count))
    message(sprintf("下调的基因-细胞类型组合: %d", downregulated_count))
    
    most_sig <- significant_combinations[1, ]
    direction <- ifelse(most_sig$Log2_FC > 0, "上调", "下调")
    message(sprintf("\n最显著的差异: %s 基因在 %s 细胞中%s", 
                    most_sig$Gene, most_sig$CellType, direction))
    message(sprintf("Log2FC: %.3f, P-value: %.2e", most_sig$Log2_FC, most_sig$P_value))
  }
  
  return(list(
    gene_summary = gene_summary, 
    celltype_summary = celltype_summary,
    significant_combinations = significant_combinations
  ))
}

# 生成报告
boxplot_report <- generate_boxplot_report(showgenes_stats_boxplot)

message("\n=== ShowGenes 箱线图分析完成！ ===")


# ==============================================================================
# 16. showGenes在UMAP中的表达分布可视化（对照组vs实验组并排比较）
# ==============================================================================

message("开始绘制showGenes在UMAP中的表达分布图（对照组vs实验组并排比较）...")

# 创建输出目录
umap_dirs <- c(  "ShowGenes_UMAP_Individual",
               "ShowGenes_UMAP_Statistics")

for (dir in umap_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(sprintf("创建目录: %s", dir))
  }
}


# 生成对照组vs实验组并排比较的UMAP基因表达图
create_split_gene_umap_plots <- function(seurat_obj, genes, group_by = "Type") {
  
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
  
  message(sprintf("开始生成 %d 个可用基因的UMAP表达图（对照组vs实验组比较）", length(available_genes)))
  
  # 为每个基因单独生成对照组vs实验组的并排比较图
  for (gene in available_genes) {
    message(sprintf("正在处理基因: %s", gene))
    
    tryCatch({
      
      # 方法1：使用split.by参数生成并排比较图（推荐方法）
      p_split <- FeaturePlot(
        object = seurat_obj, 
        features = gene, 
        cols = c("lightgrey", "#FF6B6B", "#FF0000", "#8B0000"),  # 灰色到深红色渐变
        reduction = "umap",
        pt.size = 1.2,
        split.by = group_by,  # 按照Type分组（Control vs Treat）
        order = TRUE,  # 将高表达细胞绘制在前面
        combine = TRUE  # 确保图形组合在一起
      ) + 
        plot_annotation(
          title = paste0(gene, " Expression: Control vs Treatment Comparison"),
          subtitle = "UMAP visualization showing gene expression levels across groups",
          theme = theme(
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 10)),
            plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 20), color = "gray50")
          )
        ) &
        theme_minimal() &
        theme(
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 14, face = "bold", color = "black"),
          strip.background = element_rect(fill = "gray90", color = "black", size = 0.8),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.margin = margin(15, 15, 15, 15),
          legend.position = "bottom",
          legend.box = "horizontal"
        ) &
        guides(
          colour = guide_colorbar(
            title = "Expression Level",
            title.position = "top",
            title.hjust = 0.5,
            barwidth = 15,
            barheight = 1.2
          )
        )
      
      # 保存PDF格式的分组比较图
      pdf_file <- file.path("ShowGenes_UMAP_Individual", 
                            paste0(gene, "_UMAP_CON_vs_IDD.pdf"))
      
      ggsave(
        filename = pdf_file,
        plot = p_split,
        width = 14,
        height = 7,
        dpi = 300,
        device = "pdf"
      )
      
      message(sprintf("已保存基因 %s 的CON组vsIDD组UMAP表达图: %s", gene, pdf_file))
      
      # 方法2：生成带有更多统计信息的版本
      # 首先计算每组的表达统计
      expr_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
      gene_expr <- as.numeric(expr_data[gene, ])
      
      control_expr <- gene_expr[seurat_obj$Type == "CON"]
      treat_expr <- gene_expr[seurat_obj$Type == "IDD"]
      
      control_mean <- round(mean(control_expr, na.rm = TRUE), 3)
      treat_mean <- round(mean(treat_expr, na.rm = TRUE), 3)
      control_pct <- round(sum(control_expr > 0) / length(control_expr) * 100, 1)
      treat_pct <- round(sum(treat_expr > 0) / length(treat_expr) * 100, 1)
      
      # 带统计信息的标题
      detailed_subtitle <- sprintf(
        "CON: Mean=%.3f, Expressing=%.1f%% | IDD: Mean=%.3f, Expressing=%.1f%%",
        control_mean, control_pct, treat_mean, treat_pct
      )
      
      p_split_detailed <- FeaturePlot(
        object = seurat_obj, 
        features = gene, 
        cols = viridis::viridis(100),  # 使用viridis调色板
        reduction = "umap",
        pt.size = 1.2,
        split.by = group_by,
        order = TRUE,
        combine = TRUE
      ) + 
        plot_annotation(
          title = paste0(gene, " Expression Distribution"),
          subtitle = detailed_subtitle,
          theme = theme(
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 5)),
            plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 15), color = "gray30")
          )
        ) &
        theme_minimal() &
        theme(
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 14, face = "bold", color = "white"),
          strip.background = element_rect(fill = "gray20", color = "black", size = 0.8),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.margin = margin(15, 15, 15, 15),
          legend.position = "bottom"
        ) &
        guides(
          colour = guide_colorbar(
            title = "Expression Level",
            title.position = "top",
            title.hjust = 0.5,
            barwidth = 15,
            barheight = 1.2
          )
        )
      
      # 保存带详细统计信息的版本
      pdf_file_detailed <- file.path("ShowGenes_UMAP_Individual", 
                                     paste0(gene, "_UMAP_Detailed_Stats.pdf"))
      
      ggsave(
        filename = pdf_file_detailed,
        plot = p_split_detailed,
        width = 14,
        height = 7,
        dpi = 300,
        device = "pdf"
      )
      
    }, error = function(e) {
      message(sprintf("生成基因 %s 的UMAP表达图时出错: %s", gene, e$message))
    })
  }
  
  return(available_genes)
}

# 生成所有基因的综合比较图
create_comprehensive_split_umap <- function(seurat_obj, genes, group_by = "Type") {
  
  message("生成所有基因的综合对照组vs实验组UMAP比较图...")
  
  # 检查基因是否存在
  available_genes <- genes[genes %in% rownames(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))]
  
  if (length(available_genes) == 0) {
    warning("没有可用基因用于综合图")
    return()
  }
  
  tryCatch({
    
    # 方法1：使用FeaturePlot生成所有基因的综合分组图
    n_genes <- length(available_genes)
    ncol_genes <- min(3, n_genes)  # 每行最多3个基因
    
    p_comprehensive <- FeaturePlot(
      object = seurat_obj, 
      features = available_genes, 
      cols = c("lightgrey", "#FF6B6B", "#FF0000", "#8B0000"),
      reduction = "umap",
      pt.size = 0.8,
      split.by = group_by,
      ncol = ncol_genes,  # 控制基因排列的列数
      order = TRUE,
      combine = TRUE
    ) + 
      plot_annotation(
        title = "ShowGenes Expression: Comprehensive CON vs IDD Comparison",
        subtitle = paste0("UMAP visualization of ", n_genes, " genes across experimental groups"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
          plot.subtitle = element_text(hjust = 0.5, size = 16, margin = margin(b = 20), color = "gray50")
        )
      ) &
      theme_minimal() &
      theme(
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 11, face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray90", color = "black", size = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        plot.margin = margin(10, 10, 10, 10),
        legend.position = "bottom"
      ) &
      guides(
        colour = guide_colorbar(
          title = "Expression Level",
          title.position = "top",
          title.hjust = 0.5,
          barwidth = 20,
          barheight = 1
        )
      )
    
    # 计算图形大小
    nrow_genes <- ceiling(n_genes / ncol_genes)
    plot_width <- max(18, ncol_genes * 6)
    plot_height <- max(12, nrow_genes * 4 + 4)
    
    # 保存综合比较图
   # comprehensive_file <- file.path("ShowGenes_UMAP_Expression", "Comprehensive_ShowGenes_Control_vs_Treatment.pdf")
    
    ggsave(
      filename = comprehensive_file,
      plot = p_comprehensive,
      width = plot_width,
      height = plot_height,
      dpi = 300,
      device = "pdf"
    )
    
    message(sprintf("已保存综合CON组vsIDD组UMAP比较图: %s", comprehensive_file))
    message(sprintf("图形尺寸: %.0f x %.0f", plot_width, plot_height))
    
  }, error = function(e) {
    message(sprintf("生成综合UMAP比较图时出错: %s", e$message))
  })
}

# 生成基因表达的定量比较统计
create_expression_statistics <- function(seurat_obj, genes, type_col = "Type") {
  
  message("生成基因表达的定量统计比较...")
  
  # 检查基因是否存在
  available_genes <- genes[genes %in% rownames(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))]
  
  if (length(available_genes) == 0) return()
  
  # 获取表达数据
  expr_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  
  # 创建统计表
  stats_df <- data.frame()
  
  for (gene in available_genes) {
    gene_expr <- as.numeric(expr_data[gene, ])
    
    # 分组数据
    control_expr <- gene_expr[seurat_obj@meta.data[[type_col]] == "CON"]
    treat_expr <- gene_expr[seurat_obj@meta.data[[type_col]] == "IDD"]
    
    # 计算统计量
    stats <- data.frame(
      Gene = gene,
      Control_Mean = mean(control_expr, na.rm = TRUE),
      Control_Median = median(control_expr, na.rm = TRUE),
      Control_SD = sd(control_expr, na.rm = TRUE),
      Control_Expressing_Cells = sum(control_expr > 0),
      Control_Total_Cells = length(control_expr),
      Control_Expression_Pct = round(sum(control_expr > 0) / length(control_expr) * 100, 2),
      
      Treatment_Mean = mean(treat_expr, na.rm = TRUE),
      Treatment_Median = median(treat_expr, na.rm = TRUE),
      Treatment_SD = sd(treat_expr, na.rm = TRUE),
      Treatment_Expressing_Cells = sum(treat_expr > 0),
      Treatment_Total_Cells = length(treat_expr),
      Treatment_Expression_Pct = round(sum(treat_expr > 0) / length(treat_expr) * 100, 2)
    )
    
    # 计算fold change和统计检验
    stats$Fold_Change_Mean = ifelse(stats$Control_Mean > 0, 
                                    stats$Treatment_Mean / stats$Control_Mean, 
                                    Inf)
    stats$Log2_FC = ifelse(stats$Control_Mean > 0, 
                           log2(stats$Treatment_Mean / stats$Control_Mean), 
                           NA)
    
    # Wilcoxon检验
    if (length(control_expr) >= 3 && length(treat_expr) >= 3) {
      tryCatch({
        wilcox_test <- wilcox.test(treat_expr, control_expr, alternative = "two.sided")
        stats$P_value = wilcox_test$p.value
        stats$Significant = wilcox_test$p.value < 0.05
      }, error = function(e) {
        stats$P_value = NA
        stats$Significant = FALSE
      })
    } else {
      stats$P_value = NA
      stats$Significant = FALSE
    }
    
    stats_df <- rbind(stats_df, stats)
  }
  
  # 调整p值
  if (any(!is.na(stats_df$P_value))) {
    stats_df$P_adj = p.adjust(stats_df$P_value, method = "BH")
  } else {
    stats_df$P_adj = NA
  }
  
  # 保存统计结果
  write.csv(stats_df, "ShowGenes_UMAP_Statistics/Gene_Expression_Group_Comparison.csv", row.names = FALSE)
  
  # 输出摘要到控制台
  message("\n=== ShowGenes 表达比较统计摘要 ===")
  for (i in 1:nrow(stats_df)) {
    gene_stats <- stats_df[i, ]
    direction <- ifelse(is.na(gene_stats$Log2_FC), "无变化", 
                        ifelse(gene_stats$Log2_FC > 0, "上调", "下调"))
    
    message(sprintf("%s: %s (Log2FC=%.3f), CON %.1f%% vs IDD %.1f%%, P=%.3f",
                    gene_stats$Gene,
                    direction,
                    ifelse(is.na(gene_stats$Log2_FC), 0, gene_stats$Log2_FC),
                    gene_stats$Control_Expression_Pct,
                    gene_stats$Treatment_Expression_Pct,
                    ifelse(is.na(gene_stats$P_value), 1, gene_stats$P_value)))
  }
  
  return(stats_df)
}

# 执行所有UMAP表达分析
message("开始执行showGenes的UMAP表达分析（CON组vsIDD组比较）...")

# 生成单个基因的对照组vs实验组比较图
available_showgenes <- create_split_gene_umap_plots(
  filtered_seurat, 
  showGenes, 
  group_by = "Type"
)

# 生成综合比较图
create_comprehensive_split_umap(
  filtered_seurat, 
  available_showgenes, 
  group_by = "Type"
)

# 生成定量统计比较
expression_stats <- create_expression_statistics(
  filtered_seurat, 
  showGenes, 
  type_col = "Type"
)

message("\n=== ShowGenes UMAP表达分析完成！ ===")

# ==============================================================================
# 17. ShowGenes基因和注释单细胞的气泡图（Dot Plot）- 分别保存版本
# ==============================================================================

message("开始绘制ShowGenes基因和注释单细胞的气泡图...")

# 创建输出目录
dotplot_dirs <- c("ShowGenes_DotPlots", 
                  "ShowGenes_DotPlot_Data")

for (dir in dotplot_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(sprintf("创建目录: %s", dir))
  }
}

# 生成ShowGenes的气泡图函数
create_showgenes_dotplot <- function(seurat_obj, genes, group_by = NULL) {
  
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
  
  message(sprintf("开始生成 %d 个可用基因的气泡图", length(available_genes)))
  
  # 方法1：基础气泡图（按细胞类型）
  create_basic_dotplot <- function() {
    message("生成基础气泡图（按细胞类型）...")
    
    tryCatch({
      p1 <- DotPlot(
        object = seurat_obj, 
        features = available_genes,
        cols = c("lightgrey", "#FF0000"),  # 灰色到红色渐变
        dot.scale = 8,  # 控制点的大小
        split.by = NULL  # 不分组
      ) + 
        theme_minimal() +
        theme(
          # 标题设置
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", 
                                    margin = margin(b = 20)),
          plot.subtitle = element_text(hjust = 0.5, size = 14, 
                                       margin = margin(b = 15), color = "gray50"),
          
          # 坐标轴设置
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, 
                                     face = "bold.italic", color = "black"),
          axis.text.y = element_text(size = 12, face = "bold", color = "black"),
          axis.title.x = element_text(size = 14, face = "bold", 
                                      margin = margin(t = 10)),
          axis.title.y = element_text(size = 14, face = "bold", 
                                      margin = margin(r = 10)),
          
          # 图例设置
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.margin = margin(t = 15),
          
          # 面板设置
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "gray90", size = 0.5),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          
          # 边距设置
          plot.margin = margin(25, 25, 25, 25)
        ) +
        
        # 标签设置
        labs(
          title = "ShowGenes Expression: Dot Plot Analysis",
          subtitle = paste0("Gene expression patterns across cell types (", 
                            length(available_genes), " genes)"),
          x = "Genes",
          y = "Cell Types",
          caption = "Dot size = % of cells expressing gene; Color intensity = average expression level"
        ) +
        
        # 颜色渐变设置
        scale_color_gradient2(
          low = "lightgrey", 
          mid = "#FF6B6B", 
          high = "#8B0000",
          midpoint = 0.5,
          name = "Average\nExpression",
          guide = guide_colorbar(
            title.position = "top",
            title.hjust = 0.5,
            barwidth = 8,
            barheight = 1
          )
        ) +
        
        # 点大小设置
        scale_size_continuous(
          name = "Percent\nExpressed",
          guide = guide_legend(
            title.position = "top",
            title.hjust = 0.5,
            override.aes = list(color = "black")
          )
        ) +
        
        # 坐标轴翻转（可选）
        coord_flip()
      
      return(p1)
      
    }, error = function(e) {
      message(sprintf("生成基础气泡图时出错: %s", e$message))
      return(NULL)
    })
  }
  
  # 方法2：分组气泡图（按Treatment vs Control）
  create_grouped_dotplot <- function() {
    message("生成分组气泡图（CON vs IDD）...")
    
    tryCatch({
      p2 <- DotPlot(
        object = seurat_obj, 
        features = available_genes,
        cols = viridis::viridis(100),  # 使用viridis调色板
        dot.scale = 6,
        split.by = "Type"  # 按照Type分组
      ) + 
        theme_minimal() +
        theme(
          # 标题设置
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", 
                                    margin = margin(b = 20)),
          plot.subtitle = element_text(hjust = 0.5, size = 14, 
                                       margin = margin(b = 15), color = "gray50"),
          
          # 坐标轴设置
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11, 
                                     face = "bold.italic", color = "black"),
          axis.text.y = element_text(size = 11, face = "bold", color = "black"),
          axis.title.x = element_text(size = 14, face = "bold", 
                                      margin = margin(t = 10)),
          axis.title.y = element_text(size = 14, face = "bold", 
                                      margin = margin(r = 10)),
          
          # 图例设置
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.margin = margin(t = 15),
          
          # 面板设置
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "gray90", size = 0.5),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          
          # 分面标签设置
          strip.text = element_text(size = 12, face = "bold", color = "white"),
          strip.background = element_rect(fill = "gray20", color = "black", size = 0.8),
          
          # 边距设置
          plot.margin = margin(25, 25, 25, 25)
        ) +
        
        # 标签设置
        labs(
          title = "ShowGenes Expression: CON vs IDD Comparison",
          subtitle = paste0("Dot plot showing differential expression patterns (", 
                            length(available_genes), " genes)"),
          x = "Genes",
          y = "Cell Types",
          caption = "Split by experimental groups | Dot size = % expressing; Color = avg expression"
        ) +
        
        # 坐标轴翻转
        coord_flip()
      
      return(p2)
      
    }, error = function(e) {
      message(sprintf("生成分组气泡图时出错: %s", e$message))
      return(NULL)
    })
  }
  
  # 方法3：高级定制气泡图
  create_advanced_dotplot <- function() {
    message("生成高级定制气泡图...")
    
    tryCatch({
      # 首先获取细胞类型信息并排序
      cell_types <- levels(Idents(seurat_obj))
      
      p3 <- DotPlot(
        object = seurat_obj, 
        features = available_genes,
        cols = c("#F0F0F0", "#FFE5E5", "#FFB5B5", "#FF8585", "#FF4444", "#CC0000"),
        dot.scale = 10,
        dot.min = 0.05,  # 最小点大小
        col.min = 0,     # 最小颜色值
        col.max = 3      # 最大颜色值（可调整）
      ) + 
        theme_classic() +
        theme(
          # 标题设置
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold", 
                                    margin = margin(b = 25)),
          plot.subtitle = element_text(hjust = 0.5, size = 16, 
                                       margin = margin(b = 20), color = "gray40"),
          
          # 坐标轴设置
          axis.text.x = element_text(angle = 45, hjust = 1, size = 13, 
                                     face = "bold.italic", color = "black"),
          axis.text.y = element_text(size = 13, face = "bold", color = "black"),
          axis.title.x = element_text(size = 16, face = "bold", 
                                      margin = margin(t = 15)),
          axis.title.y = element_text(size = 16, face = "bold", 
                                      margin = margin(r = 15)),
          
          # 坐标轴线条
          axis.line = element_line(color = "black", size = 1),
          axis.ticks = element_line(color = "black", size = 0.8),
          axis.ticks.length = unit(0.25, "cm"),
          
          # 图例设置
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.margin = margin(t = 20),
          legend.spacing.x = unit(1, "cm"),
          
          # 背景设置
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "gray95", size = 0.8),
          panel.grid.minor = element_line(color = "gray98", size = 0.5),
          
          # 边距设置
          plot.margin = margin(30, 30, 30, 30)
        ) +
        
        # 标签设置
        labs(
          title = "ShowGenes Expression Profile",
          subtitle = paste0("Comprehensive dot plot analysis across cell types\n",
                            "Genes analyzed: ", paste(available_genes, collapse = ", ")),
          x = "Target Genes",
          y = "Annotated Cell Types",
          caption = paste0("Analysis includes ", length(available_genes), 
                           " genes across ", length(cell_types), " cell types\n",
                           "Dot size represents percentage of expressing cells; ",
                           "Color intensity represents average expression level")
        ) +
        
        # 颜色和大小的图例设置
        guides(
          color = guide_colorbar(
            title = "Average Expression",
            title.position = "top",
            title.hjust = 0.5,
            barwidth = 12,
            barheight = 1.5,
            frame.colour = "black",
            ticks.colour = "black"
          ),
          size = guide_legend(
            title = "Percent Expressed",
            title.position = "top",
            title.hjust = 0.5,
            override.aes = list(color = "black"),
            keywidth = 1.2,
            keyheight = 1.2
          )
        )
      
      return(p3)
      
    }, error = function(e) {
      message(sprintf("生成高级定制气泡图时出错: %s", e$message))
      return(NULL)
    })
  }
  
  # 生成所有类型的气泡图
  plots_list <- list()
  
  # 基础气泡图
  p1 <- create_basic_dotplot()
  if (!is.null(p1)) {
    plots_list[["basic"]] <- p1
  }
  
  # 分组气泡图
  p2 <- create_grouped_dotplot()
  if (!is.null(p2)) {
    plots_list[["grouped"]] <- p2
  }
  
  # 高级定制气泡图
  p3 <- create_advanced_dotplot()
  if (!is.null(p3)) {
    plots_list[["advanced"]] <- p3
  }
  
  return(plots_list)
}

# 保存气泡图数据的函数
save_dotplot_data <- function(seurat_obj, genes) {
  
  message("保存气泡图的原始数据...")
  
  # 检查基因是否存在
  available_genes <- genes[genes %in% rownames(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))]
  
  if (length(available_genes) == 0) return()
  
  # 获取DotPlot的数据
  tryCatch({
    dot_data <- DotPlot(seurat_obj, features = available_genes)$data
    
    # 添加额外的信息
    dot_data$Gene <- factor(dot_data$features.plot, levels = available_genes)
    dot_data$CellType <- dot_data$id
    
    # 重新排列列
    dot_data_clean <- dot_data[, c("Gene", "CellType", "avg.exp", "pct.exp", 
                                   "avg.exp.scaled", "features.plot", "id")]
    
    colnames(dot_data_clean) <- c("Gene", "CellType", "Average_Expression", 
                                  "Percent_Expressed", "Scaled_Average_Expression",
                                  "Original_Feature", "Original_ID")
    
    # 保存数据
    write.csv(dot_data_clean, "ShowGenes_DotPlot_Data/DotPlot_Data.csv", row.names = FALSE)
    
    # 生成数据摘要
    summary_stats <- dot_data_clean %>%
      group_by(Gene) %>%
      summarise(
        Max_Expression = max(Average_Expression, na.rm = TRUE),
        Mean_Expression = mean(Average_Expression, na.rm = TRUE),
        Max_Percent = max(Percent_Expressed, na.rm = TRUE),
        Mean_Percent = mean(Percent_Expressed, na.rm = TRUE),
        Expressing_CellTypes = sum(Percent_Expressed > 0),
        Total_CellTypes = n(),
        .groups = 'drop'
      )
    
    write.csv(summary_stats, "ShowGenes_DotPlot_Data/Gene_Summary_Stats.csv", row.names = FALSE)
    
    message("气泡图数据保存完成")
    
    return(dot_data_clean)
    
  }, error = function(e) {
    message(sprintf("保存气泡图数据时出错: %s", e$message))
    return(NULL)
  })
}

# 执行气泡图分析
message("开始执行ShowGenes的气泡图分析...")

# 生成所有类型的气泡图
dotplot_results <- create_showgenes_dotplot(filtered_seurat, showGenes)

# 分别保存每个气泡图到独立的PDF文件
if (length(dotplot_results) > 0) {
  
  # 保存基础气泡图
  if ("basic" %in% names(dotplot_results)) {
    basic_pdf <- "ShowGenes_DotPlots/ShowGenes_DotPlot_Basic.pdf"
    ggsave(
      filename = basic_pdf,
      plot = dotplot_results[["basic"]],
      width = 14,
      height = 10,
      dpi = 300,
      device = "pdf"
    )
    message(sprintf("已保存基础气泡图到: %s", basic_pdf))
  }
  
  # 保存分组气泡图
  if ("grouped" %in% names(dotplot_results)) {
    grouped_pdf <- "ShowGenes_DotPlots/ShowGenes_DotPlot_Grouped.pdf"
    ggsave(
      filename = grouped_pdf,
      plot = dotplot_results[["grouped"]],
      width = 16,
      height = 10,
      dpi = 300,
      device = "pdf"
    )
    message(sprintf("已保存分组气泡图到: %s", grouped_pdf))
  }
  
  # 保存高级定制气泡图
  if ("advanced" %in% names(dotplot_results)) {
    advanced_pdf <- "ShowGenes_DotPlots/ShowGenes_DotPlot_Advanced.pdf"
    ggsave(
      filename = advanced_pdf,
      plot = dotplot_results[["advanced"]],
      width = 14,
      height = 12,
      dpi = 300,
      device = "pdf"
    )
    message(sprintf("已保存高级定制气泡图到: %s", advanced_pdf))
  }
  
  # 可选：也保存一个综合版本（如果需要的话）
  comprehensive_pdf <- "ShowGenes_DotPlots/ShowGenes_DotPlot_All_Combined.pdf"
  
  pdf(file = comprehensive_pdf, width = 14, height = 10)
  
  if ("basic" %in% names(dotplot_results)) {
    print(dotplot_results[["basic"]])
  }
  
  if ("grouped" %in% names(dotplot_results)) {
    print(dotplot_results[["grouped"]])
  }
  
  if ("advanced" %in% names(dotplot_results)) {
    print(dotplot_results[["advanced"]])
  }
  
  dev.off()
  
  message(sprintf("已保存综合版气泡图到: %s", comprehensive_pdf))
}

# 保存气泡图数据
dotplot_data <- save_dotplot_data(filtered_seurat, showGenes)

# 生成简单的统计报告
if (!is.null(dotplot_data)) {
  message("\n=== ShowGenes 气泡图分析摘要 ===")
  
  # 按基因统计
  gene_stats <- dotplot_data %>%
    group_by(Gene) %>%
    summarise(
      Max_Expression = max(Average_Expression, na.rm = TRUE),
      Max_Percent = max(Percent_Expressed, na.rm = TRUE),
      Expressing_CellTypes = sum(Percent_Expressed > 5),  # 表达>5%的细胞类型
      .groups = 'drop'
    ) %>%
    arrange(desc(Max_Expression))
  
  for (i in 1:nrow(gene_stats)) {
    gene_info <- gene_stats[i, ]
    message(sprintf("%s: 最高表达 %.3f, 最高表达细胞比例 %.1f%%, 表达细胞类型数 %d",
                    gene_info$Gene,
                    gene_info$Max_Expression,
                    gene_info$Max_Percent,
                    gene_info$Expressing_CellTypes))
  }
  
  # 按细胞类型统计
  celltype_stats <- dotplot_data %>%
    group_by(CellType) %>%
    summarise(
      Expressing_Genes = sum(Percent_Expressed > 5),
      Mean_Expression = mean(Average_Expression, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(Expressing_Genes))
  
  message(sprintf("\n细胞类型表达概况（前5名）:"))
  for (i in 1:min(5, nrow(celltype_stats))) {
    ct_info <- celltype_stats[i, ]
    message(sprintf("%s: 表达基因数 %d, 平均表达水平 %.3f",
                    ct_info$CellType,
                    ct_info$Expressing_Genes,
                    ct_info$Mean_Expression))
  }
}

message("\n=== ShowGenes 气泡图分析完成！ ===")

###########################################################
##########################################################
# 保存Seurat对象
save(filtered_seurat, cellAnn,  analysis_params, annotation_results, marker_results, 
     file = "filtered_seurat_2.rda")
load(file = "filtered_seurat_2.rda")
######################################################
