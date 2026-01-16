# ========================== 0. 设置工作目录（可根据实际情况修改） ==========================
setwd("I:\\2025教程新n\\14.Q1区7.8分，MR+单细胞\\05.对照组和实验组的差异分析")  # 


threshold_logFC <- 1       # log2折叠变化阈值
threshold_adjP  <- 0.05    # 调整后P值阈值
max_display_genes <- 50    # 热图中每侧（上调和下调）展示的基因最大数量
result_dir <- "result_csv"
if(!dir.exists(result_dir)) dir.create(result_dir)

# 记录 Step 1 开始时间
step1_start <- Sys.time()
cat("[Step 1] 开始运行...\n")

# 载入必要的 R 包 -----------------------------------------------------
library(limma)    # 用于线性模型和差异表达分析

# 检查并安装 pheatmap 包 -----------------------------------------------------
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)
library(ggrepel)
library(pheatmap) # 用于绘制热图

# 检查并安装 RColorBrewer 包（用于更丰富的颜色方案）---------------------------
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)

# 记录 Step 1 结束时间
step1_end <- Sys.time()
cat("[Step 1] 运行完毕，耗时:", 
    round(difftime(step1_end, step1_start, units="secs"), 2), "秒\n\n")



# ========================== 2. 参数设定及文件路径配置 ==========================
step2_start <- Sys.time()
cat("[Step 2] 开始参数设定和文件路径配置...\n")



# 定义输入文件路径（请根据实际情况修改）---------------------------------------
file_expr  <- "Sample Type Matrix.csv"  # 基因表达矩阵文件
 

step2_end <- Sys.time()
cat("[Step 2] 运行完毕，耗时:", 
    round(difftime(step2_end, step2_start, units="secs"), 2), "秒\n\n")
# ========================== 读取表达矩阵并保证数值型 ==========================
file_expr <- "Sample Type Matrix.csv" # 修改成你的实际数据文件名

# 建议自动检测分隔符
tmp_head <- readLines(file_expr, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")

# 读取原始表达文件，行为基因名，列为样本
expr_raw <- read.table(file_expr, header=TRUE, sep=sep, check.names=FALSE, stringsAsFactors=FALSE)
rownames(expr_raw) <- expr_raw[,1]
expr_mat <- expr_raw[,-1, drop=FALSE]

# 强制转为matrix并确保是数值型
expr_mat <- as.matrix(expr_mat)
expr_mat <- apply(expr_mat, 2, as.numeric) # 逐列强制类型
rownames(expr_mat) <- rownames(expr_raw)
colnames(expr_mat) <- colnames(expr_raw)[-1] # 保持列名
combined_expr <- expr_mat

# 检查所有元素皆为数值
if(!is.numeric(expr_mat)) stop("表达矩阵仍存在非数值列！")
cat("表达矩阵维度:", dim(expr_mat)[1], "基因 x", dim(expr_mat)[2], "样本\n")

# ========================== 自动判断分组信息 ==========================
sample_names <- colnames(expr_mat)
group_info <- ifelse(grepl("_con$", sample_names, ignore.case = TRUE), "Control",
                     ifelse(grepl("_tre$", sample_names, ignore.case = TRUE), "Treatment", "Unknown"))
if(any(group_info == "Unknown")) stop("分组未知的样本存在，请检查样本名后缀！")
num_ctrl <- sum(group_info == "Control")
num_treat <- sum(group_info == "Treatment")
cat("分组情况：Control =", num_ctrl, ", Treatment =", num_treat, "\n")

# ========================== 差异表达分析（LIMMA） ==========================
library(limma)
group_labels <- factor(group_info, levels = c("Control", "Treatment"))
design_mat <- model.matrix(~0 + group_labels)
colnames(design_mat) <- c("Control", "Treatment")
fit <- lmFit(expr_mat, design_mat)
contrast_mat <- makeContrasts(Treatment - Control, levels = design_mat)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)
all_diff_results <- topTable(fit2, adjust.method = "fdr", number = Inf)
cat("差异分析完成！发现FDR<0.05基因数：", sum(all_diff_results$adj.P.Val < 0.05), "\n")

# 可选：输出全部基因的差异分析表
write.table(cbind(Gene=rownames(all_diff_results), all_diff_results),
            file="DE_results.csv", sep=",", quote=FALSE, row.names=FALSE)


# 可选：写出到文件
# write.csv(all_diff_results, file="DEG_results.csv")

step5_end <- Sys.time()
#cat("[Step 5] 用时：", round(as.numeric(difftime(step5_end, step5_start, units = "secs")), 2), "秒\n")


# ------------------------- Step 6: 筛选显著差异表达基因并输出包含必要信息的表格 -------------------------
step6_start <- Sys.time()
cat("[Step 6] 开始筛选显著差异表达基因并输出包含必要信息的表格...\n")

# 根据设定的阈值筛选显著上调或下调的基因
significant_DEGs <- all_diff_results[with(all_diff_results,
                                          (abs(logFC) > threshold_logFC & adj.P.Val < threshold_adjP)), ]

# 整理输出格式，添加 "Gene" 列显示基因名称，
# 此外默认的 topTable 输出已经包含以下信息：
#   - logFC: 基因的对数折叠变化
#   - AveExpr: 各样本的平均表达值
#   - t: 经过经验贝叶斯调整后的 t 统计量
#   - P.Value: 原始 p 值
#   - adj.P.Val: 调整后 p 值（FDR校正）
#   - B: 差异表达的对数 odds
output_DEGs <- cbind(Gene = rownames(significant_DEGs), significant_DEGs)

# 计算标准误差 (SE)
# 根据 t 统计量公式： t = logFC / SE  =>  SE = |logFC/t|
# 当 t 不为 0 时计算，否则返回 NA
SE <- ifelse(as.numeric(output_DEGs[, "t"]) != 0,
             abs(as.numeric(output_DEGs[, "logFC"]) / as.numeric(output_DEGs[, "t"])),
             NA)
output_DEGs <- cbind(output_DEGs, SE = SE)

# 如果有额外的注释信息（如 EntrezID、GeneDescription、Chromosome 等），
# 可通过 merge() 将外部注释表与结果合并，示例如下：
# annotation_table <- read.table("gene_annotation.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
# output_DEGs <- merge(output_DEGs, annotation_table, by.x="Gene", by.y="GeneSymbol", all.x=TRUE)

# 可选：调整列顺序，使结果更直观（例如把 Gene, logFC, SE, AveExpr, t, P.Value, adj.P.Val, B 放在前面）
output_DEGs <- as.data.frame(output_DEGs, stringsAsFactors = FALSE)
desired_order <- c("Gene", "logFC", "SE", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
existing_cols <- intersect(desired_order, colnames(output_DEGs))
output_DEGs <- output_DEGs[, c(existing_cols, setdiff(colnames(output_DEGs), existing_cols))]

# 将输出表格写入文件，确保所有信息都展示
write.table(output_DEGs, file = "DE_significant_genes.xls", sep = "\t", quote = FALSE, row.names = FALSE)

step6_end <- Sys.time()
cat("[Step 6] 运行完毕，耗时:", round(difftime(step6_end, step6_start, units = "secs"), 2), "秒\n\n")



# ========================== 7. 绘制热图，在下方显示分组样本数量 ==========================
step7_start <- Sys.time()
cat("[Step 7] 开始绘制热图\n")

# 1. 准备绘图数据 ------------------------------------------------------------
ordered_DEGs <- significant_DEGs[order(as.numeric(as.vector(significant_DEGs$logFC))), ]
ordered_gene_names <- rownames(ordered_DEGs)
total_DEG_count <- length(ordered_gene_names)

if (total_DEG_count > (max_display_genes * 2)) {
  selected_gene_set <- ordered_gene_names[c(1:max_display_genes, 
                                            (total_DEG_count - max_display_genes + 1):total_DEG_count)]
} else {
  selected_gene_set <- ordered_gene_names
}
heatmap_expr <- combined_expr[selected_gene_set, ]

# 2. 构建样本注释信息，并设置分组颜色 ------------------------------------------
# 修改实验组名称为 "Treat"
sample_annotation <- data.frame(Group = factor(c(rep("Control", num_ctrl), 
                                                 rep("Treat", num_treat))))
rownames(sample_annotation) <- colnames(combined_expr)

annotation_colors <- list(
  Group = c("Control" = "#66C2A5",    # 浅绿
            "Treat"   = "#FC8D62")     # 浅橙
)

# 3. 定义更丰富的调色板，使用 RColorBrewer 的 "RdYlBu" 反转调色板 ------------------
color_palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(255)

# 4. 准备热图对象（添加主标题），silent=TRUE 不直接绘制 --------
heatmap_result <- pheatmap(
  mat               = heatmap_expr,
  annotation_col    = sample_annotation,
  annotation_colors = annotation_colors,
  color             = color_palette,
  cluster_cols      = FALSE,
  show_colnames     = FALSE,
  scale             = "row",
  fontsize          = 12,
  fontsize_row      = 7,
  fontsize_col      = 10,
  border_color      = NA,
  main              = "Differential Expression Heatmap",
  silent            = TRUE
)

# 5. 计算对照组和实验组的样本数量，并拼接文本 ------------------------------------
control_count <- num_ctrl
treatment_count <- num_treat
sample_text <- paste("Control:", control_count, "| Treat:", treatment_count)

# 6. 绘制到 PDF，采用 grid 分区，将热图和文本分上下两部分 -----------------------
library(grid)
pdf(file = "DE_heatmap.pdf", height = 8, width = 10)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, 
                                           heights = unit(c(0.93, 0.07), "npc"))))
# 第1行：绘制热图
pushViewport(viewport(layout.pos.row = 1))
grid.draw(heatmap_result$gtable)
popViewport()
# 第2行：绘制文本（分组样本数量）
pushViewport(viewport(layout.pos.row = 2))
grid.text(label = sample_text,
          x = 0.5, 
          y = 0.5,
          gp = gpar(fontsize = 13, fontface = "bold", col = "#333333"))
popViewport()
popViewport()
dev.off()

step7_end <- Sys.time()
cat("[Step 7] 运行完毕，耗时:", 
    round(difftime(step7_end, step7_start, units="secs"), 2), "秒\n\n")

# ========================== 8. 绘制火山图 ==========================
step8_start <- Sys.time()
cat("[Step 8] 开始绘制火山图...\n")

# 为每个基因添加显著性标记，输出带空格的类别名称
significance_status <- ifelse(
  (all_diff_results$adj.P.Val < threshold_adjP & abs(all_diff_results$logFC) > threshold_logFC),
  ifelse(all_diff_results$logFC > threshold_logFC, "Up regulated", "Down regulated"),
  "Not Significant"
)
all_diff_results$Significance <- significance_status

# 构造火山图，并更新颜色映射，使类别名称与 ifelse 返回一致
volcano_plot <- ggplot(all_diff_results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significance), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("Down regulated"  = "#1E90FF", 
                                "Not Significant" = "#808080", 
                                "Up regulated"    = "#FF4500")) +
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC), 
             linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = -log10(threshold_adjP), 
             linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal(base_size = 14) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, color = "#2F4F4F"),
        axis.title    = element_text(face = "bold", color = "#2F4F4F"),
        axis.text     = element_text(color = "#2F4F4F"),
        panel.grid.major = element_line(color = "#D3D3D3", linetype = "dotted"),
        panel.grid.minor = element_blank())

# 保存火山图到 PDF
pdf(file = "DE_volcano.pdf", width = 6, height = 6)
print(volcano_plot)
dev.off()

step8_end <- Sys.time()
cat("[Step 8] 运行完毕，耗时:", round(difftime(step8_end, step8_start, units="secs"), 2), "秒\n\n")

# ========================== 9. PCA 主成分分析 ==========================
step9_start <- Sys.time()
cat("[Step 9] 开始优化 PCA 主成分分析...\n")

# 确保 ggplot2 和 ggrepel 已安装（用于美化标签）
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")  # 使标签不重叠
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr")  # 用于添加置信椭圆
}
library(ggplot2)
library(ggrepel)
library(ggpubr)

# 计算 PCA
pca_result <- prcomp(t(combined_expr), scale. = TRUE)

# 构建 PCA 数据框
pca_df <- data.frame(
  Sample = colnames(combined_expr),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = factor(c(rep("Control", num_ctrl), rep("Treatment", num_treat)))  # 分组信息
)

# 计算 PCA 贡献率（用于坐标轴标签）
pca_var <- pca_result$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 1)

# 绘制 PCA 图
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  # 增加置信椭圆，提高数据分布的可视化
  stat_ellipse(level = 0.95, linetype = "dashed", size = 1) +
  # 绘制点，调整大小和透明度，避免重叠影响
  geom_point(size = 4, alpha = 0.8, shape = 16) +
  # 添加不重叠的文本标签
  geom_text_repel(aes(label = Sample), size = 4, fontface = "bold", max.overlaps = 15) +
  # 自定义颜色，提高视觉冲击力（渐变色）
  scale_color_manual(values = c("Control" = "#0072B2",  # 深蓝
                                "Treatment" = "#E69F00")) + # 橙色
  # 修改坐标轴标签，使其更专业
  labs(title = "PCA Analysis of Gene Expression",
       x = paste("PC1 (", pca_var_perc[1], "%)", sep = ""),
       y = paste("PC2 (", pca_var_perc[2], "%)", sep = ""),
       color = "") +  # 移除图例标题
  # 美化主题
  theme_classic(base_size = 16) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 18, color = "#333333"),
    axis.title    = element_text(face = "bold", size = 16, color = "#333333"),
    axis.text     = element_text(size = 14, color = "#555555"),
    legend.text   = element_text(size = 14, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "#DDDDDD", linetype = "dotted"),
    panel.grid.minor = element_blank()
  )

# 保存 PCA 图
ggsave("DE_PCA.pdf", pca_plot, width = 7, height = 7, dpi = 300)

step9_end <- Sys.time()
cat("[Step 9] PCA 分析优化完成，耗时:", 
    round(difftime(step9_end, step9_start, units="secs"), 2), "秒\n\n")

# 加载必要的库
library(ggplot2)
library(RColorBrewer)

# 假设我们已经有了一个数据框，包含基因的log2FoldChange（logFC）和p值（adj.P.Val）
# 生成一个示例数据框
set.seed(123)
gene_data <- data.frame(
  Gene = paste("Gene", 1:10000),
  logFC = c(rnorm(5000, 2, 1), rnorm(5000, -2, 1)),  # log2FoldChange
  adj.P.Val = runif(10000, 0, 1)  # 调整后的p值
)

# 根据logFC和p值排序
gene_data <- gene_data[order(gene_data$logFC, decreasing = TRUE),]

# 创建渐变色
color_palette <- scale_color_gradientn(
  colors = rev(brewer.pal(9, "RdYlBu")), 
  values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1))
)

# 绘制基因排序点图
plot <- ggplot(gene_data, aes(x = seq_along(logFC), y = logFC)) +
  geom_point(aes(color = adj.P.Val, size = abs(logFC)), alpha = 0.8) +  # 使用p值作为颜色，logFC的绝对值作为大小
  color_palette +  # 使用渐变色来表示p值
  scale_size_continuous(range = c(5, 12)) +  # 增大点的大小范围
  labs(
    title = "Gene Ranking Dotplot", 
    x = "Ranking of Differentially Expressed Genes", 
    y = "Log2 Fold Change",
    color = "Adjusted P-value",
    size = "Log2 Fold Change"
  ) +
  theme_minimal(base_size = 16) +  # 使用最简洁的背景，调整字体大小
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "#2F4F4F", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold", color = "#2F4F4F"),
    axis.text = element_text(size = 14, color = "#555555"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "#D3D3D3", linetype = "dotted"),
    panel.grid.minor = element_blank()
  )

# 保存为PDF，调整尺寸和分辨率
ggsave("123.pdf", plot, width = 10, height = 8, dpi = 300)


# 假定你的差异分析结果为 all_diff_results，行名为基因名
de_table <- all_diff_results

# 先筛选显著差异基因（可按你的阈值，也可以全表top）
sig_degs <- de_table[with(de_table, abs(logFC) > threshold_logFC & adj.P.Val < threshold_adjP), ]

# 按logFC排序，上调Top25&下调Top25
up_top <- rownames(sig_degs)[order(-sig_degs$logFC)][1:25]     # logFC最大 top25
down_top <- rownames(sig_degs)[order(sig_degs$logFC)][1:25]    # logFC最小 top25
label_genes <- unique(c(up_top, down_top))

# 添加gene列到全结果表
de_table$Gene <- rownames(de_table)

# 增加一个新列：为top50标记gene name，否则空
de_table$label <- ifelse(de_table$Gene %in% label_genes, de_table$Gene, "")

# 颜色分组（可选，可按你的Significance或重新定义）
de_table$Group <- "Not Significant"
de_table$Group[de_table$logFC > threshold_logFC & de_table$adj.P.Val < threshold_adjP] <- "Up regulated"
de_table$Group[de_table$logFC < -threshold_logFC & de_table$adj.P.Val < threshold_adjP] <- "Down regulated"

# 绘制火山图
volcano_plot <- ggplot(de_table, aes(x=logFC, y=-log10(adj.P.Val), color=Group)) +
  geom_point(alpha=0.7, size=2) +
  geom_text_repel(aes(label=label), # 只标有label的点
                  size=3,
                  max.overlaps=50,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.color="grey50",
                  show.legend = FALSE) +
  scale_color_manual(values = c("Up regulated"="#FF4500", 
                                "Down regulated"="#1E90FF", 
                                "Not Significant"="#808080")) +
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC), linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(threshold_adjP), linetype="dashed", color="black") +
  labs(title="Volcano Plot with Top  Genes Labeled",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value") +
  theme_minimal(base_size = 14) +
  theme(plot.title=element_text(face="bold", hjust=0.5, color="#2F4F4F"))

# 保存
ggsave("DE_volcano_with_labels.pdf", volcano_plot, width=8, height=8, dpi=300)
# ========================== 全部步骤完成 ==========================

# 假设你的显著性标记列为 Group, 可根据你的实际列名适配
group_counts <- table(de_table$Group)
up_count   <- ifelse("Up regulated"   %in% names(group_counts), group_counts["Up regulated"], 0)
down_count <- ifelse("Down regulated" %in% names(group_counts), group_counts["Down regulated"], 0)
not_count  <- ifelse("Not Significant"%in% names(group_counts), group_counts["Not Significant"], 0)
# 自定义图例名字，把数量拼进去
legend_labels <- c(
  paste0("Up regulated (", up_count, ")"),
  paste0("Down regulated (", down_count, ")"),
  paste0("Not Significant (", not_count, ")")
)
# labels的顺序要和scale_color_manual的values顺序一致
volcano_plot <- ggplot(de_table, aes(x=logFC, y=-log10(adj.P.Val), color=Group)) +
  geom_point(alpha=0.7, size=2) +
  geom_text_repel(aes(label=label),
                  size=3,
                  max.overlaps=50,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.color="grey50",
                  show.legend = FALSE) +
  scale_color_manual(
    values = c("Up regulated" = "#FF4500",
               "Down regulated" = "#1E90FF",
               "Not Significant" = "#808080"),
    # 顺序要和图中实际data分组匹配
    breaks = c("Up regulated", "Down regulated", "Not Significant"),
    labels = legend_labels
  ) +
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC), linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(threshold_adjP), linetype="dashed", color="black") +
  labs(title="Volcano Plot with Top  Genes Labeled",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value") +
  theme_minimal(base_size = 14) +
  theme(plot.title=element_text(face="bold", hjust=0.5, color="#2F4F4F"))
# 差异基因名，以data.frame形式输出TXT，表头为gene
diff_gene_list <- data.frame(gene = rownames(significant_DEGs))
write.table(diff_gene_list,
            file = "DEG_geneList.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
# ========================== 7. 绘制带分布密度图的热图 ==========================
step7_start <- Sys.time()
cat("[Step 7] 开始绘制带分布密度图的热图\n")

library(ComplexHeatmap)
library(circlize)

# 准备数据
ordered_DEGs <- significant_DEGs[order(as.numeric(as.vector(significant_DEGs$logFC))), ]
ordered_gene_names <- rownames(ordered_DEGs)
total_DEG_count <- length(ordered_gene_names)

if (total_DEG_count > (max_display_genes * 2)) {
  selected_gene_set <- ordered_gene_names[c(1:max_display_genes, 
                                            (total_DEG_count - max_display_genes + 1):total_DEG_count)]
} else {
  selected_gene_set <- ordered_gene_names
}
heatmap_expr <- combined_expr[selected_gene_set, ]

# 行标准化
heatmap_expr_scaled <- t(scale(t(heatmap_expr)))

# 样本注释
sample_annotation <- data.frame(Group = factor(c(rep("Control", num_ctrl), 
                                                 rep("Treat", num_treat))))
rownames(sample_annotation) <- colnames(combined_expr)

# 分组颜色
annotation_colors <- list(
  Group = c("Control" = "#66C2A5", "Treat" = "#FC8D62")
)

# 创建分布密度图注释
ha_top <- HeatmapAnnotation(
  Distribution = anno_density(heatmap_expr_scaled, type = "heatmap", 
                              which = "column", height = unit(2, "cm")),
  Group = sample_annotation$Group,
  col = annotation_colors,
  annotation_name_side = "left"
)

# 颜色方案
color_palette <- colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026"))

# 创建热图
ht <- Heatmap(
  heatmap_expr_scaled,
  name = "Expression",
  col = color_palette,
  top_annotation = ha_top,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "Differential Expression Heatmap with Distribution",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

# 保存到PDF
pdf("DE_heatmap_with_distribution.pdf", width = 12, height = 10)
draw(ht)
dev.off()

step7_end <- Sys.time()
cat("[Step 7] 运行完毕，耗时:", 
    round(difftime(step7_end, step7_start, units="secs"), 2), "秒\n\n")


cat("所有步骤已成功完成！\n")
