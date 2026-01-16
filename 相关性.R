# 加载包，用于数据处理、统计分析和绘图
library("limma")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(linkET)

# 设置工作目录（请根据实际情况修改）
setwd("D:/2.椎间盘衰老(1)/Q1区纯生信7.8分/24.训练集单基因和免疫浸润评分的相关性分析")

# =============================== #
# Step 1: 定义文件路径及检查文件存在性  #
# =============================== #

# 定义输入文件的路径
exprFilePath <- "Sample Type Matrix.csv"         # 基因表达数据文件（CSV格式）
targetGeneFile <- "Genes.csv"                   # 待分析基因列表文件
immuneDataPath <- "热图ssgseaScore.csv"             # 免疫评分数据文件

# 判断文件是否存在，若不存在则中止执行
if (!file.exists(exprFilePath)) {
  stop("错误：基因表达数据文件不存在！")
}
if (!file.exists(targetGeneFile)) {
  stop("错误：目标基因列表文件不存在！")
}
if (!file.exists(immuneDataPath)) {
  stop("错误：免疫细胞数据文件不存在！")
}

# =============================== #
# Step 2: 定义数据读取与预处理函数  #
# =============================== #

# 定义函数：读取并处理基因表达矩阵，同时筛选出多个目标基因数据
readExprData <- function(fileExpr, fileGene) {
  cat("步骤 2.1：读取基因表达数据...\n")
  exprRaw <- read.table(fileExpr, header = TRUE, sep = ",", check.names = FALSE)
  
  # ❌ 移除错误的转置操作 - 基因表达数据通常已经是正确格式
  exprMat <- as.matrix(exprRaw)
  rownames(exprMat) <- exprMat[, 1]
  exprVals <- exprMat[, -1, drop = FALSE]
  
  # 设置行和列名，转换为数值矩阵
  rowNamesTmp <- rownames(exprVals)
  colNamesTmp <- colnames(exprVals)
  numExpr <- matrix(as.numeric(as.matrix(exprVals)), 
                    nrow = nrow(exprVals), 
                    dimnames = list(rowNamesTmp, colNamesTmp))
  
  numExpr <- avereps(numExpr)  # 对重复探针做平均处理
  
  cat("步骤 2.2：读取目标基因文件...\n")
  geneList <- read.table(fileGene, header = FALSE, sep = ",", check.names = FALSE)
  targetGene <- as.vector(geneList[, 1])
  
  # 判断目标基因是否存在于表达矩阵中
  if (any(!(targetGene %in% rownames(numExpr)))) {
    missing_genes <- targetGene[!(targetGene %in% rownames(numExpr))]
    stop("错误：以下目标基因未在表达数据中找到：", paste(missing_genes, collapse = ", "))
  }
  
  # 筛选出目标基因数据，并进行转置（样本为行，基因为列）
  exprSelected <- numExpr[targetGene, , drop = FALSE]
  exprSelected <- t(exprSelected)  # 这个转置是正确的，让样本为行
  
  # 返回一个列表，包含处理后的数据和目标基因名称
  return(list(exprData = exprSelected, geneName = targetGene))
}

# 定义函数：读取免疫细胞数据并与表达数据匹配
readImmuneData <- function(fileImm, exprData) {
  cat("步骤 2.3：读取免疫细胞数据...\n")
  immData <- read.table(fileImm, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
  
  # ✅ 正确位置：在这里对GSVA_Scores.csv进行转置
  cat("转置前免疫数据维度:", dim(immData), "\n")
  immData <- t(immData)  # 对GSVA_Scores.csv进行转置
  immData <- as.data.frame(immData)
  cat("转置后免疫数据维度:", dim(immData), "\n")
  
  # 提取在表达数据和免疫数据中共有的样本
  commonSamples <- intersect(rownames(exprData), rownames(immData))
  
  # 判断是否存在共有样本
  if (length(commonSamples) == 0) {
    cat("表达数据样本名前5个:", head(rownames(exprData), 5), "\n")
    cat("免疫数据样本名前5个:", head(rownames(immData), 5), "\n")
    stop("错误：没有共同的样本在表达数据与免疫数据中！")
  }
  
  cat("找到共同样本数:", length(commonSamples), "\n")
  
  exprData <- exprData[commonSamples, , drop = FALSE]
  immData <- immData[commonSamples, , drop = FALSE]
  
  # 去除标准差为0的列（免疫细胞类型），以防止相关性计算错误
  validImmune <- immData[, apply(immData, 2, sd) > 0, drop = FALSE]
  
  # 返回匹配后的表达数据和免疫数据
  return(list(exprData = exprData, immuneData = validImmune))
}

# 定义函数：计算多个目标基因与免疫细胞的Spearman相关性
computeCorrelation <- function(exprData, immuneData, genes) {
  cat("步骤 2.4：计算相关性...\n")
  corrResults <- data.frame()
  
  # 循环遍历每个目标基因
  for (gene in genes) {
    # 循环遍历每个免疫细胞类型
    for (cellType in colnames(immuneData)) {
      # 判断当前免疫细胞数据是否有足够的变异性
      if (sd(immuneData[, cellType]) == 0) {
        next
      }
      # 将免疫细胞数据和目标基因表达数据转为数值向量
      immuneVec <- as.numeric(immuneData[, cellType])
      geneExpr  <- as.numeric(exprData[, gene])
      # 执行Spearman相关性检验
      testResult <- cor.test(immuneVec, geneExpr, method = "spearman")
      
      # 将结果存入数据框中
      tempDF <- data.frame(spec = gene,
                           env  = cellType,
                           r    = as.numeric(testResult$estimate),
                           p    = as.numeric(testResult$p.value))
      corrResults <- rbind(corrResults, tempDF)
    }
  }
  
  # 根据P值判断相关性方向，并存入新列
  corrResults$pd <- ifelse(corrResults$p < 0.05,
                           ifelse(corrResults$r > 0, "positive correlation", "negative correlation"),
                           "not significant")
  
  # 将相关系数取绝对值，方便后续分档
  corrResults$r <- abs(corrResults$r)
  
  # 为相关系数添加分档信息
  corrResults <- corrResults %>%
    mutate(rd = cut(r,
                    breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                    labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6")))
  
  return(corrResults)
}

# =============================== #
# Step 3: 主流程执行数据读取和分析    #
# =============================== #

# 进度条初始化
totalSteps <- 5
pb <- txtProgressBar(min = 0, max = totalSteps, style = 3)

# 3.1 读取表达数据并筛选目标基因
setTxtProgressBar(pb, 1)
cat("主流程 3.1：处理基因表达数据...\n")
exprOut <- readExprData(exprFilePath, targetGeneFile)
exprMat <- exprOut$exprData
targetGene <- exprOut$geneName  # 修正这里：应该是geneName而不是gene
Sys.sleep(0.2)

# 3.2 读取免疫数据并匹配样本
setTxtProgressBar(pb, 2)
cat("主流程 3.2：处理免疫细胞数据...\n")
immuneOut <- readImmuneData(immuneDataPath, exprMat)
exprMat <- immuneOut$exprData
immuneMat <- immuneOut$immuneData
Sys.sleep(0.2)

# 3.3 检查数据维度
setTxtProgressBar(pb, 3)
if (nrow(exprMat) < 5 || nrow(immuneMat) < 5) {
  warning("样本数量较少，结果可能不稳定！")
}
cat("主流程 3.3：样本数检查完毕。\n")
Sys.sleep(0.2)

# 3.4 计算相关性
setTxtProgressBar(pb, 4)
cat("主流程 3.4：计算相关性...\n")
corrDF <- computeCorrelation(exprMat, immuneMat, targetGene)
Sys.sleep(0.2)

# 3.5 保存结果
setTxtProgressBar(pb, 5)
cat("主流程 3.5：写入相关性结果到CSV文件...\n")
write.csv(corrDF, file = "gene_correlation_results.csv", row.names = FALSE)
close(pb)

# 继续执行后面的绘图和热图代码...
# (后面的代码保持不变)


# =============================== #
# Step 4: 绘图（两种风格均保留）        #
# =============================== #

# ---------- 图形风格A：Spectral配色风格 ---------- #
cat("步骤 4A：绘制图形风格A...\n")  # 提示信息
plotA <- qcorrplot(correlate(immuneMat, method = "spearman"), type = "lower", diag = FALSE) +  # 绘制相关性矩阵
  geom_tile() +          # 绘制热图基础图层
  geom_couple(aes(colour = pd, size = rd),                # 绘制耦合线（基因与免疫细胞之间的连线)
              data = corrDF, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))) +  # 设置填充色谱（反转Spectral）
  scale_size_manual(values = c(0.5, 1, 2, 3)) +             # 自定义线宽
  scale_colour_manual(values = c("positive correlation" = "#FF1493",   # 设置正相关为亮粉
                                 "negative correlation" = "#00CED1",   # 负相关为亮青
                                 "not significant"  = "#999999")) +     # 非显著为灰色
  guides(size = guide_legend(title = "abs(Cor)", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "P-value", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Cell-cell cor", order = 3)) +
  labs(x = "immune infiltrating cells", y = "immune infiltrating cells", title = "Gene-immune infiltrating cells Correlation") +  # 设置坐标轴标签和标题
  theme_minimal(base_size = 14) +                       # 使用极简主题
  theme(axis.title.x = element_text(size = 14, face = "bold", color = "black"),  # 设置X轴标题格式
        axis.title.y = element_text(size = 14, face = "bold", color = "black"),  # 设置Y轴标题格式
        axis.text.x  = element_text(size = 12, face = "bold", color = "black", angle = 45, hjust = 1),  # X轴文字旋转
        axis.text.y  = element_text(size = 12, face = "bold", color = "black"),  # Y轴文字格式
        panel.grid.major = element_blank(),             # 去除主网格线
        panel.grid.minor = element_blank(),             # 去除次网格线
        plot.background  = element_rect(fill = "white", color = NA),  # 图形背景设置为白色
        plot.title       = element_text(size = 16, face = "bold", hjust = 0.5))  # 图标题居中加粗

# 保存图形风格A到PDF文件
pdf(file = "cor_newStyle.pdf", width = 10, height = 8)  # 打开PDF设备
print(plotA)            # 打印图形到设备
dev.off()               # 关闭PDF设备

# ---------- 图形风格B：自定义主题 + Pastel2配色 ---------- #
cat("步骤 4B：绘制图形风格B...\n")  # 提示信息

# 定义自定义主题函数（theme_cute），用于图形美化
theme_cute <- function() {
  theme_minimal() +  # 使用极简主题
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),  # 图形背景白色
      panel.background = element_rect(fill = "white", colour = NA),  # 面板背景白色
      panel.grid.major = element_line(colour = "#f0f0f0", size = 0.5), # 主网格线设置
      panel.grid.minor = element_line(colour = "#f0f0f0", size = 0.25),# 次网格线设置
      axis.text        = element_text(size = 9, colour = "black"),   # 坐标轴文字设置
      axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),# X轴文字竖排
      axis.title       = element_text(size = 9, face = "plain", colour = "black"), # 坐标轴标题设置
      plot.title       = element_text(size = 18, face = "plain", colour = "#d35400", hjust = 0.5)  # 图标题设置
    )
}

plotB <- qcorrplot(correlate(immuneMat, method = "spearman"), type = "lower", diag = FALSE) +  # 绘制相关性矩阵
  geom_tile() +                                    # 绘制热图基础图层
  geom_couple(aes(colour = pd, size = rd), data = corrDF, curvature = nice_curvature()) +  # 添加耦合连线
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(9, "Pastel2")),  # 使用Pastel2色系
    name = "Cell-Cell Correlation"                        # 色条名称
  ) +
  scale_size_manual(
    name   = "abs(Cor)",
    values = c("< 0.2" = 0.5, "0.2 - 0.4" = 1, "0.4 - 0.6" = 2, ">= 0.6" = 3),
    labels = c("< 0.2" = "< 0.2", "0.2 - 0.4" = "0.2 - 0.4", "0.4 - 0.6" = "0.4 - 0.6", ">= 0.6" = "≥ 0.6")
  ) +
  scale_colour_manual(
    name   = "p-value",
    values = c("positive correlation" = "#F28C8C",   # 正相关色：亮粉
               "negative correlation" = "#8AB8FF",   # 负相关色：亮蓝
               "not significant"      = "#B2BABB"),  # 非显著色：灰蓝
    labels = c("positive correlation" = "Positive",
               "negative correlation" = "Negative",
               "not significant"      = "Not significant")
  ) +
  guides(
    size   = guide_legend(title = "abs(Cor)", override.aes = list(colour = "black"), order = 2),
    colour = guide_legend(title = "p-value", override.aes = list(size = 2), order = 1),
    fill   = guide_colorbar(title = "Cell-Cell Correlation", order = 3)
  ) +  # 设置图形标题
  theme_cute() +             # 应用自定义主题
  labs(x = NULL, y = NULL)   # 删除X/Y轴标题

# 保存图形风格B到PDF文件
pdf(file = "cor2.pdf", width = 8, height = 6)  # 打开PDF设备
print(plotB)           # 打印图形到设备
dev.off()              # 关闭PDF设备
## ============= 热图展示：目标基因与免疫细胞相关性（P值版本） =============
library(pheatmap)

expr <- exprMat
immune <- immuneMat

# 创建相关性系数矩阵和P值矩阵
cor_mat <- matrix(NA, nrow=ncol(expr), ncol=ncol(immune))
p_mat <- matrix(NA, nrow=ncol(expr), ncol=ncol(immune))
rownames(cor_mat) <- colnames(expr)
colnames(cor_mat) <- colnames(immune)
rownames(p_mat) <- colnames(expr)
colnames(p_mat) <- colnames(immune)

# 计算相关性和P值
for (i in seq_len(ncol(expr))) {
  for (j in seq_len(ncol(immune))) {
    cor_test_result <- cor.test(expr[, i], immune[, j], method = "spearman", use = "pairwise.complete.obs")
    cor_mat[i, j] <- cor_test_result$estimate
    p_mat[i, j] <- cor_test_result$p.value
  }
}

# 创建显著性星号标签矩阵
create_significance_labels <- function(p_values) {
  labels <- matrix("", nrow=nrow(p_values), ncol=ncol(p_values))
  labels[p_values < 0.001] <- "***"
  labels[p_values >= 0.001 & p_values < 0.01] <- "**"
  labels[p_values >= 0.01 & p_values < 0.05] <- "*"
  labels[p_values >= 0.05] <- "ns"  # 不显著用ns表示，也可以用""表示空白
  return(labels)
}

# 创建P值标签矩阵（显示P值和显著性）
labels_mat <- matrix("", nrow=nrow(p_mat), ncol=ncol(p_mat), dimnames=dimnames(p_mat))

for (i in 1:nrow(p_mat)) {
  for (j in 1:ncol(p_mat)) {
    p_val <- p_mat[i, j]
    if (p_val < 0.001) {
      labels_mat[i, j] <- "***"
    } else if (p_val < 0.01) {
      labels_mat[i, j] <- "**"
    } else if (p_val < 0.05) {
      labels_mat[i, j] <- "*"
    } else {
      labels_mat[i, j] <- ""  # 不显著的不显示任何标记
    }
  }
}

# 计算图形尺寸
gene_num <- nrow(cor_mat)
immune_num <- ncol(cor_mat)
pdf_width <- 3 + 0.4 * immune_num
pdf_height <- 2 + 0.4 * gene_num

# 使用相关性系数作为颜色，P值星号作为标签
my_col_fun <- colorRampPalette(c("#6b8e23", "white", "#ba6262"))(100)

pdf("Gene_Immune_Correlation_Heatmap.pdf", width=10, height=4)
pheatmap(
  cor_mat,  # 仍然用相关系数矩阵作为颜色基础
  color = my_col_fun,
  display_numbers = labels_mat,  # 显示显著性星号
  number_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_number = 15,  # 增大字体以便清楚看到星号
  fontsize_row = 12,
  fontsize_col = 13,
  border_color = NA,
  angle_col = 90,  # x轴标签旋转90度
  main = "Correlation between genes and immune cells\n(* p<0.05, ** p<0.01, *** p<0.001)"
)
dev.off()

# 可选：同时生成一个显示P值数值的版本
labels_p_values <- matrix(
  sprintf("%.3f", p_mat),
  nrow=nrow(p_mat), ncol=ncol(p_mat),
  dimnames=dimnames(p_mat)
)

pdf("Gene_Immune_Pvalue_Heatmap.pdf", width=10, height=4)
pheatmap(
  cor_mat,  # 仍然用相关系数矩阵作为颜色基础
  color = my_col_fun,
  display_numbers = labels_p_values,  # 显示P值数值
  number_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_number = 8,
  fontsize_row = 12,
  fontsize_col = 13,
  border_color = NA,
  angle_col = 90,  # x轴标签旋转90度
  main = "Correlation between genes and immune cells\n(P-values shown)"
)
dev.off()

cat("热图已生成：\n")
cat("1. Gene_Immune_Correlation_Heatmap.pdf - 显示显著性星号\n")
cat("2. Gene_Immune_Pvalue_Heatmap.pdf - 显示P值数值\n")

