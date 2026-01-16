#---------------------#
options(timeout=600000)
# 1. 包加载函数 #
# 分析包（Bioconductor）
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
#------------------------------#
# 2. 初始化参数/环境与进度条    #
#------------------------------#
wdPath <- "D:/2.椎间盘衰老/Q1区纯生信7.8分/07.交集基因的GO富集分析"
if (!dir.exists(wdPath)) dir.create(wdPath, recursive = TRUE)
setwd(wdPath)

updateProgress <- function(pb, value, message) {
  setTxtProgressBar(pb, value)
  cat(message, "\n")
}
totalSteps <- 9
pb <- txtProgressBar(min = 0, max = totalSteps, style = 3)

pvalThreshold <- 0.05
padjThreshold <- 1
colorParameter <- if (padjThreshold > 0.05) "pvalue" else "p.adjust"
updateProgress(pb, 1, "第1步：设置参数和工作环境完成")

#------------------------------#
# 3. 读取基因数据               #
#------------------------------#
geneFile <- "DEG_geneList.txt"
if (!file.exists(geneFile)) stop("错误：基因文件 不存在")
geneData <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
updateProgress(pb, 2, "第2步：读取基因数据完成")

geneSymbols <- unique(as.vector(geneData[, 1]))
if (length(geneSymbols) == 0) stop("错误：未找到有效的基因符号")

#------------------------------#
# 4. 基因符号转Entrez ID       #
#------------------------------#
entrezMapping <- mget(geneSymbols, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezMapping)
validGenes <- entrezIDs[!is.na(entrezIDs) & entrezIDs!="NA"]
if (length(validGenes) == 0) stop("错误：无有效的 Entrez 基因 ID")
updateProgress(pb, 3, "第3步：基因符号转换为 EntrezID完成")

#------------------------------#
# 5. GO富集与过滤              #
#------------------------------#
goAnalysis <- enrichGO(
  gene = validGenes,
  OrgDb = org.Hs.eg.db,
  pvalueCutoff = 1, qvalueCutoff = 1,
  ont = "all",
  readable = TRUE
)
goResult <- as.data.frame(goAnalysis)
if (nrow(goResult)==0) stop("警告：未检测到任何富集结果")
filteredGO <- goResult[goResult$pvalue < pvalThreshold & goResult$p.adjust < padjThreshold, ]
updateProgress(pb, 4, "第4步：GO 富集分析及结果过滤完成")

#------------------------------#
# 6. 结果输出                  #
#------------------------------#
outputFile <- "GO_results.txt"
write.table(filteredGO, file = outputFile, sep = "\t", quote = FALSE, row.names = FALSE)
updateProgress(pb, 5, "第5步：富集结果写入文件完成")

#------------------------------#
# 7. 柱状图/气泡图 PDF绘制     #
#------------------------------#
# 柱状图
pdf("GO_barplot.pdf", width=8, height=10)
barPlot <- barplot(
  goAnalysis,
  drop=TRUE, showCategory=10, label_format=50,
  split="ONTOLOGY", color=colorParameter
) +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_fill_gradientn(colors = c("#FF6666", "#FFB266", "#FFFF99", "#99FF99", "#6666FF", "#7F52A0", "#B266FF"))
print(barPlot)
dev.off()

# 气泡图
pdf("GO_bubble.pdf", width=8, height=10)
bubblePlot <- dotplot(
  goAnalysis,
  showCategory = 10,
  orderBy = "GeneRatio",
  label_format = 50,
  split = "ONTOLOGY",
  color = colorParameter
) +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_color_gradientn(colors = c("#FFB266", "#FFFF99", "#99FF99", "#6666FF", "#7F52A0", "#B266FF"))
print(bubblePlot)
dev.off()
updateProgress(pb, 6, "第6步：柱状图和气泡图生成完成")

#------------------------------#
# 8. 分组条形图 ggbarplot       #
#------------------------------#
topGO <- filteredGO %>% group_by(ONTOLOGY) %>% slice_head(n = 10)
pdf("GO_grouped_barplot.pdf", width=11, height=8)
groupBarPlot <- ggbarplot(
  topGO,
  x="Description", y="Count", fill="ONTOLOGY", color="white",
  xlab="", palette="aaas",
  legend="right", sort.val="desc", sort.by.groups=TRUE,
  position=position_dodge(0.9)
) +
  rotate_x_text(75) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size=10, color="black")) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  geom_text(
    aes(label=Count),
    position=position_dodge(0.9), vjust=-0.3, size=3
  )
print(groupBarPlot)
dev.off()
updateProgress(pb, 7, "第7步：分组条形图生成完成")

#------------------------------#
# 9. chord弦图 PDF绘制         #
#------------------------------#
pdf("GO_chord_diagram.pdf", width=12, height=12)
go <- read.delim("GO_results.txt", header=TRUE, stringsAsFactors=FALSE)
top_terms <- go %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup()
insert_linebreak <- function(text, line_length=35) {
  if(nchar(text) <= line_length) return(text)
  paste(strwrap(text, width=line_length), collapse = "\n")
}
top_terms$Description_new <- sapply(top_terms$Description, insert_linebreak)
mat <- table(
  factor(top_terms$ONTOLOGY, levels=c("BP","CC","MF")),
  factor(top_terms$Description_new, levels=unique(top_terms$Description_new))
)
n_ont <- 3
n_term <- ncol(mat)
gap.deg <- c(rep(1, n_ont-1), 10, rep(1, n_term-1), 10)
grid.col <- c(BP="#E69F00", CC="#56B4E9", MF="#009E73",
              setNames(rep("#BBBBBB", n_term), colnames(mat)))
circos.clear()
circos.par(gap.degree = gap.deg, start.degree = 90)
chordDiagram(
  mat,
  grid.col = grid.col,
  transparency = 0.4,
  annotationTrack = c("", "grid"),
  preAllocateTracks = list(track.height = 0.2)
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name,
              facing="clockwise", niceFacing=TRUE,
              adj=c(0,0.5), cex=0.60)
}, bg.border=NA)
title("GO Ontology")
circos.clear()
dev.off()
updateProgress(pb, 8, "第8步：chord弦图生成完成")

# 载入必要的包
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 读取GO富集分析的结果（你的filteredGO）

# TRUE则全部三类（BP, CC, MF），如果只要其中某几类，改为只留需要的类别名称
filteredGO_sub <- filteredGO %>% filter(ONTOLOGY %in% c("BP", "CC", "MF"))

# 分别对每个大类取Top10最显著的GO term
topGO <- filteredGO_sub %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust)

# 高级：保证每一类内，Description是从最显著到最不显著排序，显示时不乱
topGO <- topGO %>% 
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  ungroup()

# GeneRatio转数值（主要用于气泡图横轴）
topGO$GeneRatio_num <- sapply(topGO$GeneRatio, function(x) {
  sp <- unlist(strsplit(as.character(x), "/"))
  as.numeric(sp[1]) / as.numeric(sp[2])
})

# GO富集条形图（分面，和KEGG一样）
bar_colors <- brewer.pal(7, "YlOrRd")
p_bar <- ggplot(topGO, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = bar_colors) +
  labs(x = "GO Term", y = "-log10(Adjusted p-value)", title = "GO Enrichment Analysis (BP/CC/MF)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11)
  )
ggsave("GO_barplot_custom.pdf", p_bar, width = 12, height = 12)

# GO富集气泡图（分面，和KEGG一样）
dot_colors <- brewer.pal(7, "Spectral")
p_dot <- ggplot(topGO, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_color_gradientn(colors = dot_colors) +
  labs(x = "Gene Ratio", y = "GO Term", title = "GO Dotplot (BP/CC/MF)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkred"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11)
  )
ggsave("GO_dotplot_custom.pdf", p_dot, width = 12, height = 12)

# ---------------------------
# GO富集气泡图/条形图 分面 分别显示BP/CC/MF 通路名字加粗
# ---------------------------


# 2. 读取数据（用你自己的GO富集结果路径）

# 3. 只保留BP,CC,MF三类（如想只看部分，改filter内容即可）
filteredGO_sub <- filteredGO %>% filter(ONTOLOGY %in% c("BP", "CC", "MF"))

# 4. 每类top10最显著GO term
topGO <- filteredGO_sub %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust)

# 5. 保证每个类别内部按显著性排序
topGO <- topGO %>% 
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  ungroup()

# 6. GeneRatio数值化
topGO$GeneRatio_num <- sapply(topGO$GeneRatio, function(x) {
  sp <- unlist(strsplit(as.character(x), "/"))
  as.numeric(sp[1]) / as.numeric(sp[2])
})

# 7. 条形图 (KEGG style, 通路名字加粗)
bar_colors <- brewer.pal(7, "YlOrRd")
p_bar <- ggplot(topGO, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = bar_colors) +
  labs(x = "GO Term", y = "-log10(adjusted p-value)", title = "GO Enrichment Analysis (BP/CC/MF)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold") # 通路名字加粗
  )
ggsave("GO_barplot_custom.pdf", p_bar, width = 12, height = 12)

# 8. 气泡图 (KEGG style, 通路名字加粗)
dot_colors <- brewer.pal(7, "Spectral")
p_dot <- ggplot(topGO, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_color_gradientn(colors = dot_colors) +
  labs(x = "Gene Ratio", y = "GO Term", title = "GO Dotplot (BP/CC/MF)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkred"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold") # 通路名字加粗
  )
ggsave("GO_dotplot_custom.pdf", p_dot, width = 12, height = 12)


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# 读取GO富集分析结果
filteredGO <- read.table("GO_results.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

# 筛选三大类，每类取top10最显著的条目
topGO <- filteredGO %>%
  filter(ONTOLOGY %in% c("BP", "CC", "MF")) %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup()

# 按ONTOLOGY分组，然后按Count排序（同组的颜色相同，聚集在一起）
topGO <- topGO %>% 
  arrange(ONTOLOGY, desc(Count))

# 设置Description为因子，控制显示顺序（保持分组效果）
topGO$Description <- factor(topGO$Description, levels = topGO$Description)

# 创建横向点线图，颜色相同的为一组
p <- ggplot(topGO, aes(x = Description, y = Count, color = ONTOLOGY)) +
  # 绘制垂直线段（从0到Count值）
  geom_segment(aes(x = Description, xend = Description, y = 0, yend = Count), 
               size = 1.2, alpha = 0.8) +
  # 在线段顶端添加圆点
  geom_point(size = 4, alpha = 0.9) +
  # 设置三个分类的颜色（颜色相同的为一组）
  scale_color_manual(
    values = c("BP" = "#4472C4",    # 蓝色 - 生物过程
               "CC" = "#70AD47",    # 绿色 - 细胞组分  
               "MF" = "#FF6B35"),   # 橙色 - 分子功能
    labels = c("BP" = "Biological Process", 
               "CC" = "Cellular Component", 
               "MF" = "Molecular Function")
  ) +
  # 移除coord_flip()以实现横向布局
  labs(
    x = "",
    y = "Count",
    color = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(
    # 调整x轴文本角度以适应横向布局
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.5),
    # 调整页边距
    plot.margin = margin(t = 20, r = 20, b = 80, l = 20, unit = "pt")
  ) +
  # 在圆点上方添加数值标签
  geom_text(aes(label = Count), 
            vjust = -0.5, size = 3, color = "black") +
  # 调整y轴范围，给标签留出空间
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# 保存图片
ggsave("GO_dotline_plot_horizontal.pdf", p, width = 16, height = 8, dpi = 300)

# 如果GO term名称太长，可以考虑换行处理的版本
library(stringr)

# 添加换行处理函数
wrap_text <- function(text, width = 30) {
  str_wrap(text, width = width)
}

# 应用文本换行
topGO$Description_wrapped <- sapply(topGO$Description, wrap_text)
topGO$Description_wrapped <- factor(topGO$Description_wrapped, levels = topGO$Description_wrapped)

# 带换行的版本
p_wrapped <- ggplot(topGO, aes(x = Description_wrapped, y = Count, color = ONTOLOGY)) +
  geom_segment(aes(x = Description_wrapped, xend = Description_wrapped, y = 0, yend = Count), 
               size = 1.2, alpha = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(
    values = c("BP" = "#4472C4", "CC" = "#70AD47", "MF" = "#FF6B35"),
    labels = c("BP" = "Biological Process", 
               "CC" = "Cellular Component", 
               "MF" = "Molecular Function")
  ) +
  labs(x = "", y = "Count", color = "") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.5),
    plot.margin = margin(t = 20, r = 20, b = 100, l = 20, unit = "pt")
  ) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("GO_dotline_plot_horizontal_wrapped.pdf", p_wrapped, width = 18, height = 8, dpi = 300)

#------------------------------#
# 10. 关闭进度条               #
#------------------------------#
close(pb)
cat('全部分析及绘图完成！\n')
