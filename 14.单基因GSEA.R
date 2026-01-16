############################################################################===
# 在我们使用R软件之前，我们需要对R软件进行一些设置，这些设置是流程化的。我们只需要适当的修改就可以了。
setwd("D:/2.椎间盘衰老(1)/Q1区纯生信7.8分/22.单基因的GSEA通路分析-KEGG")   #设置工作路径
##只需要表达矩阵文件
load(file = "GSE_expr_targets.RData")
#table(DEG$change)
library(tidyverse)
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(readxl)

# 检查数据格式
# 确保GSE15227_expr是数据框或矩阵，且行名为基因名
if (!exists("GSE_expr")) stop("表达矩阵GSE_expr未找到")

# 读取Hub基因列表（假设基因在"A"列）
hub_genes <- read_excel("genes.xlsx", sheet = "Sheet1") %>%
  pull(genes) %>%        # 提取包含基因名的列（根据实际列名调整）
  na.omit() %>%      # 移除空值
  as.character()     # 转换为字符向量

# 检查Hub基因是否存在于表达矩阵中
missing_genes <- setdiff(hub_genes, rownames(GSE_expr))
if (length(missing_genes) > 0) {
  warning("以下基因在表达矩阵中缺失：", paste(missing_genes, collapse = ", "))
  hub_genes <- intersect(hub_genes, rownames(GSE110359_expr)) # 只保留存在的基因
}
#######基于基因相关性排序
# 循环处理每个Hub基因
for (gene in hub_genes) {
  # 提取目标基因表达值（确保转换为数值）
  target_expr <- as.numeric(GSE_expr[gene, ])
  
  # 计算所有基因与目标基因的相关性（处理缺失值）
  cor_values <- apply(GSE_expr, 1, function(x) {
    cor(as.numeric(x), target_expr, method = "spearman", use = "complete.obs")
  })
  
  # 转换为数据框并排序
  ranked_df <- data.frame(
    gene_symbol = names(cor_values),
    correlation = cor_values
  ) %>%
    arrange(desc(correlation))  # 按相关性从高到低排序
  
  # 保存结果（包含基因名和相关性值）
  write.csv(
    ranked_df,
    file = paste0("GSEA_", gene, "_ranked_genes.csv"),
    row.names = FALSE
  )
}


###基于样本分组（高/低表达）
library(limma)
load(file = "data/GSE110359_expr&targets.rda")
data <- as.matrix(GSE110359_expr)
mode(data) <- "numeric"
sgene <- "PTGS1"    #输入进行单基因GSEA的基因名称

group <- ifelse(data[c(sgene),]>median(data[c(sgene),]), "High", "Low")   
group <- factor(group,levels = c("High","Low"))

#差异分析
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(High-Low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg=topTable(fit2,adjust='fdr',number=nrow(data))

#保存单基因分组的所有基因差异结果
write.csv(deg, file = "PTGS1_gsea.csv")
