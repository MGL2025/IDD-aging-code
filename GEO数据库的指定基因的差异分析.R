library(limma)
library(ggpubr)
library(pROC)
library(ggplot2)
library(gridExtra)
library(rms)
library(regplot)
library(rmda)
library(tools)

setwd("D:/2.椎间盘衰老(1)/Q1区纯生信7.8分/17.诊断模型基因的ROC曲线（验证集gse34095合176205)")

# 读取样本分组信息
targets_file = "GSE34095+176205_targets.csv"
targets_data = read.csv(targets_file, header = TRUE, sep = ",", check.names = FALSE)

# 创建样本ID到分组的映射
sample_to_group = setNames(targets_data$group, targets_data$sample_id)

expFile = "geneexp.csv"
rt = read.csv(expFile, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)

# 从表达矩阵的列名中提取样本ID（假设列名格式为"GSMXXXXX_xxx"）
sample_ids = gsub("(.*)\\_.*", "\\1", colnames(rt))

# 使用映射获取分组信息
Type = as.character(sample_to_group[sample_ids])
# 确保分组格式正确
Type = ifelse(Type == "CON", "CON", "IDD")
y = ifelse(Type == "CON", 0, 1)

conNum = sum(y == 0)
treatNum = sum(y == 1)

# 使用所有基因，不读取基因列表文件
selectedGenes = rownames(rt)
rt_filtered = rt

colors <- c("IDD" = "#FF6347", "CON" = "#4682B4")

result_df <- data.frame(
  Gene = character(),
  P_Value = numeric(),
  AUC = numeric(),
  AUC_Lower_CI = numeric(),
  AUC_Upper_CI = numeric(),
  SE = numeric(),
  Mean_CON = numeric(),
  Mean_IDD = numeric(),
  stringsAsFactors = FALSE
)

cat("开始基因差异分析...\n")
total_genes = nrow(rt_filtered)
pb <- txtProgressBar(min = 0, max = total_genes, style = 3)

for (i in rownames(rt_filtered)) {
  rt1 = data.frame(expression = as.numeric(rt_filtered[i, ]), Type = Type)
  
  t_test_result = t.test(as.numeric(rt_filtered[i, ]) ~ Type)
  p_value = t_test_result$p.value
  SE = t_test_result$stderr
  
  # 修改为使用CON和IDD
  mean_CON = mean(rt1[rt1$Type == "CON", "expression"])
  mean_IDD = mean(rt1[rt1$Type == "IDD", "expression"])
  
  roc1 = roc(y, as.numeric(rt_filtered[i, ]))
  ci1 = ci.auc(roc1, method = "bootstrap")
  ciVec = as.numeric(ci1)
  
  result_df = rbind(result_df, data.frame(
    Gene = i,
    P_Value = p_value,
    AUC = ciVec[2],
    AUC_Lower_CI = ciVec[1],
    AUC_Upper_CI = ciVec[3],
    SE = SE,
    Mean_CON = mean_CON,
    Mean_IDD = mean_IDD
  ))
  
  boxplot = ggplot(rt1, aes(x = Type, y = expression, fill = Type)) +
    geom_boxplot(outlier.colour = "red", outlier.size = 3, width = 0.5, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    geom_jitter(color = "black", size = 2, width = 0.2) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      legend.title = element_blank(),
      legend.position = "none"
    ) +
    labs(y = paste(i, "expression")) +
    stat_compare_means(method = "t.test")
  
  pdf(file = paste0("boxplot.", i, ".pdf"), width = 3.4, height = 4.5)
  print(boxplot)
  dev.off()
  
  pdf(file = paste0("ROC.", i, ".pdf"), width = 5, height = 5)
  plot(roc1, print.auc = TRUE, col = "orange", legacy.axes = TRUE, main = i)
  polygon(c(roc1$specificities, rev(roc1$specificities)), 
          c(roc1$sensitivities, rep(0, length(roc1$specificities))),
          col = rgb(1, 0.647, 0, 0.3), border = NA)
  text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "orange")
  dev.off()
  
  density_plot = ggplot(rt1, aes(x = expression, fill = Type)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    labs(title = paste(i, "Density Plot"), x = "Expression", y = "Density")
  
  pdf(file = paste0("density_plot.", i, ".pdf"), width = 6, height = 6)
  print(density_plot)
  dev.off()
  
  setTxtProgressBar(pb, which(rownames(rt_filtered) == i))
}

close(pb)

write.csv(result_df, file = "gene_analysis_results.csv", row.names = FALSE)
cat("差异分析结果已保存\n")

roc_plots <- list()
boxplot_plots <- list()

for (i in rownames(rt_filtered)) {
  rt1 = data.frame(expression = as.numeric(rt_filtered[i, ]), Type = Type)
  
  boxplot = ggplot(rt1, aes(x = Type, y = expression, fill = Type)) +
    geom_boxplot(outlier.colour = "red", outlier.size = 3, width = 0.5, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    geom_jitter(color = "black", size = 2, width = 0.2) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      legend.title = element_blank(),
      legend.position = "none"
    ) +
    labs(y = paste(i, "expression")) +
    stat_compare_means(method = "t.test")
  
  boxplot_plots[[i]] <- boxplot
  
  roc1 = roc(y, as.numeric(rt_filtered[i, ]))
  roc_plot = ggroc(roc1) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_minimal() +
    ggtitle(paste(i, "AUC = ", round(auc(roc1), 3))) +
    theme(legend.position = "none")
  
  roc_plots[[i]] <- roc_plot
}

combined_plots <- list()
current_directory <- getwd()
output_folder <- file.path(current_directory, "output")

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

for (i in 1:length(boxplot_plots)) {
  combined_plots[[i]] <- grid.arrange(
    boxplot_plots[[i]], 
    roc_plots[[i]], 
    ncol = 2, 
    top = paste("Gene:", names(boxplot_plots)[i])
  )
  
  ggsave(
    filename = file.path(output_folder, paste0("combined_", names(boxplot_plots)[i], ".pdf")),
    plot = combined_plots[[i]],
    width = 10,
    height = 6
  )
}

cat("所有基因的合并图已保存！\n")

combined_data <- data.frame(expression = numeric(), Type = character(), Gene = character())

for (i in rownames(rt_filtered)) {
  rt1 = data.frame(expression = as.numeric(rt_filtered[i, ]), Type = Type)
  rt1$Gene = i
  combined_data = rbind(combined_data, rt1)
}

# 修改图例标签为CON和IDD
combined_boxplot = ggplot(combined_data, aes(x = Gene, y = expression, fill = Type)) +
  geom_boxplot(outlier.colour = "red", outlier.size = 3, width = 0.5, alpha = 0.7) +
  scale_fill_manual(values = colors, labels = c("CON", "IDD")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  labs(
    y = "Expression", 
    x = "Gene", 
    title = "Diagnostic Model Differential Analysis"
  ) +
  stat_compare_means(
    method = "t.test", 
    label = "p.signif",
    label.x = 1.5
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("combined_boxplot_all_genes_with_labels_and_title.pdf", plot = combined_boxplot, width = 15, height = 8)
cat("所有基因的合并箱线图已保存！\n")

roc_list <- list()
auc_list <- c()
leglab <- c()
my_cols <- c("red", "deepskyblue", "forestgreen", "orange", 
             "purple", "gray40", "black", "magenta", "gold", "brown")
while(length(my_cols) < length(selectedGenes)) {
  my_cols <- c(my_cols, rainbow(length(selectedGenes) - length(my_cols)))
}

pdf("All_Genes_Combined_ROC.pdf", width = 8, height = 8)
par(mar = c(5, 6, 4, 2)+0.1, cex = 1.3)

plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "All Genes Combined ROC", 
     cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45)
abline(0, 1, lty = 2, col = "gray70", lwd = 2)

for(i in seq_along(selectedGenes)) {
  gene <- selectedGenes[i]
  gene_exp <- as.numeric(rt_filtered[gene, ])
  roc1 <- roc(y, gene_exp, direction="auto")
  auc1 <- as.numeric(auc(roc1))
  if (auc1 < 0.5) {
    auc1 <- 1-auc1
    roc1 <- roc(y, gene_exp, direction=ifelse(roc1$direction==">","<",">"))
  }
  lines(1 - roc1$specificities, roc1$sensitivities, col = my_cols[i], lwd = 2.5)
  leglab <- c(leglab, paste0(gene, "  AUC=", sprintf("%.3f", auc1)))
}
legend("bottomright", legend = leglab, col = my_cols[1:length(selectedGenes)], 
       lwd = 2.5, lty = 1, bty = "n", cex = 0.8)
dev.off()
cat("所有基因合并ROC图已保存 All_Genes_Combined_ROC.pdf\n")

print(result_df)

cat("==== 基因分析pipeline启动 ====\n")
steps <- c(
  "创建输出目录", "读取表达数据", "筛选基因",
  "样本结构化", "特征准备", "建模/列线图", "校准曲线", "决策曲线"
)
pb <- txtProgressBar(min=0, max=length(steps), width=40, style=3)

result_dir <- "Analysis_Result"
if (!dir.exists(result_dir)) dir.create(result_dir)
setwd(result_dir)
setTxtProgressBar(pb, 1)

input_expr <- file.path("..", "geneexp.csv")
if (!file.exists(input_expr)) stop(sprintf("未找到文件: %s", input_expr))
dat_expr <- tryCatch(
  read.table(input_expr, sep = ",", header = TRUE, row.names = 1, check.names = FALSE), 
  error = function(e) stop("表达文件读取失败")
)
setTxtProgressBar(pb, 2)

feature_genes <- rownames(dat_expr)
dat_expr_filt <- dat_expr
setTxtProgressBar(pb, 3)

df_expr <- as.data.frame(t(dat_expr_filt))
sample_names <- rownames(df_expr)
# 使用从GSE_targets.csv中获取的分组信息
groups <- as.character(sample_to_group[gsub("(.*)\\_.*", "\\1", sample_names)])
df_expr$GroupType <- groups
if (length(unique(groups)) < 2) warning("分组数不足二，分析或有问题")
setTxtProgressBar(pb, 4)

ddinfo <- datadist(df_expr)
options(datadist="ddinfo")
setTxtProgressBar(pb, 5)

model_vars <- setdiff(colnames(df_expr), "GroupType")
reg_formula <- as.formula(paste("GroupType ~", paste(model_vars, collapse=" + ")))
setTxtProgressBar(pb, 6)

lrm_fit <- lrm(reg_formula, data=df_expr, x=TRUE, y=TRUE)
nomo_obj <- nomogram(
  lrm_fit, fun = plogis,
  fun.at = c(0.001,0.1,0.3,0.5,0.7,0.9,0.99), lp=FALSE, funlabel="Disease Risk"
)
pdf("Nomogram_Plot.pdf", width=11, height=6)
plot(nomo_obj)
dev.off()
cat(">> 列线图输出至: Nomogram_Plot.pdf\n")
setTxtProgressBar(pb, 7)

calibrate_obj <- calibrate(lrm_fit, method="boot", B=1000, group=df_expr$GroupType)
pdf("Calibration_Curve.pdf", width=5.5, height=5.5)
plot(calibrate_obj, xlab="Predicted", ylab="Observed", sub=FALSE)
dev.off()
cat(">> 校准曲线输出至: Calibration_Curve.pdf\n")

# 将分组转换为数值格式
df_expr$GroupType <- ifelse(df_expr$GroupType == "CON", 0, 1)
if(!all(df_expr$GroupType %in% c(0, 1))) stop("GroupType不是严格的0/1!")
set.seed(123)
dca_obj <- decision_curve(
  formula = reg_formula,
  data = df_expr,
  thresholds = seq(0, 1, by = 0.01),
  family = binomial(link = "logit"),
  bootstraps = 100
)
if(is.null(dca_obj$derived.data)) stop("dca_obj生成失败，请检查数据！")
pdf("DCA.pdf", width=6, height=6)
plot_decision_curve(
  dca_obj,
  xlab = "Threshold Probability",
  col = "orange",
  confidence.intervals = TRUE,
  standardize = TRUE,
  cost.benefit.axis = TRUE
)
dev.off()

cat(">> DCA决策曲线输出至 DCA.pdf\n")

write.csv(df_expr, file="Filtered_Expression_Matrix.csv", quote=F)
write.csv(as.data.frame(dca_obj$derived.data), file="Decision_Curve_Data.csv")
pred_probs <- predict(lrm_fit, newdata=df_expr, type="fitted")
df_out <- data.frame(Sample=rownames(df_expr), GroupType=df_expr$GroupType, Predicted_Prob=pred_probs)
write.csv(df_out, file="Sample_Predicted_Probabilities.csv")

coef_df <- as.data.frame(summary(lrm_fit))
if("Effect" %in% names(coef_df) & "S.E." %in% names(coef_df)) {
  coef_df$Effect <- as.numeric(as.character(coef_df$Effect))
  coef_df$`S.E.` <- as.numeric(as.character(coef_df$`S.E.`))
  coef_df$OR <- exp(coef_df$Effect)
  coef_df$OR_low <- exp(coef_df$Effect - 1.96 * coef_df$`S.E.`)
  coef_df$OR_high <- exp(coef_df$Effect + 1.96 * coef_df$`S.E.`)
  write.csv(coef_df, file="Model_Coefficients.csv", row.names=FALSE)
}

roc_y <- as.numeric(df_expr$GroupType)
roc_pred <- pred_probs
roc_obj <- roc(roc_y, roc_pred, levels = c(0,1), direction = "<")

my_cols <- c("red", "deepskyblue", "forestgreen", "orange", 
             "purple", "gray40", "black", "magenta", "gold", "brown")
while(length(my_cols) < length(model_vars) + 1) { 
  my_cols <- c(my_cols, rainbow(length(model_vars) + 1 - length(my_cols)))
}

auc_val <- as.numeric(auc(roc_obj))
auc_ci <- ci(roc_obj)

pdf("CombinedGenes_ROC_SCI.pdf", width = 8, height = 8)
par(mar = c(5,6,4,2)+0.1, cex = 1.3)

plot(1 - roc_obj$specificities, roc_obj$sensitivities, type = "l",
     col = my_cols[1], lwd = 4, lty = 1,
     xlab = expression("1 - Specificity"), ylab = "Sensitivity",
     main = "Model ROC Curve", 
     cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45,
     xlim = c(0, 1), ylim = c(0, 1))

abline(0, 1, lty = 2, col = "gray70", lwd = 2)

legend("bottomright", legend = sprintf("AUC = %.3f (95%% CI: %.3f, %.3f)", auc_val, auc_ci[1], auc_ci[2]), 
       col = my_cols[1], lwd = 4, lty = 1, bty = "n", cex = 1.2)

dev.off()

cat(">> 合并主模型 ROC 曲线 (红色) 已输出至 CombinedGenes_ROC_SCI.pdf\n")

# Step 13：增强标注列线图（高亮中位数虚拟患者）
median_point <- sapply(df_expr[,model_vars,drop=FALSE], median)
median_df <- as.data.frame(t(median_point)); median_df$GroupType <- 1
rownames(median_df) <- "AllMedian"
regplot(
  lrm_fit, showP=TRUE, rank="sd", distribution=TRUE,
  observation=median_df, title="Prediction Nomogram"
)

# Step X: 决策曲线分析 (DCA)
cat(">> 开始进行决策曲线分析 (DCA)...\n")

# 构建联合模型的 formula
dca_model_formula <- as.formula(paste("GroupType ~", paste(model_vars, collapse = " + ")))

# 联合模型 DCA
dca_combined <- decision_curve(
  formula = dca_model_formula,
  data = df_expr,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = FALSE,
  bootstraps = 100
)

# 针对每一个单独的基因构建 DCA
dca_list <- list(Combined_Model = dca_combined)  # 初始化结果列表
for (gene in model_vars) {
  formula_gene <- as.formula(paste("GroupType ~", gene))
  dca_gene <- decision_curve(
    formula = formula_gene,
    data = df_expr,
    family = binomial(link = "logit"),
    thresholds = seq(0, 1, by = 0.01),
    confidence.intervals = FALSE,
    bootstraps = 100
  )
  dca_list[[gene]] <- dca_gene
}

# 保存所有 DCA 数据到 CSV（可选：合并所有 derived.data）
all_dca_data <- do.call(rbind, lapply(names(dca_list), function(name) {
  df <- as.data.frame(dca_list[[name]]$derived.data)
  df$Model <- name
  return(df)
}))
write.csv(all_dca_data, file = "DCA_All_Models_Data.csv", row.names = FALSE)

# 绘图
pdf("DCA_Curve_All_Models.pdf", width = 6, height = 6)
plot_decision_curve(
  x = dca_list,
  curve.names = names(dca_list),
  xlab = "Threshold Probability",
  ylab = "Net Benefit",
  col = c("#DE3E54", rainbow(length(model_vars))),  # 设置颜色，第一个是联合模型
  cost.benefit.axis = TRUE,
  confidence.intervals = FALSE,
  standardize = FALSE,
  legend.position = "bottomleft"
)
dev.off()
cat("==== 全流程执行完毕！====\n")