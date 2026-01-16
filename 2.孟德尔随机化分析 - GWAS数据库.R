library(VariantAnnotation)   # 用于处理变异注释数据
library(gwasglue)            # 用于 GWAS 数据转换
library(TwoSampleMR)         # 用于 Two-Sample Mendelian Randomization 分析
library(R.utils)             # 提供工具函数
library(qqman)               # 用于绘制 Manhattan 图和 QQ 图
library(circular)            # 用于绘制圆形图
library(ggplot2)             # 用于高级数据可视化
library(RadialMR)            # 加载用于Radial MR分析的包


# 设置工作目录
workDir <- "D:/2.椎间盘衰老/Q1区纯生信7.8分/15.有意义的基因单独MR分析出图表"  # 工作目录
# 设置工作目录
setwd(workDir)  # 切换到指定工作目录
outcomeFiles <- list.files(pattern = "\\.vcf.gz$")          # 自动读取此目录下的.gz文件

# 定义输入文件、输出文件夹和工作目录
exposureFile <- "eqtl_filtered.txt"                   # 暴露数据文件路径
outcomeFile <- outcomeFiles[1]         # 结果数据文件路径
diseaseName <- tools::file_path_sans_ext(outcomeFile)   # 从结果文件名中提取疾病名称（去除后缀）
resultDir <- "results_output1"                           # 结果输出文件夹名称



# 检查输入文件是否存在
if (!file.exists(exposureFile)) {
  stop("Error: Exposure data file not found. Please check the path!")
}
if (!file.exists(outcomeFile)) {
  stop("Error: Outcome data file not found. Please check the path!")
}

# 读取暴露数据
expData <- read_exposure_data(
  filename = exposureFile,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  pval_col = "pval",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  phenotype_col = "exposure",  # 基因名称在此列
  id_col = "id",
  samplesize_col = "samplesize",
  chr_col = "chr",
  pos_col = "pos",
  clump = FALSE
)

#读取结局数据的vcf文件,并对数据进行格式转换
vcf3=readVcf(outcomeFiles)
outcomeData=gwasvcf_to_TwoSampleMR(vcf=vcf3, type="outcome")

#从结局数据中提取工具变量
outcomedata2=merge(expData, outcomeData, by.x="SNP", by.y="SNP")
#  outcome.csv
write.csv(outcomedata2[ , -(2:ncol(expData)) ],
          file = "outcome_instruments.csv",
          row.names = FALSE)

# 然后读的时候，一定要把文件名用引号包起来：
outData <- read_outcome_data(
  snps = expData$SNP,
  filename = "outcome_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "pval.outcome",
  eaf_col = "eaf.outcome"
)


# 提取所有唯一的暴露 ID（基因名字）
uniqueExp <- unique(expData$exposure)
numExp <- length(uniqueExp)
progBar <- txtProgressBar(min = 0, max = numExp, style = 3)
stepCounter <- 1

# 定义一个函数来替换特殊字符为空格
clean_filename <- function(filename) {
  filename <- gsub("[^[:alnum:] ]", " ", filename)  # 替换非字母数字字符为空格
  filename <- gsub("\\s+", " ", filename)  # 替换多个空格为一个空格
  filename <- trimws(filename)  # 去除两端的空格
  return(filename)
}

# 创建一个空的数据框，用于存储所有odds ratios结果
all_odds_ratios <- data.frame()


# 添加代码以跳过SNP不足的暴露因素
for (i in seq_along(uniqueExp)) {
  currentID <- uniqueExp[i]  # 使用基因名称作为暴露 ID
  currentID_clean <- clean_filename(currentID)  # 清理基因名称中的特殊字符
  
  cat("Step", stepCounter, ": Processing exposure", currentID, "(Progress:", i, "/", numExp, ")\n")
  stepCounter <- stepCounter + 1
  
  # 创建针对当前暴露的文件夹
  exposureDir <- file.path(resultDir, currentID_clean)
  if (!dir.exists(exposureDir)) {
    dir.create(exposureDir, recursive = TRUE)
    cat("Exposure folder", exposureDir, "has been created.\n")
  } else {
    cat("Exposure folder", exposureDir, "already exists; some files may be overwritten.\n")
  }
  
  # 筛选当前暴露的数据子集
  currentSubset <- expData[expData$exposure == currentID, ]
  
  # 若数据为空，则跳过
  if (nrow(currentSubset) == 0) {
    warning(paste("Warning: Exposure", currentID, "has no data. Skipping!"))
    setTxtProgressBar(progBar, i)
    next
  }
  
  # 为结果数据添加疾病名称字段，便于后续合并
  outData$outcome <- diseaseName
  
  # 合并暴露与结果数据（对齐等位基因方向）
  mergedData <- harmonise_data(currentSubset, outData)
  
  # 若合并后数据为空，则跳过
  if (nrow(mergedData) == 0) {
    warning(paste("Warning: Merged data for exposure", currentID, "is empty. Skipping!"))
    setTxtProgressBar(progBar, i)
    next
  }
  
  # 检查是否有足够的SNP进行MR分析
  if (nrow(mergedData) < 2) {  # 如果合并后的数据行数小于2，说明SNP太少
    warning(paste("Warning: Not enough SNPs for MR analysis of", currentID, ". Skipping!"))
    setTxtProgressBar(progBar, i)
    next
  }
  
  # 根据结果数据的 p 值过滤（保留 p > 5e-06 的记录）  
  filteredData <- mergedData[mergedData$pval.outcome > 5e-06, ]
  if (nrow(filteredData) < 1) {
    warning(paste("Warning: No data left after filtering for exposure", currentID, ". Skipping!"))
    setTxtProgressBar(progBar, i)
    next
  }
  
  # 保存合并且过滤后的数据到 CSV 文件，文件名使用暴露基因名字
  write.csv(filteredData, file = file.path(exposureDir, paste0(currentID_clean, "_harmonised_data.csv")), row.names = FALSE)
  
  # 筛选满足 MR 要求的有效工具变量
  validIV <- filteredData[filteredData$mr_keep == "TRUE", ]
  write.csv(validIV, file = file.path(exposureDir, paste0(currentID_clean, "_valid_SNPs.csv")), row.names = FALSE)
  
  # 运行 MR 分析，使用 tryCatch 捕获可能错误（如非数值参数错误）
  mrResult <- tryCatch({
    mr(filteredData)
  }, error = function(e) {
    warning(paste("MR 分析出错，暴露", currentID, "跳过：", e$message))
    return(NULL)
  })
  
  # 若 MR 分析返回 NULL，则跳过当前暴露
  if (is.null(mrResult)) {
    setTxtProgressBar(progBar, i)
    next
  }
  
  # 生成 OR 结果并保存，文件名使用暴露基因名字
  orResult <- generate_odds_ratios(mrResult)
  write.csv(orResult, file = file.path(exposureDir, paste0(currentID_clean, "_MR_results.csv")), row.names = FALSE)
  
  # 保存 OR 结果到所有汇总数据框中
  all_odds_ratios <- rbind(all_odds_ratios, orResult)
  
  # 进行异质性检验并保存结果，文件名使用暴露基因名字
  heteroResult <- tryCatch({
    mr_heterogeneity(filteredData)
  }, error = function(e) {
    warning(paste("异质性检验出错，暴露", currentID, "跳过：", e$message))
    return(NULL)
  })
  
  if (!is.null(heteroResult)) {
    write.csv(heteroResult, file = file.path(exposureDir, paste0(currentID_clean, "_heterogeneity.csv")), row.names = FALSE)
  }
  
  # 进行多效性检验并保存结果，文件名使用暴露基因名字
  pleiotropyResult <- tryCatch({
    mr_pleiotropy_test(filteredData)
  }, error = function(e) {
    warning(paste("多效性检验出错，暴露", currentID, "跳过：", e$message))
    return(NULL)
  })
  
  if (!is.null(pleiotropyResult)) {
    write.csv(pleiotropyResult, file = file.path(exposureDir, paste0(currentID_clean, "_pleiotropy.csv")), row.names = FALSE)
  }
  
  # 进行单 SNP 分析并保存结果，文件名使用暴露基因名字
  singleSNPResult <- tryCatch({
    mr_singlesnp(filteredData)
  }, error = function(e) {
    warning(paste("单 SNP 分析出错，暴露", currentID, "跳过：", e$message))
    return(NULL)
  })
  
  if (!is.null(singleSNPResult)) {
    write.csv(singleSNPResult, file = file.path(exposureDir, paste0(currentID_clean, "_singlesnp.csv")), row.names = FALSE)
  }
  
  # 进行留一法敏感性分析并保存结果，文件名使用暴露基因名字
  leaveOneOutResult <- tryCatch({
    mr_leaveoneout(filteredData)
  }, error = function(e) {
    warning(paste("留一法分析出错，暴露", currentID, "跳过：", e$message))
    return(NULL)
  })
  
  if (!is.null(leaveOneOutResult)) {
    write.csv(leaveOneOutResult, file = file.path(exposureDir, paste0(currentID_clean, "_leaveoneout.csv")), row.names = FALSE)
  }
  
  # 计算每个工具变量的 F 统计量（F = (beta/se)^2）并保存，文件名使用暴露基因名字
  currentSubset$F_stat <- tryCatch({
    (currentSubset$beta.exposure / currentSubset$se.exposure)^2
  }, error = function(e) {
    warning(paste("计算 F 统计量出错，暴露", currentID, "跳过：", e$message))
    return(NULL)
  })
  
  if (!is.null(currentSubset$F_stat)) {
    write.csv(currentSubset, file = file.path(exposureDir, paste0(currentID_clean, "_instrument_strength.csv")), row.names = FALSE)
  }
  
  # ------------------------------
  # 绘制散点图
  tryCatch({
    pdf(file = file.path(exposureDir, paste0(currentID_clean, "_scatter_plot.pdf")), width = 7.5, height = 7)
    scatterPlot <- mr_scatter_plot(mrResult, filteredData)
    print(scatterPlot)
    dev.off()
  }, error = function(e) {
    warning(paste("绘制散点图时出错，暴露", currentID, "跳过：", e$message))
  })
  
  # ------------------------------
  # 绘制森林图
  tryCatch({
    pdf(file = file.path(exposureDir, paste0(currentID_clean, "_forest_plot.pdf")), width = 7, height = 5.5)
    forestPlot <- mr_forest_plot(singleSNPResult)
    print(forestPlot)
    dev.off()
  }, error = function(e) {
    warning(paste("绘制森林图时出错，暴露", currentID, "跳过：", e$message))
  })
  
  # ------------------------------
  # 绘制漏斗图
  tryCatch({
    pdf(file = file.path(exposureDir, paste0(currentID_clean, "_funnel_plot.pdf")), width = 7, height = 6.5)
    funnelPlot <- mr_funnel_plot(singlesnp_results = singleSNPResult)
    print(funnelPlot)
    dev.off()
  }, error = function(e) {
    warning(paste("绘制漏斗图时出错，暴露", currentID, "跳过：", e$message))
  })
  
  # ------------------------------
  # 绘制留一法敏感性图
  tryCatch({
    pdf(file = file.path(exposureDir, paste0(currentID_clean, "_leaveoneout_plot.pdf")), width = 7, height = 5.5)
    leaveOneOutPlot <- mr_leaveoneout_plot(leaveOneOutResult)
    print(leaveOneOutPlot)
    dev.off()
  }, error = function(e) {
    warning(paste("绘制留一法图时出错，暴露", currentID, "跳过：", e$message))
  })
  
  # 构造曼哈顿图数据
  manhattanData <- data.frame(
    SNP = currentSubset$SNP,
    CHR = currentSubset$chr.exposure,
    BP  = currentSubset$pos.exposure,
    P   = currentSubset$pval.exposure
  )
  # 检查是否成功构建数据
  if (nrow(manhattanData) > 0) {
    # 打开PDF设备
    pdf(file = file.path(exposureDir, paste0(currentID_clean, "_manhattan_plot.pdf")), width = 7, height = 5)
    
    # 绘制曼哈顿图
    tryCatch({
      # 使用 plot() 绘制曼哈顿图
      # 对于曼哈顿图，x 轴是 SNP 位置，y 轴是 -log10(p-value)
      plot(manhattanData$BP, -log10(manhattanData$P),
           xlab = "SNP Position", ylab = "-log10(p-value)",
           main = paste("Manhattan Plot -", currentID_clean),
           pch = 19, col = ifelse(manhattanData$P < 1e-5, "red", "blue"), # 高显著性的 SNP 用红色标记
           cex = 0.5, las = 1)  # cex 为点的大小，las 为轴标签的方向
      
      # 绘制基因组的染色体分隔线（根据您的数据）
      abline(h = -log10(5e-8), col = "green", lty = 2)  # 设置显著性水平为 5e-8
      
      dev.off()  # 关闭图形设备并保存文件
    }, error = function(e) {
      warning(paste("绘制曼哈顿图时出错，暴露", currentID, "跳过：", e$message))
      dev.off()  # 确保即使发生错误也关闭图形设备
    })
  } else {
    warning(paste("曼哈顿图数据为空，无法绘制图形：", currentID))
  }
  
  
  # ------------------------------
  # 绘制 QQ 图
  tryCatch({
    pdf(file = file.path(exposureDir, paste0(currentID_clean, "_qq_plot.pdf")), width = 7, height = 7)
    qq(p = currentSubset$pval.exposure, main = paste("QQ Plot - p-value distribution for", currentID_clean))
    dev.off()
  }, error = function(e) {
    warning(paste("绘制 QQ 图时出错，暴露", currentID, "跳过：", e$message))
  })
  
  # ------------------------------
  # 绘制 SNP 位置与效应大小的散点图
  tryCatch({
    pdf(file = file.path(exposureDir, paste0(currentID_clean, "_effect_size_plot.pdf")), width = 7, height = 7)
    plot(currentSubset$pos.exposure, currentSubset$beta.exposure,
         main = paste("SNP Position vs. Effect Size -", currentID_clean),
         xlab = "SNP Position",
         ylab = "Effect Size",
         pch = 19,
         col = "blue")
    dev.off()
  }, error = function(e) {
    warning(paste("绘制 SNP 位置与效应大小的图时出错，暴露", currentID, "跳过：", e$message))
  })
  
  # ------------------------------
  # Radial MR 分析及图形生成
  tryCatch({
    radial_data <- format_radial(
      mergedData$beta.exposure, 
      mergedData$beta.outcome,    
      mergedData$se.exposure,     
      mergedData$se.outcome,      
      mergedData$SNP              
    )
    radial_data <- na.omit(radial_data)
    radial_results <- ivw_radial(radial_data, alpha = 0.05, weights = 3)
    
    # 生成 Radial 图并保存
    radial_prefix <- paste0(currentID_clean, "_radial_analysis")
    pdf(file = file.path(exposureDir, paste0(radial_prefix, "_radial_plot.pdf")), width = 7, height = 6.5)
    radial_plot_obj <- plot_radial(radial_results)
    if (!is.null(radial_plot_obj)) {
      print(radial_plot_obj)
    } else {
      message("Radial图生成失败或返回NULL。")
    }
    dev.off()
    
    # 保存Radial分析中识别的异常值（outliers）
    radial_outliers <- radial_results$outliers  # 获取异常值
    if (!is.null(radial_outliers) && length(radial_outliers) > 0 && !any(is.na(radial_outliers))) {
      write.csv(radial_outliers, file = file.path(exposureDir, paste0(currentID_clean, "_radial_outliers.csv")), row.names = FALSE)  # 保存异常值
    } else {
      message("Radial分析中未检测到异常值或异常值数据为空。")
    }
    
    # 保存Radial分析中的SNP变异数据信息
    variant_info <- radial_results$data  # 获取SNP变异数据
    if (!is.null(variant_info) && nrow(variant_info) > 0) {
      write.csv(variant_info, file = file.path(exposureDir, paste0(currentID_clean, "_radial_variant_data.csv")), row.names = FALSE)  # 保存变异数据
    } else {
      message("Radial分析中未检测到变异数据或数据为空。")
    }
  }, error = function(e) {
    warning(paste("Radial MR分析时出错，暴露", currentID, "跳过：", e$message))
  })
  
  # 更新进度条
  setTxtProgressBar(progBar, i)
}


# 关闭进度条并输出完成信息
close(progBar)

# 将所有的 OR 结果保存为汇总表格
write.csv(all_odds_ratios, file = "odds_ratios.csv", row.names = FALSE)

cat("Processing and analysis for all exposures are complete! Results are saved in folder:", resultDir, "\n")
