library(VariantAnnotation)   # 用于处理变异注释数据
library(gwasglue)            # 用于 GWAS 数据转换
library(TwoSampleMR)         # 用于 Two-Sample Mendelian Randomization 分析
library(R.utils)             # 提供工具函数
library(ggplot2)             # 用于高级数据可视化
library(RadialMR)            # 加载用于Radial MR分析的包

# 设置工作目录
workDir <- "D:/2.椎间盘衰老(1)/Q1区纯生信7.8分/12.eqtl数据和疾病的MR分析"  # 工作目录
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
vcf3 = readVcf(outcomeFiles)
outcomeData = gwasvcf_to_TwoSampleMR(vcf = vcf3, type = "outcome")

#从结局数据中提取工具变量
outcomedata2 = merge(expData, outcomeData, by.x = "SNP", by.y = "SNP")
#  outcome.csv
write.csv(outcomedata2[ , -(2:ncol(expData))],
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

# 创建空的数据框用于存储结果
all_odds_ratios <- data.frame()
all_heterogeneity_results <- data.frame()
all_pleiotropy_results <- data.frame()

# 添加代码以跳过SNP不足的暴露因素
for (i in seq_along(uniqueExp)) {
  currentID <- uniqueExp[i]  # 使用基因名称作为暴露 ID
  currentID_clean <- clean_filename(currentID)  # 清理基因名称中的特殊字符
  
  cat("Step", stepCounter, ": Processing exposure", currentID, "(Progress:", i, "/", numExp, ")\n")
  stepCounter <- stepCounter + 1
  
  # 筛选当前暴露的数据子集
  currentSubset <- expData[expData$exposure == currentID, ]
  
  # 若数据为空，则跳过
  if (nrow(currentSubset) == 0) {
    warning(paste("Warning: Exposure", currentID, "has no data. Skipping!"))
    setTxtProgressBar(progBar, i)
    next
  }
  
  # 合并暴露与结果数据（对齐等位基因方向）
  outData$outcome <- diseaseName
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
  
  # ------------------------------
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
  
  # 生成 OR 结果并保存到汇总数据框中
  orResult <- generate_odds_ratios(mrResult)
  all_odds_ratios <- rbind(all_odds_ratios, orResult)
  
  # 进行异质性检验并保存结果到汇总数据框中
  heteroResult <- tryCatch({
    mr_heterogeneity(filteredData)
  }, error = function(e) {
    warning(paste("异质性检验出错，暴露", currentID, "跳过：", e$message))
    return(NULL)
  })
  
  if (!is.null(heteroResult)) {
    all_heterogeneity_results <- rbind(all_heterogeneity_results, heteroResult)
  }
  
  # 进行多效性检验并保存结果到汇总数据框中
  pleiotropyResult <- tryCatch({
    mr_pleiotropy_test(filteredData)
  }, error = function(e) {
    warning(paste("多效性检验出错，暴露", currentID, "跳过：", e$message))
    return(NULL)
  })
  
  if (!is.null(pleiotropyResult)) {
    all_pleiotropy_results <- rbind(all_pleiotropy_results, pleiotropyResult)
  }
  
  # 更新进度条
  setTxtProgressBar(progBar, i)
}

# 关闭进度条并输出完成信息
close(progBar)

# 将所有的 OR 结果保存为汇总表格
write.csv(all_odds_ratios, file = "all_odds_ratios.csv", row.names = FALSE)
# 在处理所有暴露后，保存异质性和多效性结果的汇总表格
write.csv(all_heterogeneity_results, file = "all_heterogeneity_results.csv", row.names = FALSE)
write.csv(all_pleiotropy_results, file = "all_pleiotropy_results.csv", row.names = FALSE)

cat("Processing and analysis for all exposures are complete! Results are saved in folder:", resultDir, "\n")
