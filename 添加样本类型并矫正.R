# ---------------------------- 包管理函数 ----------------------------
# 定义函数：检查并安装指定包（支持CRAN和Bioconductor）
check_install_packages <- function(pkg, bioc = FALSE) {                 
  if (!requireNamespace(pkg, quietly = TRUE)) {                          
    if (bioc) {                                                         
      BiocManager::install(pkg)                                         
    } else {                                                            
      install.packages(pkg)                                             
    }                                                                   
  }                                                                     
  library(pkg, character.only = TRUE)                                   
}

# -------------------------- 包安装与加载 --------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {                 
  install.packages("BiocManager")                                        
}
check_install_packages("limma", bioc = TRUE)        # 差异表达、重复值平均、归一化
check_install_packages("data.table", bioc = FALSE)  # 高效读取文本数据

# -------------------------- 设置工作目录 --------------------------
setwd("I:\\2025教程新n\\14.Q1区7.8分，MR+单细胞\\19.验证集添加样本类型")  # 根据实际情况修改

# -------------------------- 初始化进度条 --------------------------
totalSteps <- 8                                                      # 定义总共8个主要步骤
progressBar <- txtProgressBar(min = 0, max = totalSteps, style = 3)  # 初始化文本进度条
currentStep <- 0                                                     # 初始化当前步骤计数

# -------------------------- 函数定义 --------------------------
# 读取表达矩阵数据（geneMatrix.txt，制表符分隔），并转换为数值矩阵
readExpressionData <- function(filePath) {
  message("步骤1：读取基因表达矩阵数据...")
  dt <- data.table::fread(filePath, header = TRUE, sep = "\t")  # 读取制表符分隔文件
  mat <- as.matrix(dt)                                          # 转换为矩阵
  rownames(mat) <- mat[, 1]                                     # 将第一列设为行名（基因名）
  expr <- mat[, -1, drop = FALSE]                               # 提取表达数据（除第一列）
  numericMat <- matrix(as.numeric(as.matrix(expr)),             # 转换为数值型矩阵
                       nrow = nrow(expr),
                       dimnames = list(rownames(expr), colnames(expr)))
  return(numericMat)
}

# 对表达矩阵进行重复值平均、判断是否需要log2转换以及归一化处理
processExpressionData <- function(exprMat) {
  message("步骤2：处理表达数据：重复值平均、判断log2转换...")
  avgMat <- avereps(exprMat)                                     # 平均重复基因表达值
  qx <- as.numeric(quantile(avgMat, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  needLog <- (qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0)
  if (needLog) {
    avgMat[avgMat < 0] <- 0                                      # 负值设为0，防止log2出错
    avgMat <- log2(avgMat + 1)                                   # 进行log2转换
  }
  message("步骤2：进行数组间归一化...")
  normMat <- normalizeBetweenArrays(avgMat)                      # 数组间归一化
  if (any(is.na(normMat))) {
    normMat[is.na(normMat)] <- 0
  }
  return(normMat)
}

# -------------------------- 主流程 --------------------------
# 步骤1：读取基因表达矩阵数据（geneMatrix.txt）
expMatrix <- readExpressionData("geneMatrix.txt")
currentStep <- currentStep + 1
setTxtProgressBar(progressBar, currentStep)

# 步骤2：预处理表达数据
processedMatrix <- processExpressionData(expMatrix)
currentStep <- currentStep + 1
setTxtProgressBar(progressBar, currentStep)

# 步骤3：读取对照组与处理组样本信息
message("步骤3：读取样本信息文件...")
controlInfo <- tryCatch({
  read.table("control.txt", header = FALSE, sep = "\t", check.names = FALSE)
}, error = function(e) {
  stop("读取 control.txt 文件出错！")
})
treatmentInfo <- tryCatch({
  read.table("treat.txt", header = FALSE, sep = "\t", check.names = FALSE)
}, error = function(e) {
  stop("读取 treat.txt 文件出错！")
})
currentStep <- currentStep + 1
setTxtProgressBar(progressBar, currentStep)

# 步骤4：根据样本信息提取数据
controlData <- processedMatrix[, as.vector(controlInfo[, 1]), drop = FALSE]
treatmentData <- processedMatrix[, as.vector(treatmentInfo[, 1]), drop = FALSE]

# 步骤5：合并样本数据并生成样本类型标签
combinedData <- cbind(controlData, treatmentData)               # 合并对照组与处理组数据
numControl <- ncol(controlData)                                 # 对照组样本数
numTreatment <- ncol(treatmentData)                             # 处理组样本数
sampleTypes <- c(rep("con", numControl), rep("tre", numTreatment))
newColNames <- paste0(colnames(combinedData), "_", sampleTypes)  # 为列添加类型后缀
colnames(combinedData) <- newColNames
currentStep <- currentStep + 1
setTxtProgressBar(progressBar, currentStep)

# 步骤6：生成最终结果矩阵（将基因名作为第一列）
message("步骤6：生成最终样本类型矩阵...")
finalMatrix <- cbind(GeneName = rownames(combinedData), combinedData)
currentStep <- currentStep + 1
setTxtProgressBar(progressBar, currentStep)

# 步骤7：构建样本数量统计信息
message("步骤7：构建样本数量统计信息...")
sampleInfoLines <- c(
  paste("Number of control samples (con):", numControl),
  paste("Number of treatment samples (tre):", numTreatment)
)
currentStep <- currentStep + 1
setTxtProgressBar(progressBar, currentStep)

# 步骤8：创建输出文件夹，并将结果分别写入到 "Sample Type Matrix.csv" 和 "Sample_Summary.txt"
message("步骤8：写出结果到新建文件夹中...")
outputFolder <- "Final_Results"
if (!dir.exists(outputFolder)) {
  dir.create(outputFolder)
}
# 定义输出文件路径
matrixFile <- file.path(outputFolder, "Sample Type Matrix.csv")
summaryFile <- file.path(outputFolder, "Sample_Summary.txt")

# 写出样本类型矩阵 CSV 文件
write.table(finalMatrix,
            file = matrixFile,
            sep = ",",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)

# 写出样本统计信息 TXT 文件
writeLines(sampleInfoLines, con = summaryFile)

currentStep <- currentStep + 1
setTxtProgressBar(progressBar, currentStep)
# -------------------------- 步骤9：将 "Sample Type Matrix.csv" 转换为 TXT 格式 --------------------------
message("步骤9：将样本类型矩阵保存为 TXT 文件...")
txtFile <- file.path(outputFolder, "Sample Type Matrix.txt")

# 读取 CSV 文件并保存为 TXT 文件
write.table(finalMatrix,
            file = txtFile,
            sep = "\t",   # 使用制表符分隔
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)

message("样本类型矩阵已成功保存为 TXT 文件：", txtFile)

# -------------------------- 完成提示 --------------------------
close(progressBar)
message("所有数据处理步骤已成功完成！请查看文件夹 '", outputFolder, "' 下的：\n",
        "1) 'Sample Type Matrix.csv'\n",
        "2) 'Sample_Summary.txt'")
