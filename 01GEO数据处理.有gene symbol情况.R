# 1. 设置工作目录
setwd("I:\\2025教程新n\\14.Q1区7.8分，MR+单细胞\\03.GEO数据整理数据整理")  # 设置工作目录

# 2. 设置参数：用户指定基因信息所在行数（例如第11行）并转换为0起始索引
inputRow <- 13                              # 用户指定使用第?行数据（1起始）
targetColIdx <- inputRow - 1                  # 转换为0起始索引，用于平台文件中目标列定位

# 3. 定义数据文件及输出文件的路径
exprFilePath <- "GSE.txt"                     # 表达数据文件路径
platformFilePath <- "GPL.txt"                 # 平台文件路径（含探针与基因信息）
outputFilePath <- "geneMatrix.txt"            # 最终输出文件路径

# 4. 读取表达数据文件（GSE.txt）
cat("正在加载表达数据文件：", exprFilePath, "...\n")  # 输出提示信息
exprData <- read.delim(exprFilePath,          # 读取文件
                       header = TRUE,          # 第一行为表头
                       sep = "\t",             # 使用制表符分隔
                       quote = "\"",           # 使用双引号包裹字符型数据
                       comment.char = "!")     # 忽略以"!"开头的注释行
colnames(exprData)[1] <- "ProbeID"            # 将第一列重命名为 "ProbeID"
cat("表达数据加载完成，共", ncol(exprData)-1, "个样本。\n")  # 输出样本数（减去探针ID列）

# 5. 读取平台文件（GPL.txt）
cat("正在加载平台文件：", platformFilePath, "...\n")  # 输出提示信息
platformData <- read.delim(platformFilePath,    # 读取平台文件
                           header = FALSE,        # 无表头
                           sep = "\t",            # 以制表符分隔
                           quote = "\"",          # 使用双引号包裹字符型数据
                           comment.char = "#",    # 忽略以"#"开头的注释行
                           stringsAsFactors = FALSE)  # 不将字符转换为因子

# 6. 构建探针与基因符号映射（逐行处理并显示进度）
cat("开始处理平台文件数据...\n")         # 输出提示信息
totalRows <- nrow(platformData)               # 获取平台文件总行数
geneMapping <- list()                         # 初始化空列表用于存储映射：探针ID -> 基因符号
pbPlat <- txtProgressBar(min = 0, max = totalRows, style = 3)  # 创建进度条，范围从0到总行数

for (row in 1:totalRows) {                     # 遍历平台文件的每一行
  # 若整行全为空，或第一列以 "ID" 或 "!" 开头，则跳过该行
  if (all(platformData[row, ] == "") || grepl("^(ID|\\!)", platformData[row, 1])) {
    setTxtProgressBar(pbPlat, row)             # 更新进度条
    next
  }
  
  # 检查当前行是否包含目标列（targetColIdx + 1）
  if (ncol(platformData) >= (targetColIdx + 1)) {
    rawGene <- platformData[row, targetColIdx + 1]  # 提取目标列的原始基因信息
    # 如果原始信息非空且不包含空格（只包含一个单词）
    if (rawGene != "" && !grepl(".+\\s+.+", rawGene)) {
      # 若包含"///"，只取左侧部分；同时去除引号
      cleanGene <- sub("(.+?)///(.+)", "\\1", rawGene)
      cleanGene <- gsub('"', '', cleanGene)
      # 建立映射：将当前行第一列（探针ID）映射到处理后的基因符号
      geneMapping[[ as.character(platformData[row, 1]) ]] <- cleanGene
    }
  }
  setTxtProgressBar(pbPlat, row)             # 更新进度条
}
close(pbPlat)                                 # 关闭进度条
cat("\n平台文件处理完成，成功建立映射数：", length(geneMapping), "\n")  # 输出映射建立完成信息

# 7. 合并表达数据与平台映射数据
cat("正在合并表达数据与平台映射数据...\n")  # 输出提示信息
# 从表达数据中提取探针ID列并保存
probeIDs <- exprData[["ProbeID"]]
# 取出除探针ID以外的表达数据部分
exprDataCore <- exprData[, -1, drop = FALSE]
# 构造数据框，用于后续合并：将ProbeID与表达数据一起使用
exprDataWithID <- exprData
# 将平台映射转换为数据框，包含探针ID和对应的基因符号
mappingDF <- data.frame(ProbeID = names(geneMapping),
                        geneSymbol = unlist(geneMapping),
                        stringsAsFactors = FALSE)
# 利用 merge 函数按"ProbeID"匹配，合并表达数据与基因映射信息
mergedData <- merge(exprDataWithID, mappingDF, by = "ProbeID")
cat("合并完成，共获得", nrow(mergedData), "行有效数据。\n")  # 输出合并后数据行数

# 8. 按基因分组，计算各样本的平均表达值（逐基因处理显示进度）
cat("正在按基因分组并计算各样本平均表达值...\n")  # 输出提示信息
# 获取所有样本列名称（排除ProbeID和新合并的geneSymbol列）
sampleCols <- setdiff(colnames(mergedData), c("ProbeID", "geneSymbol"))
# 根据基因符号对数据进行分组
geneGroups <- split(mergedData, mergedData$geneSymbol)
# 获取所有不同基因名称
geneNamesUnique <- names(geneGroups)
totalGenes <- length(geneNamesUnique)         # 总基因数
resultList <- vector("list", totalGenes)        # 初始化列表保存计算结果
pbAgg <- txtProgressBar(min = 0, max = totalGenes, style = 3)  # 创建聚合进度条

for (i in seq_along(geneNamesUnique)) {         # 遍历每个基因组
  geneName <- geneNamesUnique[i]                # 当前基因名称
  groupData <- geneGroups[[i]]                  # 取出该基因对应的所有探针数据
  # 对样本列（除ProbeID和geneSymbol）计算均值，忽略缺失值
  means <- colMeans(groupData[, sampleCols, drop = FALSE], na.rm = TRUE)
  # 将当前基因名称与计算的均值合并为一行
  resultList[[i]] <- c(geneSymbol = geneName, means)
  setTxtProgressBar(pbAgg, i)                   # 更新进度条
}
close(pbAgg)                                   # 关闭聚合进度条
# 将结果列表合并为数据框
aggData <- do.call(rbind, resultList)
# 转换为数据框，并保持字符型基因名称
aggData <- as.data.frame(aggData, stringsAsFactors = FALSE)
# 对除geneSymbol外的数值型列转换为数值型（因do.call返回字符型）
for (col in colnames(aggData)[-1]) {
  aggData[[col]] <- as.numeric(aggData[[col]])
}
# 对结果按照基因名称排序
aggData <- aggData[order(aggData$geneSymbol), ]
cat("数据聚合完成，共计算出", nrow(aggData), "个基因的平均表达值。\n")  # 输出聚合结果

# 9. 写入输出文件（geneMatrix.txt）
cat("正在写入输出文件：", outputFilePath, "...\n")  # 输出写入提示信息
# 利用 write.table 写出数据，使用制表符分隔、不写行号、不加引号
write.table(aggData, file = outputFilePath, sep = "\t", row.names = FALSE, quote = FALSE)
cat("输出文件生成成功：", outputFilePath, "\n")  # 输出完成提示
