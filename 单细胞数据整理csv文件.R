# 安装所需的包
#if (!requireNamespace("progress", quietly = TRUE)) {
#  install.packages("progress")
#}

# 导入必要的库
library(progress)      # 用于创建进度条
library(Seurat)        # 用于单细胞分析
library(limma)         # 用于差异分析
library(dplyr)         # 用于数据处理
library(magrittr)      # 用于管道操作符（%>%）
library(celldex)       # 用于细胞类型注释
library(SingleR)       # 用于单细胞注释
library(monocle)       # 用于单细胞轨迹分析
library(ggplot2)       # 用于数据可视化
library(RColorBrewer)  # 用于生成调色板
# 2. 设置工作目录
workDir <- "D:/2.椎间盘衰老(1)/Q1区纯生信7.8分/30.单细胞数据整理"  # 设置工作目录为存储单细胞数据的文件夹
setwd(workDir)  # 更改当前工作目录为上述指定的路径

# 3. 读取目录并预处理数据
dirs <- list.dirs(workDir)  # 获取工作目录下的所有子目录
dirs_sample <- dirs[-1]  # 排除第一个父级目录（即工作目录本身）

# 检查是否有样本目录，如果没有则报错
if (length(dirs_sample) == 0) {
  stop("没有找到样本目录，请确认数据路径是否正确！")  # 如果没有找到样本目录，停止程序并报错
}

# 提取子目录名称作为样本的名称
names(dirs_sample) <- gsub(".+\\/(.+)", "\\1", dirs_sample)

# 4. 创建进度条并初始化
progress_bar <- progress::progress_bar$new(
  format = " 进度 [:bar] :percent 完成 时间: :elapsedfull 剩余时间: :eta",  # 设置进度条的格式，显示进度条、百分比、已用时间、剩余时间
  total = length(dirs_sample),  # 设置进度条的总步骤数
  clear = FALSE,  # 完成后保持进度条
  width = 60       # 设置进度条的宽度
)

# 5. 逐一读取每个样本目录的数据
counts <- list()  # 创建空列表用于存储每个样本的基因表达数据

for (i in seq_along(dirs_sample)) {  # 循环处理每个样本的目录
  message(paste("正在读取目录：", dirs_sample[i]))  # 输出当前正在处理的目录
  tryCatch({
    counts[[i]] <- Read10X(data.dir = dirs_sample[i])  # 读取10X单细胞数据
  }, error = function(e) {  # 如果读取过程中发生错误，捕获并输出错误信息
    message(paste("读取目录", dirs_sample[i], "时发生错误：", e$message))
    next  # 跳过当前目录，继续处理下一个目录
  })
  
  progress_bar$tick()  # 每处理完一个目录，更新一次进度条
}

# 6. 合并所有读取的数据
if (length(counts) == 0) {
  stop("没有成功读取任何数据，请检查输入路径和数据格式！")  # 如果没有成功读取任何数据，停止程序并报错
}

counts_combined <- do.call(cbind, counts)  # 将所有样本的基因表达数据按列合并

# 7. 创建Seurat对象，用于后续分析
Peripheral_Blood_Mononuclear_Cells <- CreateSeuratObject(counts_combined, min.cells = 10, min.features = 40)  # 创建Seurat对象，设置最小细胞数和最小基因数

# 8. 获取基因表达数据
single_cell_data <- GetAssayData(Peripheral_Blood_Mononuclear_Cells, slot = "counts")  # 获取原始计数数据

# 9. 输出数据到CSV文件
message("正在导出数据为 CSV 文件...")  # 输出导出提示信息

tryCatch({
  write.csv(as.data.frame(single_cell_data), file = "single_cell_data.csv")  # 将数据框转换为CSV文件
}, error = function(e) {  # 如果导出过程发生错误，捕获并输出错误信息
  message(paste("导出CSV文件时发生错误：", e$message))
})

# 10. 输出处理完成的消息
message("数据已成功导出为 single_cell_data.csv 文件")

# 11. 输出数据矩阵的维度信息（行数和列数）
cat("数据矩阵的维度为：", dim(single_cell_data)[1], "行（基因数）", dim(single_cell_data)[2], "列（细胞数）\n")
