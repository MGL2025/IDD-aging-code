
setwd("D:/2.椎间盘衰老/Q1区纯生信7.8分/14.提取有意义的暴露因素")  # 设置工作目录路径

# 1. 设置文件路径
ivw_filter_file <- "IVWfilter.csv"  # IVW筛选结果文件
pqtl_file <- "eqtl.txt"             # pqtl文件

# 2. 读取IVWfilter.csv数据，提取exposure列的基因名字
ivw_data <- read.csv(ivw_filter_file, header = TRUE, sep = ",", check.names = FALSE)
genes_in_ivw <- unique(ivw_data$exposure)  # 提取唯一的基因名

# 3. 读取pqtl.txt文件
pqtl_data <- read.table(pqtl_file, header = TRUE, sep = "\t", check.names = FALSE)

# 4. 筛选pqtl文件中exposure列为IVWfilter.csv基因名字的行
pqtl_filtered <- pqtl_data[pqtl_data$exposure %in% genes_in_ivw, ]

# 5. 输出筛选后的结果到新的文件
write.table(pqtl_filtered, file = "eqtl_filtered.txt", sep = "\t", row.names = FALSE, quote = FALSE)
