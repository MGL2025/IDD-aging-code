# 1. 设置参数
filter_pvalue = 0.05    # 设置筛选用的p值阈值 (0.05, 0.01, 0.001)
mr_data_file = "all_odds_ratios.csv"        # 存储MR分析结果的文件
pleiotropy_file = "all_pleiotropy_results.csv"     # 存储多效性分析结果的文件
heterogeneity_file = "all_heterogeneity_results.csv"  # 存储异质性分析结果的文件
output_file = "IVWfilter.csv"  # 输出筛选结果的文件

# 2. 加载所需的库
library(progress)  # 用于显示进度条
library(dplyr)     # 用于数据操作

# 3. 设置工作目录
setwd("D:/2.椎间盘衰老(1)/Q1区纯生信7.8分/13.筛选有意义的暴露因素")  # 设置工作目录路径

# 4. 读取MR结果数据文件
mr_results = read.csv(mr_data_file, header = TRUE, sep = ",", check.names = FALSE)

# 5. 只选择逆方差加权（IVW）法的结果
ivw_results = mr_results[mr_results$method == "Inverse variance weighted", ]

# 6. 根据设定的p值阈值筛选IVW法结果
ivw_results_filtered = ivw_results[ivw_results$pval < filter_pvalue, ]

# 7. 打印经过IVW筛选后的数据行数
print(paste("IVW Method, rows after p-value filter: ", nrow(ivw_results_filtered)))

# 8. 根据OR值筛选数据，保留OR值均大于1或均小于1的暴露变量
ivw_exposures = data.frame()  # 用来存储筛选后的IVW数据
for (exposure in unique(ivw_results_filtered$exposure)) {
  exposure_data = mr_results[mr_results$exposure == exposure, ]
  if (sum(exposure_data$or > 1) == nrow(exposure_data) | sum(exposure_data$or < 1) == nrow(exposure_data)) {
    ivw_exposures = rbind(ivw_exposures, ivw_results_filtered[ivw_results_filtered$exposure == exposure, ])
  }
}

# 9. 打印筛选OR值后的数据行数
print(paste("Rows after OR filter: ", nrow(ivw_exposures)))

# 10. 读取多效性结果数据
pleiotropy_results = read.csv(pleiotropy_file, header = TRUE, sep = ",", check.names = FALSE)

# 11. 根据p值筛选多效性结果，选择p值大于0.05的暴露变量
pleiotropy_results_filtered = pleiotropy_results[pleiotropy_results$pval > 0.05, ]
exposure_list = as.vector(pleiotropy_results_filtered$exposure)

# 12. 打印经过多效性筛选后的数据行数
print(paste("Rows in pleiotropy results after filtering p-value > 0.05: ", nrow(pleiotropy_results_filtered)))

# 13. 读取异质性分析结果数据
heterogeneity_results = read.csv(heterogeneity_file, header = TRUE, sep = ",", check.names = FALSE)

# 14. 根据Q_pval筛选异质性分析结果，保留p值大于0.05的结果
heterogeneity_results_filtered = heterogeneity_results[heterogeneity_results$Q_pval > 0.05, ]

# 15. 打印经过异质性筛选后的数据行数
print(paste("Rows in heterogeneity results after filtering Q_pval > 0.05: ", nrow(heterogeneity_results_filtered)))

# 16. 提取同时出现在多效性和异质性分析结果中的暴露变量
matching_exposures = exposure_list[exposure_list %in% heterogeneity_results_filtered$exposure]

# 17. 打印筛选出符合条件的暴露变量数量
print(paste("Number of exposures matching between pleiotropy and heterogeneity results: ", length(matching_exposures)))

# 18. 筛选MR结果中，符合异质性条件的暴露变量
final_results = ivw_exposures[ivw_exposures$exposure %in% matching_exposures, ]

# 19. 打印最终输出结果的行数
print(paste("Rows in final output table: ", nrow(final_results)))

# 20. 判断如果结果为空，给出提示
if (nrow(final_results) == 0) {
  print("没有找到符合条件的结果，请检查暴露变量是否匹配。")
} else {
  # 如果结果不为空，保存筛选后的结果到CSV文件
  write.csv(final_results, file = output_file, row.names = FALSE)
  print("筛选后的结果已保存。")
}

# 21. 添加进度条以显示数据处理进度
pb <- progress_bar$new(
  format = "处理进度 [:bar] :percent",
  total = 100, clear = FALSE, width = 60
)

for (i in 1:100) {
  Sys.sleep(0.04)  
  pb$tick()
}

# 23. 完成所有处理
print("所有步骤完成。")
