# ----------------------------- #
# 第1步：加载所需 R 包
# ----------------------------- #
library(limma)     # 用于 avereps 等去重函数
library(GSVA)      # 用于 ssGSEA 分析
library(GSEABase)  # 用于读取 GMT 格式的基因集
library(ggpubr)    # 美观的统计作图
library(reshape2)  # 数据重塑工具，如 melt/dcast
# ----------------------------- #
# 第2步：设置文件路径与工作目录
# ----------------------------- #
expr_csv  <- "Sample Type Matrix.csv"         
gmt_path  <- "immune.gmt"                      
wd_target <- "D:/2.椎间盘衰老(1)/Q1区纯生信7.8分/23.免疫浸润评分GSVA"

if (!dir.exists(wd_target)) {
  stop("工作目录不存在：", wd_target)
}
setwd(wd_target)

if (!file.exists(expr_csv) || !file.exists(gmt_path)) {
  stop("未找到表达矩阵或 GMT 文件：", expr_csv, gmt_path)
}

# ----------------------------- #
# 第3步：初始化进度条及步骤计数
# ----------------------------- #
total_steps  <- 5                              # 总共 5 个步骤
step_counter <- 0
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)

# ----------------------------- #
# 步骤1：读取并整理表达数据
# ----------------------------- #
step_counter <- step_counter + 1
setTxtProgressBar(pb, step_counter)

raw_df <- read.table(expr_csv, header = TRUE, sep = ",",
                     check.names = FALSE,
                     stringsAsFactors = FALSE)     # 读入原始表
rownames(raw_df) <- raw_df[[1]]                   # 第一列设为行名
expr_df <- raw_df[, -1, drop = FALSE]             # 删除 ID 列
expr_mat <- apply(expr_df, 2, as.numeric)         # 强制数值矩阵
rownames(expr_mat) <- rownames(raw_df)            # 恢复行名

# ----------------------------- #
# 步骤2：去重 & 过滤低表达
# ----------------------------- #
step_counter <- step_counter + 1
setTxtProgressBar(pb, step_counter)

expr_numeric <- as.matrix(expr_mat)                # 确保矩阵
if (is.null(rownames(expr_numeric))) {
  stop("表达矩阵缺少行名")
}

# 优先使用 limma::avereps 去重
if ("avereps" %in% ls("package:limma")) {
  dedup_mat <- limma::avereps(expr_numeric)        # 去重均值
} else {
  dedup_simple <- function(mtx) {
    gids <- rownames(mtx)
    uids <- unique(gids)
    mat_list <- lapply(uids, function(g) {
      sub <- mtx[gids == g, , drop = FALSE]
      colMeans(sub, na.rm = TRUE)
    })
    out <- do.call(rbind, mat_list)
    rownames(out) <- uids
    out
  }
  dedup_mat <- dedup_simple(expr_numeric)
}

# 过滤：保留行平均 > 0 的基因
keep       <- rowMeans(dedup_mat, na.rm = TRUE) > 0
expr_filtered <- dedup_mat[keep, , drop = FALSE]

# ----------------------------- #
# 步骤3：加载基因集并运行 ssGSEA
# ----------------------------- #
step_counter <- step_counter + 1
setTxtProgressBar(pb, step_counter)

gene_sets <- getGmt(gmt_path, geneIdType = SymbolIdentifier())  # 载入 GMT

# 新版 GSVA：直接用最简 ssgseaParam 构造，默认 minSize/maxSize/alpha/normalize
ssgsea_prm <- ssgseaParam(
  exprData = expr_filtered,
  geneSets = gene_sets
)  # :contentReference[oaicite:0]{index=0}

# 并行参数（此处串行示例）
bp_param <- BiocParallel::SerialParam()

# 调用 gsva，触发 ssGSEA 方法
ssgsea_res <- gsva(
  ssgsea_prm,
  verbose = TRUE,
  BPPARAM = bp_param
)  # :contentReference[oaicite:1]{index=1}

stopifnot(is.matrix(ssgsea_res))  # 强制检查输出类型

# ----------------------------- #
# 步骤4：0–1 归一化
# ----------------------------- #
step_counter <- step_counter + 1
setTxtProgressBar(pb, step_counter)

zero_one_norm <- function(v) {
  rng <- range(v, na.rm = TRUE)
  (v - rng[1]) / (rng[2] - rng[1])
}

norm_scores <- apply(ssgsea_res, 2, zero_one_norm)              # 按列归一化
norm_scores <- rbind(DUMMY = rep(NA, ncol(norm_scores)), norm_scores) # 冗余示例

# ----------------------------- #
# 步骤5：保存结果 & 收尾
# ----------------------------- #
step_counter <- step_counter + 1
setTxtProgressBar(pb, step_counter)

output_mat <- rbind(SampleID = colnames(norm_scores), norm_scores)
write.table(output_mat,
            file      = "GSVA_Scores.csv",
            sep       = ",",
            quote     = FALSE,
            col.names = FALSE)

close(pb)
message("全部 ", total_steps, " 步骤完成，结果已保存至 ssGSEA_Scores.csv")
