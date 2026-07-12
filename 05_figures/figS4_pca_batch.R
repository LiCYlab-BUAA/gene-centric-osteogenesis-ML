# ============================================================
# 合并表达矩阵的 PCA，按“数据集来源”上色 (Supplementary Fig. S5A)
# 用于展示特征空间中确实存在批次结构
# ============================================================
library(ggplot2)

## 1. 读取矩阵：基因 × 样本 (9661 × 223) ---------------------
expr <- read.csv("log2_dataset_frame.csv", header = TRUE,
                 row.names = 1, check.names = FALSE)   # 保留GSM列名原样
expr <- as.matrix(expr)

## 2. ★样本 → 数据集 映射（必须自备）------------------------
# 做一张 sample_dataset_map.csv，两列：Sample,Dataset
#   Sample  = GSM编号（与矩阵列名完全一致）
#   Dataset = 所属数据集（取自 Table S1，如 GSE1367…）
map <- openxlsx::read.xlsx('./sample_dataset_map.xlsx')
dataset_vec <- map$Dataset[match(colnames(expr), map$Sample)]
if (any(is.na(dataset_vec)))
  stop("未匹配到数据集的样本: ",
       paste(colnames(expr)[is.na(dataset_vec)], collapse = ", "))

## 3. 计算 PCA（样本为观测，基因为变量）--------------------
X   <- t(expr)                        # 转置为 样本 × 基因
X   <- X[, apply(X, 2, sd) > 0]       # 去掉零方差基因（scale 要求）
pca <- prcomp(X, center = TRUE, scale. = TRUE)
ve  <- pca$sdev^2 / sum(pca$sdev^2)   # 各主成分方差占比

## 4. 组装绘图数据 ------------------------------------------
plot_df <- data.frame(PC1 = pca$x[, 1],
                      PC2 = pca$x[, 2],
                      Dataset = factor(dataset_vec))

## 5. 34 个数据集需要一个大调色板 --------------------------
n_ds <- nlevels(plot_df$Dataset)
pal  <- if (requireNamespace("pals", quietly = TRUE)){
  unname(pals::polychrome(n_ds))            # 推荐：颜色区分度高
}else{grDevices::hcl.colors(n_ds, "Dark 3")}  # 备用

## 6. 绘图 --------------------------------------------------
p <- ggplot(plot_df, aes(PC1, PC2, colour = Dataset)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_colour_manual(values = pal) +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * ve[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * ve[2]),
       colour = "Dataset") +
  theme_bw(base_size = 12) +
  theme(legend.key.size = unit(0.35, "cm"))

## 7. 输出 --------------------------------------------------
ggsave("PCA_by_dataset.pdf", p, width = 8, height = 5.5)
ggsave("PCA_by_dataset.png", p, width = 8, height = 5.5, dpi = 600)
print(p)