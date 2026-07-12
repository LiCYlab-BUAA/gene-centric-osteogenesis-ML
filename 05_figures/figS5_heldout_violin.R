# ============================================================
# 标签基因的“留出预测概率”小提琴图  —— 纯后处理版
# ============================================================

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
})

# ------------------------------------------------------------
# 0. 参数：结果目录 + 正文报告的 15 个模型标签
# ------------------------------------------------------------
# setwd(".")  # run from the repo root (original absolute path removed)
model_labels <- c("SVM_Linear","SVM_Radial","RandomForest","BaggedTree","RotationForest",
                  "AdaBoost","GBM","LogitBoost","ParallelRF","Ranger","Rborist","PLS",
                  "Multinomial","WSRF","PCA_NNet")   # 与 Fig 2/3j 一致的 15 个

# ------------------------------------------------------------
# 1. 解析单个模型 CSV -> 每基因平均留出概率 + 真实标签
# ------------------------------------------------------------
parse_model <- function(lab) {
  f <- sprintf("OB_Model_result_frame_%s.csv", lab)
  if (!file.exists(f)) { message("  缺文件，跳过: ", f); return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  
  # 按行拆分(行内对齐有保证)，再 unlist
  gl <- strsplit(as.character(d$Test_genes),  ", ", fixed = TRUE)
  pl <- lapply(strsplit(as.character(d$Probs), ", ", fixed = TRUE), as.numeric)
  yl <- lapply(strsplit(as.character(d$Test_labels), ", ", fixed = TRUE), as.numeric)
  ok <- lengths(gl) == lengths(pl) & lengths(gl) == lengths(yl)   # 丢弃异常行
  if (any(!ok)) message(sprintf("  %s: 跳过 %d 个长度不齐的行", lab, sum(!ok)))
  g <- unlist(gl[ok]); p <- unlist(pl[ok]); y <- unlist(yl[ok])
  
  list(prob = tapply(p, g, mean, na.rm = TRUE),           # 每基因平均留出概率
       lab  = tapply(y, g, function(x) x[1]))             # 每基因真实标签(0/1)
}

res <- setNames(lapply(model_labels, parse_model), model_labels)
res <- res[!sapply(res, is.null)]
message(sprintf("成功解析 %d 个模型: %s", length(res), paste(names(res), collapse = ", ")))

# ------------------------------------------------------------
# 2. 组装 基因 × 模型 的留出概率矩阵
# ------------------------------------------------------------
all_genes <- sort(unique(unlist(lapply(res, function(x) names(x$prob)))))
mean_hold <- sapply(res, function(x) x$prob[all_genes])
rownames(mean_hold) <- all_genes

# 真实标签(取自 Test_labels，与实际评分一致)
lab_vec  <- res[[1]]$lab[all_genes]
true_lab <- factor(ifelse(lab_vec == 1, "POSITIVE", "NEGATIVE"),
                   levels = c("POSITIVE", "NEGATIVE"))
message(sprintf("标签基因: %d POSITIVE + %d NEGATIVE = %d",
                sum(true_lab == "POSITIVE"), sum(true_lab == "NEGATIVE"), length(true_lab)))

# ------------------------------------------------------------
# 3. 长表 + pooled(跨模型平均)
# ------------------------------------------------------------
long <- data.frame(
  Gene  = rep(all_genes, times = ncol(mean_hold)),
  Model = rep(colnames(mean_hold), each = nrow(mean_hold)),
  Prob  = as.vector(mean_hold),
  stringsAsFactors = FALSE
)
long$Label <- true_lab[match(long$Gene, all_genes)]
long <- long[!is.na(long$Prob), ]

gene_mean <- rowMeans(mean_hold, na.rm = TRUE)
pooled <- data.frame(Gene = all_genes, Prob = gene_mean, Label = true_lab)

# 导出概率表(可并入给审稿人的 comprehensive table)
write.csv(data.frame(Gene = all_genes, Label = as.character(true_lab),
                     round(mean_hold, 3), Ensemble_mean = round(gene_mean, 3),
                     check.names = FALSE),
          "heldout_prob_labeled_genes.csv", row.names = FALSE)

# ------------------------------------------------------------
# 4. 分离度统计(回答“分数是否重叠严重”)
# ------------------------------------------------------------
w <- wilcox.test(Prob ~ Label, data = pooled)
auc_pool <- as.numeric(pROC::auc(pROC::roc(
  response = pooled$Label, predictor = pooled$Prob,
  levels = c("NEGATIVE", "POSITIVE"), direction = "<", quiet = TRUE)))
message(sprintf("Pooled held-out separation: AUC=%.3f, Wilcoxon p=%.2e", auc_pool, w$p.value))

# ------------------------------------------------------------
# 5. 绘图(红=Positive, 灰=Negative，与 Fig3 一致)
# ------------------------------------------------------------
cols <- c(POSITIVE = "#E8756A", NEGATIVE = "#BEBEBE")

p_pool <- ggplot(pooled, aes(Label, Prob, fill = Label)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "Held-out P(OB)",
       subtitle = sprintf("Ensemble of %d models | Wilcoxon p = %.1e",
                          ncol(mean_hold), w$p.value)) +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave("FigS_heldout_violin_pooled.pdf", p_pool, width = 4.2, height = 4.5)
ggsave("FigS_heldout_violin_pooled.svg", p_pool, width = 4.2, height = 4.5)

message("Done: FigSX_heldout_violin_pooled.pdf, FigSX_heldout_violin_byModel.pdf, heldout_prob_labeled_genes.csv")