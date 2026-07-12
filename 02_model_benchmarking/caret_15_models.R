# setwd('.')  # run from the repo root (original absolute path removed)
# ============================================================
# OB 基因分类模型 —— caret 多模型版本
# 原始结构完全保留，仅替换建模核心部分
# ============================================================

# 0. 加载包 ---------------------------------------------------
library(stringr)
library(dplyr)
library(pROC)
library(caret)
library(parallel)
library(doParallel)

# 需要的模型依赖包（首次运行会自动安装缺失的）
required_pkgs <- c(
  "randomForest",   # rf
  "gbm",            # gbm
  "e1071",          # svmRadial, svmLinear, naiveBayes
  "kernlab",        # svmRadial2
  "rpart",          # rpart
  "C50",            # C5.0
  "xgboost",        # xgbTree
  "earth",          # bagEarth, earth
  "mda",            # mda
  "klaR",           # nb
  "pls",            # pls
  "elasticnet",     # enet
  "LogicReg",       # logreg
  "ada",            # ada
  "plyr"            # needed by several
)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch(install.packages(pkg, repos = "https://cloud.r-project.org",
                              quiet = TRUE),
             error = function(e) message("  跳过安装: ", pkg))
  }
}

# ============================================================
# 准备工作（原始代码不变）
# ============================================================
{
  # 读取数据
  exp_change <- read.csv('log2_dataset_frame.csv', header = TRUE)
  colnames(exp_change)[1] <- 'Symbol'
  all_genes <- unique(exp_change$Symbol)
  
  # 建立标签集
  exp_change$Label <- 'Unknown'
  Labels <- read.csv('Label_genes.csv', header = TRUE)
  OB_gene_all             <- Labels$Symbol[which(Labels$Label == 1)]
  raw_negative_genes      <- Labels$Symbol[which(Labels$Label == 0)]
  exp_change$Label[which(exp_change$Symbol %in% OB_gene_all)] <- 1
  ava_OB_genes            <- all_genes[which(all_genes %in% OB_gene_all)]
  ava_selected_negative_genes <- intersect(all_genes, toupper(raw_negative_genes))
  exp_change$Label[which(exp_change$Symbol %in% ava_selected_negative_genes)] <- 0
  
  # 可以抽取的训练测试矩阵
  Train_test_frame          <- exp_change[which(exp_change$Label != 'Unknown'), ]
  Train_test_frame$Label    <- factor(Train_test_frame$Label, levels = c(0, 1))
  ava_genes                 <- Train_test_frame$Symbol
  rownames(Train_test_frame) <- Train_test_frame$Symbol
  Train_test_frame          <- Train_test_frame[, -1]
  
  # 训练集和测试集的基因数量
  OB_train_num  <- round(length(ava_OB_genes) * 0.8)
  OB_test_num   <- length(ava_OB_genes) - OB_train_num
  Fat_train_num <- round(length(ava_selected_negative_genes) * 0.8)
  Fat_test_num  <- length(ava_selected_negative_genes) - Fat_train_num
  
  # 构建必要的其他函数
  {
    numbers_to_char <- function(x) {
      temp_char <- ''
      for (i1 in 1:length(x)) {
        temp_char <- paste0(temp_char, ', ', x[i1])
      }
      final_result <- str_sub(temp_char, 3, nchar(temp_char))
      return(final_result)
    }
    char_to_numbers <- function(x) {
      temp_numbers <- c()
      temp_sep     <- str_locate_all(x, ', ')
      temp_start   <- c(1, c(as.numeric(temp_sep[[1]][, 2]) + 1))
      temp_end     <- c((as.numeric(temp_sep[[1]][, 1]) - 1), nchar(x))
      for (i2 in 1:length(temp_start)) {
        temp_numbers <- c(
          temp_numbers,
          str_sub(x, temp_start[i2], temp_end[i2])
        )
      }
      temp_numbers <- as.numeric(temp_numbers)
      return(temp_numbers)
    }
    judge_anwser <- function(answer, to_check) {
      assign('result', list())
      FP <- length(which(to_check == 1 & answer == 0))
      TP <- length(which(to_check == 1 & answer == 1))
      FN <- length(which(to_check == 0 & answer == 1))
      TN <- length(which(to_check == 0 & answer == 0))
      result[[1]] <- FP / (TN + FP)   # FPR
      result[[2]] <- TP / (TP + FN)   # TPR
      result[[3]] <- FN / (TP + FN)   # FNR
      result[[4]] <- TN / (TN + FP)   # TNR
      names(result) <- c('FPR', 'TPR', 'FNR', 'TNR')
      return(result)
    }
    MAE <- function(actual, predicted) {
      mean(abs(actual - predicted))
    }
  }
}

# ============================================================
# 模型列表（caret method 名称）
# 每个元素：list(method, needs_prob)
# ============================================================
model_list <- list(
  list(method = "rf",          label = "RandomForest"),
  list(method = "gbm",         label = "GBM"),
  list(method = "svmRadial",   label = "SVM_Radial"),
  list(method = "svmLinear",   label = "SVM_Linear"),
  list(method = "rpart",       label = "DecisionTree"),
  list(method = "C5.0",        label = "C50"),
  list(method = "xgbTree",     label = "XGBoost"),
  list(method = "earth",       label = "MARS"),
  list(method = "qda",         label = "QDA"),
  list(method = "nb",          label = "NaiveBayes"),
  list(method = "knn",         label = "KNN_caret"),
  list(method = "pls",         label = "PLS"),
  list(method = "ada",         label = "AdaBoost"),
  list(method = "bayesglm",    label = "BayesGLM"),
  list(method = "glm",         label = "LogisticReg"),
  list(method = "multinom",    label = "Multinomial"),
  list(method = "bagEarth",    label = "BaggedEarth"),
  list(method = "treebag",     label = "BaggedTree"),
  list(method = "parRF",       label = "ParallelRF"),
  list(method = "LogitBoost",  label = "LogitBoost"),
  list(method = "enet",        label = "ElasticNet"),
  list(method = "pcaNNet",     label = "PCA_NNet"),
  list(method = "hdda",        label = "HDDA"),
  list(method = "wsrf",        label = "WSRF"),
  list(method = "ranger",      label = "Ranger"),
  list(method = "Rborist",     label = "Rborist"),
  list(method = "rotationForest", label = "RotationForest")
)

# ============================================================
# 输出目录
# ============================================================
cycles_dir <- "cycles"
if (!dir.exists(cycles_dir)) dir.create(cycles_dir)

# ============================================================
# caret 训练控制参数
# ============================================================
train_ctrl <- trainControl(
  method          = "none",        # 不做内部CV，加速循环
  classProbs      = TRUE,          # 输出概率
  summaryFunction = twoClassSummary,
  savePredictions = FALSE
)

# 把 0/1 factor 转为 caret 喜欢的字母标签
relabel <- function(fac) {
  factor(ifelse(as.character(fac) == "1", "OB", "N"), levels = c("N", "OB"))
}

# ============================================================
# 建立空白结果表（每个模型一张，存在同名 CSV）
#
# 新增列说明：
#   Train_genes      —— 训练集基因名（逗号分隔字符串）
#   Test_genes       —— 测试集基因名（逗号分隔字符串）
#   Train_rows       —— 训练集行索引（与原始代码一致）
#   Probs            —— 测试集预测概率
#   Test_labels      —— 测试集真实标签
#   Model            —— 模型名称
#   Cycle            —— 循环编号
#   AUC              —— ROC-AUC
#   Accuracy         —— 整体准确率
#   Sensitivity/TPR  —— 真阳性率（召回率）
#   Specificity/TNR  —— 真阴性率
#   FPR              —— 假阳性率
#   FNR              —— 假阴性率
#   Precision        —— 精确率 TP/(TP+FP)
#   F1               —— F1 分数
#   MCC              —— Matthews相关系数
#   Kappa            —— Cohen's Kappa
#   number_1         —— 训练集正样本数
#   number_0         —— 训练集负样本数
#   Cor              —— 预测概率与真实标签的相关系数
#   MAE              —— 均方误差
#   ENPEP_Result     —— ENPEP 基因的预测概率
# ============================================================
result_cols <- c(
  'Cycle',
  'Train_genes', 'Test_genes', 'Train_rows',
  'Probs', 'Test_labels',
  'Model',
  'AUC', 'Accuracy',
  'Sensitivity', 'Specificity',
  'FPR', 'FNR',
  'Precision', 'F1', 'MCC', 'Kappa',
  'number_1', 'number_0',
  'Cor', 'MAE',
  'ENPEP_Result'
)

# MCC 计算函数
calc_MCC <- function(TP, TN, FP, FN) {
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (denom == 0) return(NA_real_)
  (TP * TN - FP * FN) / denom
}

# Cohen's Kappa
calc_Kappa <- function(TP, TN, FP, FN) {
  n    <- TP + TN + FP + FN
  po   <- (TP + TN) / n
  pe   <- ((TP + FP) / n) * ((TP + FN) / n) +
    ((TN + FN) / n) * ((TN + FP) / n)
  if (pe == 1) return(NA_real_)
  (po - pe) / (1 - pe)
}

# 把基因名向量转为逗号分隔字符串
genes_to_char <- function(gene_vec) {
  paste(gene_vec, collapse = ", ")
}

# 每个模型独立 data.frame，存在一个命名 list 中
all_results <- lapply(model_list, function(m) {
  df <- data.frame(matrix(nrow = 0, ncol = length(result_cols)))
  colnames(df) <- result_cols
  df
})
names(all_results) <- sapply(model_list, `[[`, "label")

# ============================================================
# 提取 ENPEP 数据（原始代码不变）
# ============================================================
ENPEP_data <- exp_change[which(exp_change$Symbol == 'ENPEP'), ]
rownames(ENPEP_data) <- 'ENPEP'
ENPEP_data <- ENPEP_data[, -1]
ENPEP_data <- ENPEP_data[, -which(colnames(ENPEP_data) == 'Label')]

# ============================================================
# 大循环：随机采样
# ============================================================
All_cycle_times <- 1000   # 与原始代码一致

# 可选：多核加速（去掉注释启用）
# cl <- makeCluster(max(1, detectCores() - 2))
# registerDoParallel(cl)

message("开始训练，共 ", All_cycle_times, " 次循环 × ",
        length(model_list), " 个模型...")

for (cycle_times in 1:All_cycle_times) {
  
  # —— 随机选取训练/测试基因（与原始代码完全相同）——
  OB_rows <- sample(which(Train_test_frame$Label == 1),
                    OB_train_num, replace = FALSE)
  N_rows  <- sample(which(Train_test_frame$Label == 0),
                    Fat_train_num, replace = FALSE)
  
  train_frame     <- Train_test_frame[c(OB_rows, N_rows), ]
  ob_train_label  <- train_frame$Label
  # 记录训练集基因名
  train_gene_names <- rownames(train_frame)
  train_frame     <- train_frame[, -which(colnames(train_frame) == 'Label')]
  
  test_frame      <- Train_test_frame[-c(OB_rows, N_rows), ]
  ob_test_label   <- test_frame$Label
  # 记录测试集基因名
  test_gene_names  <- rownames(test_frame)
  test_frame      <- test_frame[, -which(colnames(test_frame) == 'Label')]
  
  # caret 标签（字母）
  train_label_chr <- relabel(ob_train_label)
  test_label_chr  <- relabel(ob_test_label)
  
  # —— 遍历每个模型 ——
  for (mi in seq_along(model_list)) {
    method_name <- model_list[[mi]]$method
    label_name  <- model_list[[mi]]$label
    
    # 跳过依赖包未安装的模型
    pkg_ok <- tryCatch({
      getModelInfo(method_name)[[1]]  # 检测 caret 是否认识这个 method
      TRUE
    }, error = function(e) FALSE)
    if (!pkg_ok) next
    
    # 训练
    fitted_model <- tryCatch({
      suppressWarnings(
        train(
          x         = train_frame,
          y         = train_label_chr,
          method    = method_name,
          trControl = train_ctrl,
          metric    = "ROC",
          # 关闭冗余输出（gbm 等）
          verbose   = FALSE
        )
      )
    }, error = function(e) NULL)
    
    if (is.null(fitted_model)) next  # 训练失败则跳过
    
    # 预测测试集概率
    prob_pred <- tryCatch({
      p <- predict(fitted_model, newdata = test_frame, type = "prob")
      p[, "OB"]   # 取 OB 类的概率
    }, error = function(e) NULL)
    
    if (is.null(prob_pred)) next
    
    # AUC
    ROC_frame <- data.frame(
      pre    = prob_pred,
      answer = as.numeric(as.vector(ob_test_label))
    )
    auc_val <- tryCatch({
      p1 <- roc(answer ~ pre, ROC_frame, levels = c(0, 1), direction = '<')
      as.numeric(p1$auc)
    }, error = function(e) NA_real_)
    
    # ENPEP 预测
    ENPEP_prob <- tryCatch({
      pe <- predict(fitted_model, newdata = ENPEP_data, type = "prob")
      pe[, "OB"]
    }, error = function(e) {
      tryCatch({
        as.character(predict(fitted_model, newdata = ENPEP_data))
      }, error = function(e2) NA)
    })
    
    # —— 效能评价（基于最优阈值 0.5）——
    pred_class   <- ifelse(prob_pred >= 0.5, 1, 0)
    true_class   <- as.numeric(as.vector(ob_test_label))
    
    TP  <- sum(pred_class == 1 & true_class == 1)
    TN  <- sum(pred_class == 0 & true_class == 0)
    FP  <- sum(pred_class == 1 & true_class == 0)
    FN  <- sum(pred_class == 0 & true_class == 1)
    
    Accuracy    <- (TP + TN) / (TP + TN + FP + FN)
    Sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_  # TPR
    Specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_  # TNR
    FPR_val     <- if ((TN + FP) > 0) FP / (TN + FP) else NA_real_
    FNR_val     <- if ((TP + FN) > 0) FN / (TP + FN) else NA_real_
    Precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    F1_val      <- if (!is.na(Precision) & !is.na(Sensitivity) &
                       (Precision + Sensitivity) > 0)
      2 * Precision * Sensitivity / (Precision + Sensitivity) else NA_real_
    MCC_val     <- calc_MCC(TP, TN, FP, FN)
    Kappa_val   <- calc_Kappa(TP, TN, FP, FN)
    
    # 写入该模型的结果 data.frame
    new_row <- c(
      cycle_times,                                          # Cycle
      genes_to_char(train_gene_names),                     # Train_genes
      genes_to_char(test_gene_names),                      # Test_genes
      numbers_to_char(c(OB_rows, N_rows)),                 # Train_rows
      numbers_to_char(prob_pred),                          # Probs
      numbers_to_char(true_class),                         # Test_labels
      label_name,                                          # Model
      auc_val,                                             # AUC
      Accuracy,                                            # Accuracy
      Sensitivity,                                         # Sensitivity/TPR
      Specificity,                                         # Specificity/TNR
      FPR_val,                                             # FPR
      FNR_val,                                             # FNR
      Precision,                                           # Precision
      F1_val,                                              # F1
      MCC_val,                                             # MCC
      Kappa_val,                                           # Kappa
      length(which(ob_train_label == 1)),                  # number_1
      length(which(ob_train_label == 0)),                  # number_0
      tryCatch(cor(prob_pred, true_class),
               error = function(e) NA),                   # Cor
      MAE(true_class, prob_pred),                         # MAE
      as.character(ENPEP_prob)                             # ENPEP_Result
    )
    
    nr <- nrow(all_results[[label_name]]) + 1
    all_results[[label_name]][nr, ] <- new_row
    
    # —— 保存模型到 cycles/ 文件夹 ——
    model_save_path <- file.path(
      cycles_dir,
      sprintf("%s_cycle%05d.rds", label_name, cycle_times)
    )
    tryCatch(saveRDS(fitted_model, model_save_path), error = function(e) NULL)
    
  }  # end model loop
  
  # 每 100 次循环输出进度 + 写中间结果
  if (cycle_times %% 100 == 0) {
    message(Sys.time(), "  完成 ", cycle_times, " / ", All_cycle_times, " 次循环")
    for (lab in names(all_results)) {
      if (nrow(all_results[[lab]]) > 0) {
        write.csv(all_results[[lab]],
                  paste0("OB_Model_result_frame_", lab, ".csv"),
                  row.names = FALSE)
      }
    }
  }
  
}  # end cycle loop

# stopCluster(cl)  # 若启用了并行，在此关闭

# ============================================================
# 最终输出
# ============================================================
for (lab in names(all_results)) {
  if (nrow(all_results[[lab]]) > 0) {
    write.csv(all_results[[lab]],
              paste0("OB_Model_result_frame_", lab, ".csv"),
              row.names = FALSE)
    message("已输出: OB_Model_result_frame_", lab, ".csv  (",
            nrow(all_results[[lab]]), " 行)")
  }
}

message("全部完成。")