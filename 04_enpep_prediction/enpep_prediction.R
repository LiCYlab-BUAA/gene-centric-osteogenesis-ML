# ============================================================
# 调取各模型最优版本，对 ENPEP 进行预测并可视化
# 依赖：
#   - log2_dataset_frame.csv          原始表达数据
#   - OB_Model_result_frame_*.csv     训练结果（含最优 Cycle）
#   - cycles/                         保存的 RDS 模型文件
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(caret)

# ============================================================
# 1. 读取原始数据，重建 ENPEP_data
# ============================================================
exp_change <- read.csv("log2_dataset_frame.csv", header = TRUE,
                       stringsAsFactors = FALSE)
colnames(exp_change)[1] <- "Symbol"

ENPEP_data <- exp_change[which(exp_change$Symbol == "ENPEP"), ]
if (nrow(ENPEP_data) == 0) stop("在 log2_dataset_frame.csv 中未找到 ENPEP 基因行。")
rownames(ENPEP_data) <- "ENPEP"
ENPEP_data <- ENPEP_data[, -1]                            # 去掉 Symbol 列
if ("Label" %in% colnames(ENPEP_data)) {
  ENPEP_data <- ENPEP_data[, -which(colnames(ENPEP_data) == "Label")]
}
message("ENPEP 特征维度：", ncol(ENPEP_data), " 列")

# ============================================================
# 2. 读取所有模型结果 CSV，挑出每个模型 AUC 最优的那行
# ============================================================
csv_files <- list.files(pattern = "^OB_Model_result_frame_.*\\.csv$",
                        full.names = TRUE)
if (length(csv_files) == 0) stop("未找到 OB_Model_result_frame_*.csv，请确认工作目录。")

all_df <- lapply(csv_files, function(f) {
  tryCatch(read.csv(f, header = TRUE, stringsAsFactors = FALSE),
           error = function(e) NULL)
})
all_df    <- Filter(Negate(is.null), all_df)
combined  <- bind_rows(all_df)
combined$AUC      <- suppressWarnings(as.numeric(combined$AUC))
combined$Accuracy <- suppressWarnings(as.numeric(combined$Accuracy))
combined$Cycle    <- suppressWarnings(as.integer(combined$Cycle))

best_rows <- combined %>%
  filter(!is.na(AUC), !is.na(Cycle)) %>%
  group_by(Model) %>%
  arrange(desc(AUC), desc(Accuracy)) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(desc(AUC))

message("\n共找到 ", nrow(best_rows), " 个模型的最优行")

# ============================================================
# 3. 逐模型：加载 RDS → predict → 记录结果
# ============================================================
results <- list()

for (i in seq_len(nrow(best_rows))) {
  
  model_name <- best_rows$Model[i]
  cycle_id   <- best_rows$Cycle[i]
  rds_path   <- file.path("cycles",
                          sprintf("%s_cycle%05d.rds", model_name, cycle_id))
  
  # —— 3a. 从 CSV 直接读取已记录的 ENPEP 预测值（快速路径）——
  enpep_csv <- suppressWarnings(as.numeric(best_rows$ENPEP_Result[i]))
  
  # —— 3b. 若 RDS 存在，重新预测（精确路径）——
  enpep_rds_prob <- NA_real_
  enpep_rds_class <- NA_character_
  rds_loaded <- FALSE
  
  if (file.exists(rds_path)) {
    fitted <- tryCatch(readRDS(rds_path), error = function(e) NULL)
    if (!is.null(fitted)) {
      # 概率预测
      prob_out <- tryCatch({
        p <- predict(fitted, newdata = ENPEP_data, type = "prob")
        p[1, "OB"]
      }, error = function(e) NA_real_)
      
      # 类别预测
      class_out <- tryCatch({
        as.character(predict(fitted, newdata = ENPEP_data, type = "raw"))
      }, error = function(e) NA_character_)
      
      enpep_rds_prob  <- prob_out
      enpep_rds_class <- class_out
      rds_loaded <- TRUE
    }
  } else {
    message("  [跳过 RDS] 文件不存在：", rds_path)
  }
  
  # —— 3c. 合并：优先用 RDS 重新预测的概率，否则用 CSV 记录值 ——
  final_prob <- if (!is.na(enpep_rds_prob)) enpep_rds_prob else enpep_csv
  
  # 类别判断（threshold = 0.5）
  final_class <- if (!is.na(final_prob)) {
    ifelse(final_prob >= 0.5, "OB (Osteogenesis-related)", "N (Negative)")
  } else if (!is.na(enpep_rds_class)) {
    enpep_rds_class
  } else {
    NA_character_
  }
  
  results[[i]] <- data.frame(
    Model         = model_name,
    Cycle         = cycle_id,
    Best_AUC      = round(best_rows$AUC[i], 4),
    Best_Accuracy = round(best_rows$Accuracy[i], 4),
    RDS_loaded    = rds_loaded,
    ENPEP_Prob_OB = round(final_prob, 4),
    ENPEP_Class   = final_class,
    Source        = ifelse(rds_loaded, "RDS重新预测", "CSV记录值"),
    stringsAsFactors = FALSE
  )
  
  message(sprintf("  [%2d/%d] %-20s  Prob=%.4f  Class=%s  (%s)",
                  i, nrow(best_rows), model_name,
                  ifelse(is.na(final_prob), -1, final_prob),
                  ifelse(is.na(final_class), "NA", final_class),
                  ifelse(rds_loaded, "RDS", "CSV")))
}

enpep_df_raw <- bind_rows(results) %>%
  filter(!is.na(ENPEP_Prob_OB)) %>%
  arrange(desc(ENPEP_Prob_OB)) %>%
  mutate(
    Prediction = factor(
      ifelse(ENPEP_Prob_OB >= 0.5, "Predicted OB (>=0.5)", "Predicted N (<0.5)"),
      levels = c("Predicted OB (>=0.5)", "Predicted N (<0.5)")
    ),
    # Source 统一简写（柱内标注）
    Source_short = ifelse(Source == "RDS重新预测", "RDS", "CSV"),
    # factor 排序：概率降序 → 图上从上到下由高到低
    Model_f = factor(Model, levels = rev(Model))
  )

message("\n========== ENPEP 预测结果汇总 ==========")
print(as.data.frame(enpep_df_raw %>% select(-Model_f)), row.names = FALSE)


write.csv(enpep_df_raw %>% select(-Model_f),
          "ENPEP_prediction_results.csv", row.names = FALSE)
message("\n已保存：ENPEP_prediction_results.csv")

# ============================================================
# 4. 绘图
# ============================================================
pred_colors <- c("Predicted OB (>=0.5)" = "#C0392B",
                 "Predicted N (<0.5)"   = "#2980B9")

# ---- 图1：ENPEP 预测概率横向柱状图 -------------------------
p1 <- ggplot(enpep_df_raw,
             aes(x = Model_f, y = ENPEP_Prob_OB, fill = Prediction)) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey40", linewidth = 0.7) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.3f", ENPEP_Prob_OB)),
            hjust = -0.1, size = 3, color = "grey15") +
  geom_text(aes(y = 0.02, label = Source_short),
            hjust = 0, size = 2.3, color = "white", fontface = "bold") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = pred_colors, name = "ENPEP Prediction") +
  scale_y_continuous(limits  = c(0, 1.15),
                     breaks  = seq(0, 1, 0.1),
                     labels  = label_number(accuracy = 0.1),
                     expand  = expansion(mult = c(0, 0))) +
  annotate("text", x = -Inf, y = 0.5,
           label = "threshold = 0.5", vjust = -0.5, hjust = 0.5,
           size = 3, color = "grey40", fontface = "italic") +
  labs(
    title    = "Predicted OB Probability of ENPEP by Each Optimal Model",
    subtitle = "White text in bars: prediction source (RDS = re-inferred / CSV = recorded during training) | Sorted by probability (high to low)",
    x        = NULL,
    y        = "P(ENPEP = OB)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(color = "grey40", size = 9),
    legend.position    = c(0.80, 0.08),
    legend.background  = element_rect(fill = alpha("white", 0.85),
                                      color = "grey80"),
    legend.key.size    = unit(0.45, "cm"),
    legend.text        = element_text(size = 8.5),
    legend.title       = element_text(size = 9, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 9),
    plot.margin        = margin(10, 35, 10, 10)
  )

n_m  <- nrow(enpep_df_raw)
figh <- max(5, n_m * 0.40 + 2)
ggsave("ENPEP_prediction_barplot.pdf",
       p1, width = 10, height = figh, device = cairo_pdf)
ggsave("ENPEP_prediction_barplot.png",
       p1, width = 10, height = figh, dpi = 600, bg = "white")
message("已保存：ENPEP_prediction_barplot.pdf / .png")

# ---- 图2：预测概率 vs 模型 AUC 散点图 ----------------------
p2 <- ggplot(enpep_df_raw,
             aes(x = Best_AUC, y = ENPEP_Prob_OB,
                 color = Prediction, label = Model)) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  geom_vline(xintercept = 0.8, linetype = "dotted",
             color = "grey60", linewidth = 0.5) +
  geom_point(size = 3.5, alpha = 0.85) +
  {if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 30,
                             segment.color = "grey70", segment.size = 0.3)
    else
      geom_text(size = 2.5, vjust = -0.8, hjust = 0.5, color = "grey20")} +
  scale_color_manual(values = pred_colors, name = "ENPEP Classification") +
  scale_x_continuous(breaks = seq(0.5, 1, 0.05),
                     labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     labels = label_number(accuracy = 0.1)) +
  labs(
    title    = "Model AUC vs. ENPEP Predicted Probability",
    subtitle = "Horizontal line: decision threshold 0.5   Vertical line: AUC reference line 0.80",
    x        = "Optimal Model AUC",
    y        = "P(ENPEP = OB)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 9),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("ENPEP_AUC_scatter.pdf",
       p2, width = 9, height = 7, device = cairo_pdf)
ggsave("ENPEP_AUC_scatter.png",
       p2, width = 9, height = 7, dpi = 600, bg = "white")
message("已保存：ENPEP_AUC_scatter.pdf / .png")

p2.1 <- ggplot(enpep_df_raw,
               aes(x = Best_AUC, y = ENPEP_Prob_OB,
                   color = Prediction, label = Model)) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  geom_vline(xintercept = 0.8, linetype = "dotted",
             color = "grey60", linewidth = 0.5) +
  geom_point(size = 3.5, alpha = 0.85) +
  {if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_text_repel(size = 3, max.overlaps = 30,
                             segment.color = "grey70", segment.size = 0.3)
    else
      geom_text(size = 3, vjust = -0.8, hjust = 0.5, color = "grey20")} +
  scale_color_manual(values = pred_colors, name = "ENPEP Classification") +
  scale_x_continuous(breaks = seq(0.5, 1, 0.05),
                     labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     labels = label_number(accuracy = 0.1)) +
  labs(
    # title    = "Model AUC vs. ENPEP Predicted Probability",
    # subtitle = "Horizontal line: decision threshold 0.5   Vertical line: AUC reference line 0.80",
    x        = "Optimal Model AUC",
    y        = "P(ENPEP = OB)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "grey40"),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )
ggsave("ENPEP_AUC_scatter1.pdf",
       p2.1, width = 4, height = 4, device = cairo_pdf)
ggsave("ENPEP_AUC_scatter1.png",
       p2.1, width = 4, height = 4, dpi = 600, bg = "white")

# ---- 图3：预测一致性统计饼图 --------------------------------
vote_tbl <- enpep_df_raw %>%
  count(Prediction) %>%
  mutate(
    pct   = n / sum(n),
    label = sprintf("%s\n%d models (%.0f%%)",
                    as.character(Prediction), n, pct * 100)
  )

p3 <- ggplot(vote_tbl, aes(x = "", y = pct, fill = Prediction)) +
  geom_col(width = 1, color = "white", linewidth = 1) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4, color = "white", fontface = "bold",
            lineheight = 1.4) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = pred_colors, guide = "none") +
  labs(
    title    = "Directional Consistency of ENPEP Predictions",
    subtitle = sprintf("Total %d optimal models voted", sum(vote_tbl$n))
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 10, hjust = 0.5)
  )

ggsave("ENPEP_vote_pie.pdf",
       p3, width = 6, height = 6, device = cairo_pdf)
ggsave("ENPEP_vote_pie.png",
       p3, width = 6, height = 6, dpi = 600, bg = "white")
message("已保存：ENPEP_vote_pie.pdf / .png")

