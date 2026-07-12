# ============================================================
# 从所有模型结果 CSV 中筛选最优模型，并用 ggplot2 绘制效能比较图
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(scales)

# ============================================================
# 1. 读取所有模型结果 CSV
# ============================================================
csv_files <- list.files(pattern = "^OB_Model_result_frame_.*\\.csv$",
                        full.names = TRUE)

if (length(csv_files) == 0) {
  stop("未找到 OB_Model_result_frame_*.csv 文件，请确认工作目录正确。")
}

message("找到 ", length(csv_files), " 个结果文件：")
message(paste(" -", csv_files, collapse = "\n"))

all_df <- lapply(csv_files, function(f) {
  df <- tryCatch(read.csv(f, header = TRUE, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
})
all_df <- Filter(Negate(is.null), all_df)
combined <- bind_rows(all_df)

message("\n合并后总行数：", nrow(combined))

# ============================================================
# 2. 数值化效能列
# ============================================================
numeric_cols <- c("AUC", "Accuracy", "Sensitivity", "Specificity",
                  "FPR", "FNR", "Precision", "F1", "MCC", "Kappa",
                  "Cor", "MAE")

for (col in numeric_cols) {
  if (col %in% colnames(combined)) {
    combined[[col]] <- suppressWarnings(as.numeric(combined[[col]]))
  }
}

# ============================================================
# 3. 每个模型选出最优一行
#    优先级：AUC 降序 → Accuracy 降序
# ============================================================
best_models <- combined %>%
  filter(!is.na(AUC), !is.na(Accuracy)) %>%
  group_by(Model) %>%
  arrange(desc(AUC), desc(Accuracy)) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(desc(AUC))                  # 最终按 AUC 从高到低排列

message("\n最优模型汇总（AUC 降序）：")
print(best_models %>% select(Model, Cycle, AUC, Accuracy, Sensitivity,
                             Specificity, F1, MCC, Kappa) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))),
      n = Inf)

# 保存最优模型汇总表
write.csv(best_models, "Best_model_per_algorithm.csv", row.names = FALSE)
message("\n最优模型汇总已保存：Best_model_per_algorithm.csv")

# ============================================================
# 4. 构造绘图长格式数据
#    展示指标：AUC / Accuracy / Sensitivity / Specificity /
#              Precision / F1 / MCC / Kappa
# ============================================================
plot_metrics <- c("AUC", "Accuracy", "Sensitivity", "Specificity",
                  "Precision", "F1", "MCC", "Kappa")

# 仅保留实际存在的列
plot_metrics <- intersect(plot_metrics, colnames(best_models))

plot_df <- best_models %>%
  select(Model, all_of(plot_metrics)) %>%
  pivot_longer(cols = all_of(plot_metrics),
               names_to  = "Metric",
               values_to = "Value") %>%
  mutate(
    # 固定模型顺序（按 AUC 降序）
    Model  = factor(Model, levels = best_models$Model),
    Metric = factor(Metric, levels = plot_metrics)
  )

# ============================================================
# 5. 配色方案（每个指标一个色系）
# ============================================================
metric_colors <- c(
  AUC         = "#2166AC",
  Accuracy    = "#4DAF4A",
  Sensitivity = "#FF7F00",
  Specificity = "#984EA3",
  Precision   = "#E41A1C",
  F1          = "#A65628",
  MCC         = "#F781BF",
  Kappa       = "#999999"
)[plot_metrics]

# ============================================================
# 6. 图1：分面柱状图 —— 每个指标一个面板，X 轴为模型
# ============================================================
p1 <- ggplot(plot_df, aes(x = Model, y = Value, fill = Metric)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(!is.na(Value),
                               sprintf("%.3f", Value), "NA")),
            hjust = -0.08, size = 2.5, color = "grey20") +
  coord_flip(clip = "off") +
  facet_wrap(~ Metric, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = metric_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x        = NULL,
    y        = "Performance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle    = element_text(color = "grey40", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    axis.text.y      = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin      = margin(10, 20, 10, 10)
  )

# 动态调整图高（模型越多越高）
n_models  <- nrow(best_models)
fig_h1    <- max(8, n_models * 0.35 * ceiling(length(plot_metrics) / 2))
fig_w1    <- 14

ggsave("Model_performance_facet.pdf",
       plot = p1, width = fig_w1, height = fig_h1,
       device = cairo_pdf)
ggsave("Model_performance_facet.png",
       plot = p1, width = fig_w1, height = fig_h1,
       dpi = 600, bg = "white")
message("已保存：Model_performance_facet.pdf / .png")

# ============================================================
# 7. 图2：热图式综合比较 —— X 模型 × Y 指标，填色为值
# ============================================================
p2 <- ggplot(plot_df, aes(x = Metric, y = Model, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(!is.na(Value),
                               sprintf("%.3f", Value), "—")),
            size = 2.8, color = "black") +
  # scale_fill_gradient2(
  #   low      = "#D73027",
  #   mid      = "#FFFFBF",
  #   high     = "#1A9850",
  #   midpoint = 0.5,
  #   na.value = "grey85",
  #   name     = "Performance",
  #   limits   = c(0, 1),
  #   oob      = squish
  # ) +
  scale_fill_gradient2(
    low      = "white",
    high     = "darkblue",
    midpoint = 0.5,
    na.value = "grey85",
    name     = "Performance",
    limits   = c(0, 1),
    oob      = squish
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(plot_df$Model))) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x        = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(color = "grey40", size = 9),
    axis.text.x     = element_text(face = "bold", size = 10, angle = 30,
                                   hjust = 1),
    axis.text.y     = element_text(size = 9),
    panel.grid      = element_blank(),
    legend.position = "right",
    plot.margin     = margin(10, 10, 10, 10)
  )

fig_h2 <- max(5, n_models * 0.5 + 2)
fig_w2 <- 10

ggsave("Model_performance_heatmap.pdf",
       plot = p2, width = fig_w2, height = fig_h2,
       device = cairo_pdf)
ggsave("Model_performance_heatmap.png",
       plot = p2, width = fig_w2, height = fig_h2,
       dpi = 600, bg = "white")
message("已保存：Model_performance_heatmap.pdf / .png")

# ============================================================
# 8. 图3：AUC 单指标排名柱状图（突出主效能）
# ============================================================
auc_df <- best_models %>%
  mutate(
    Model     = factor(Model, levels = rev(best_models$Model)),
    auc_label = sprintf("%.4f", AUC),
    # 按 AUC 分三档着色
    tier      = case_when(
      AUC >= 0.9              ~ "Rank1 (≥0.90)",
      AUC >= 0.8 & AUC < 0.9 ~ "Rank2 (0.80–0.90)",
      TRUE                    ~ "Rank3 (<0.80)"
    ),
    tier = factor(tier, levels = c("Rank1 (≥0.90)",
                                   "Rank2 (0.80–0.90)",
                                   "Rank3 (<0.80)"))
  )

tier_colors <- c(
  "Rank1 (≥0.90)"   = "#1A9850",
  "Rank2 (0.80–0.90)" = "#74ADD1",
  "Rank3 (<0.80)"    = "#F46D43"
)

p3 <- ggplot(auc_df, aes(x = Model, y = AUC, fill = tier)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  geom_hline(yintercept = c(0.8, 0.9), linetype = "dashed",
             color = c("#74ADD1", "#1A9850"), linewidth = 0.6) +
  geom_text(aes(label = auc_label),
            hjust = -0.1, size = 3, color = "grey20") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = tier_colors, name = "AUC") +
  scale_y_continuous(
    limits = c(0, 1.12),
    breaks = seq(0, 1, 0.1),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  annotate("text", x = Inf, y = 0.9, label = "0.90",
           hjust = 1.15, vjust = -0.4, size = 3,
           color = "#1A9850", fontface = "italic") +
  annotate("text", x = Inf, y = 0.8, label = "0.80",
           hjust = 1.15, vjust = -0.4, size = 3,
           color = "#74ADD1", fontface = "italic") +
  labs(
    title    = NULL,
    subtitle = NULL,
    x        = NULL,
    y        = "AUC"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 9),
    legend.position  = c(0.82, 0.08),
    legend.background = element_rect(fill = alpha("white", 0.8),
                                     color = "grey80"),
    legend.key.size  = unit(0.45, "cm"),
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 9),
    plot.margin        = margin(10, 30, 10, 10)
  )

fig_h3 <- max(5, n_models * 0.38 + 1.5)

ggsave("Model_AUC_ranking.pdf",
       plot = p3, width = 9, height = fig_h3,
       device = cairo_pdf)
ggsave("Model_AUC_ranking.png",
       plot = p3, width = 9, height = fig_h3,
       dpi = 600, bg = "white")
message("已保存：Model_AUC_ranking.pdf / .png")

# ============================================================
# 输出最终 Cycle 编号（用于定位对应 cycles/ 下的 RDS 文件）
# ============================================================
message("\n========== 最优模型对应的 RDS 文件（cycles/） ==========")
rds_paths <- best_models %>%
  transmute(
    Model    = Model,
    Cycle    = Cycle,
    AUC      = round(AUC, 4),
    RDS_file = sprintf("cycles/%s_cycle%05d.rds", Model, as.integer(Cycle))
  )
print(rds_paths, n = Inf)

write.csv(rds_paths, "Best_model_RDS_paths.csv", row.names = FALSE)
message("RDS 路径列表已保存：Best_model_RDS_paths.csv")

message("\n全部完成。共输出：")
message("  Best_model_per_algorithm.csv  —— 最优行的完整效能数据")
message("  Best_model_RDS_paths.csv      —— 对应 RDS 文件路径")
message("  Model_performance_facet.pdf/png   —— 分面柱状图")
message("  Model_performance_heatmap.pdf/png —— 热图")
message("  Model_AUC_ranking.pdf/png         —— AUC 排名图")