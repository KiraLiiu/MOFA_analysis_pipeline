# =============================================================================
# Omics Feature Selection & Visualization Utilities
# =============================================================================
# 通用高变异特征选择、cumulative variance (elbow plot)、分布分析等函数
# 可用于多组学数据预处理和可视化
# Author: Multi-omics Analysis Pipeline
# Date: 2024

library(ggplot2)
library(ggpubr)
library(tidyverse)

# 高变异特征选择
select_high_variance_features <- function(data_matrix, variance_threshold = 0.8) {
  feature_variance <- apply(data_matrix, 1, var, na.rm = TRUE)
  sorted_var <- sort(feature_variance, decreasing = TRUE)
  cumulative_var <- cumsum(sorted_var) / sum(sorted_var)
  num_features <- which(cumulative_var >= variance_threshold)[1]
  top_features <- names(sorted_var)[1:num_features]
  selected_features <- data_matrix[top_features, , drop = FALSE]
  return(as.data.frame(selected_features))
}

# Cumulative Variance 曲线（elbow plot）
plot_cumulative_variance <- function(data_matrix, layer_name, save_plot = TRUE, output_dir = NULL) {
  feature_variance <- apply(data_matrix, 1, var, na.rm = TRUE)
  sorted_var <- sort(feature_variance, decreasing = TRUE)
  cumulative_var <- cumsum(sorted_var) / sum(sorted_var)
  plot_data <- data.frame(
    feature_rank = 1:length(cumulative_var),
    cumulative_variance = cumulative_var,
    num_features = 1:length(cumulative_var)
  )
  diff1 <- diff(cumulative_var)
  diff2 <- diff(diff1)
  elbow_idx <- which.max(abs(diff2)) + 2
  p <- ggplot(plot_data, aes(x = feature_rank, y = cumulative_variance)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(data = plot_data[elbow_idx, ], color = "red", size = 3) +
    geom_vline(xintercept = elbow_idx, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_hline(yintercept = cumulative_var[elbow_idx], linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = paste(layer_name, "- Cumulative Variance Curve (Elbow Plot)"),
      x = "Number of Features (Ranked by Variance)",
      y = "Cumulative Variance Explained (%)",
      subtitle = paste("Elbow point at", elbow_idx, "features (", 
                      round(cumulative_var[elbow_idx] * 100, 1), "% variance)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "red"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::comma)
  thresholds <- c(0.5, 0.75, 0.8, 0.9, 0.95)
  threshold_info <- data.frame()
  for (thresh in thresholds) {
    idx <- which(cumulative_var >= thresh)[1]
    if (!is.na(idx)) {
      threshold_info <- rbind(threshold_info, 
                             data.frame(threshold = thresh, 
                                      features = idx, 
                                      variance = cumulative_var[idx]))
    }
  }
  if (nrow(threshold_info) > 0) {
    p <- p + geom_point(data = threshold_info, 
                        aes(x = features, y = variance), 
                        color = "darkgreen", size = 2, alpha = 0.7) +
      geom_text(data = threshold_info, 
                aes(x = features, y = variance, 
                    label = paste0(threshold * 100, "%")), 
                hjust = -0.2, vjust = 0.5, size = 3, color = "darkgreen")
  }
  cat(sprintf("\n%s Cumulative Variance Analysis:\n", layer_name))
  cat(sprintf("Total features: %d\n", length(cumulative_var)))
  cat(sprintf("Elbow point: %d features (%.1f%% variance)\n", elbow_idx, cumulative_var[elbow_idx] * 100))
  cat("Common thresholds:\n")
  for (i in 1:nrow(threshold_info)) {
    cat(sprintf("  %.0f%% variance: %d features\n", 
                threshold_info$threshold[i] * 100, threshold_info$features[i]))
  }
  if (save_plot && !is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    ggsave(file.path(output_dir, paste0(layer_name, "_cumulative_variance.pdf")), 
           plot = p, width = 10, height = 6)
    cat(sprintf("Plot saved: %s\n", file.path(output_dir, paste0(layer_name, "_cumulative_variance.pdf"))))
  }
  return(list(plot = p, elbow_idx = elbow_idx, cumulative_var = cumulative_var))
}

# 分布分析（均值/标准差）
analyze_distributions <- function(data, layer_name, dist_dir = NULL) {
  feature_stats <- data.frame(
    feature = rownames(data),
    mean = apply(data, 1, mean, na.rm = TRUE),
    sd = apply(data, 1, sd, na.rm = TRUE)
  )
  p1 <- ggplot(feature_stats, aes(x = mean)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    labs(title = paste(layer_name, "- Mean Distribution (Post-HVF)"), x = "Mean", y = "Count") +
    theme_minimal()
  p2 <- ggplot(feature_stats, aes(x = sd)) +
    geom_histogram(bins = 50, fill = "darkred", alpha = 0.7) +
    labs(title = paste(layer_name, "- Standard Deviation Distribution (Post-HVF)"), x = "Standard Deviation", y = "Count") +
    theme_minimal()
  p3 <- ggplot(feature_stats, aes(x = mean, y = sd)) +
    geom_point(alpha = 0.6, color = "darkgreen") +
    labs(title = paste(layer_name, "- Mean vs Standard Deviation (Post-HVF)"), x = "Mean", y = "Standard Deviation") +
    theme_minimal()
  combined_plot <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
  if (!is.null(dist_dir)) {
    if (!dir.exists(dist_dir)) {
      dir.create(dist_dir, recursive = TRUE)
    }
    ggsave(file.path(dist_dir, paste0(layer_name, "_distributions_post_HVF.pdf")), 
           plot = combined_plot, width = 15, height = 5)
  }
  cat(sprintf("\n%s Summary Statistics (Post-HVF):\n", layer_name))
  cat(sprintf("Mean - Min: %.3f, Max: %.3f, Median: %.3f\n", 
              min(feature_stats$mean, na.rm = TRUE), 
              max(feature_stats$mean, na.rm = TRUE), 
              median(feature_stats$mean, na.rm = TRUE)))
  cat(sprintf("SD - Min: %.3f, Max: %.3f, Median: %.3f\n", 
              min(feature_stats$sd, na.rm = TRUE), 
              max(feature_stats$sd, na.rm = TRUE), 
              median(feature_stats$sd, na.rm = TRUE)))
  return(feature_stats)
} 