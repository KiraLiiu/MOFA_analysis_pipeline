# Simple Multi-omics Pathway Enrichment Heatmap (Standalone)
# Separated by pathway source (KEGG/REACTOME) with grey for NA p-values

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

prepare_heatmap_data <- function(data) {
  data$transcriptome_log10p <- ifelse(is.na(data$transcriptome_pval), NA, -log10(data$transcriptome_pval))
  data$proteome_log10p <- ifelse(is.na(data$proteome_pval), NA, -log10(data$proteome_pval))
  data$metabolome_log10p <- ifelse(is.na(data$metabolome_pval), NA, -log10(data$metabolome_pval))
  return(data)
}

create_source_separated_heatmap <- function(data, source_name, 
                                            fontsize = 12, fontsize_row = 12, fontsize_col = 10, 
                                            fontsize_number = 8, fontface = "plain", show_values = TRUE,
                                            show_pvalues = FALSE) {
  source_data <- data %>% filter(source == source_name)
  if (nrow(source_data) == 0) {
    warning(paste("No data found for source:", source_name))
    return(NULL)
  }
  mat <- source_data %>%
    select(pathway, transcriptome_log10p, proteome_log10p, metabolome_log10p) %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  mat_t <- t(mat)
  rownames(mat_t) <- c("Transcriptome", "Proteome", "Metabolome")
  pval_mat_t <- NULL
  if (show_pvalues) {
    pval_mat <- source_data %>%
      select(pathway, transcriptome_pval, proteome_pval, metabolome_pval) %>%
      column_to_rownames("pathway") %>%
      as.matrix()
    pval_mat_t <- t(pval_mat)
    rownames(pval_mat_t) <- c("Transcriptome", "Proteome", "Metabolome")
  }
  can_cluster <- TRUE
  if (ncol(mat_t) < 2) {
    can_cluster <- FALSE
  } else {
    col_issues <- apply(mat_t, 2, function(x) {
      all(is.na(x)) || (length(unique(x[!is.na(x)])) <= 1)
    })
    if (any(col_issues) || any(is.infinite(mat_t), na.rm = TRUE)) {
      can_cluster <- FALSE
    }
  }
  create_custom_colors <- function(mat) {
    non_na_values <- mat[!is.na(mat) & is.finite(mat)]
    if (length(non_na_values) == 0) {
      return(list(breaks = c(0, 1), colors = c("#CCCCCC")))
    }
    min_val <- max(0, min(non_na_values, na.rm = TRUE))
    max_val <- max(non_na_values, na.rm = TRUE)
    
    if (min_val == max_val) {
      breaks <- c(min_val - 0.1, min_val + 0.1)
      colors <- c("white")
    } else {
      # 非线性色阶映射：在低值区间更密集，高值区间更稀疏
      # 这样可以更好地显示大部分数据的差异
      
      # 计算分位数来创建非线性断点
      quantiles <- quantile(non_na_values, probs = seq(0, 1, 0.05), na.rm = TRUE)
      
      # 确保包含最小值和最大值
      breaks <- unique(c(min_val, quantiles, max_val))
      
      # 如果断点太少，添加更多中间点
      if (length(breaks) < 50) {
        # 在低值区间（前50%）添加更多断点
        low_breaks <- seq(min_val, quantile(non_na_values, 0.5, na.rm = TRUE), length.out = 50)
        high_breaks <- seq(quantile(non_na_values, 0.5, na.rm = TRUE), max_val, length.out = 30)
        breaks <- unique(c(low_breaks, high_breaks))
      }
      
      # 确保断点数量合理
      if (length(breaks) > 100) {
        breaks <- breaks[seq(1, length(breaks), length.out = 100)]
      }
      
      # 创建颜色渐变
      colors <- colorRampPalette(c("white", "#F5F0F3", "#E6D5DC", "#D7BAC5", "#C89FAE", "#B98497", "#AA6980", "#9B4E69", "#8C3352", "#8f2953"))(length(breaks) - 1)
    }
    return(list(breaks = breaks, colors = colors))
  }
  color_info <- create_custom_colors(mat_t)
  p <- pheatmap(
    mat_t,
    color = color_info$colors,
    breaks = color_info$breaks,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_cols = can_cluster,
    main = paste(source_name, "Pathways - Multi-omics Enrichment"),
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    angle_col = "315",
    border_color = "white",
    cellwidth = 25,
    cellheight = 25,
    display_numbers = show_values,
    number_format = if(show_pvalues) "%.2e" else "%.1f",
    fontsize_number = fontsize_number,
    cell_fun = if(show_pvalues && show_values) {
      function(j, i, x, y, width, height, fill) {
        if (!is.null(pval_mat_t) && !is.na(pval_mat_t[i, j])) {
          grid.text(sprintf("%.2e", pval_mat_t[i, j]), x, y, 
                    gp = gpar(fontsize = fontsize_number, fontface = fontface))
        }
      }
    } else NULL,
    treeheight_row = 0,
    treeheight_col = if(can_cluster) 50 else 0
  )
  return(p)
}

create_selected_sources_heatmap <- function(data, selected_sources = NULL, show_source_annotation = TRUE, 
                                            fontsize = 11, fontsize_row = 12, fontsize_col = 9, 
                                            fontsize_number = 7, fontface = "plain", show_values = TRUE, 
                                            show_pvalues = FALSE) {
  if (is.null(selected_sources)) {
    selected_sources <- unique(data$source)
  }
  filtered_data <- data %>% filter(source %in% selected_sources)
  if (nrow(filtered_data) == 0) {
    warning("No data found for selected sources:", paste(selected_sources, collapse = ", "))
    return(NULL)
  }
  mat <- filtered_data %>%
    select(pathway, transcriptome_log10p, proteome_log10p, metabolome_log10p) %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  mat_t <- t(mat)
  rownames(mat_t) <- c("Transcriptome", "Proteome", "Metabolome")
  pval_mat <- NULL
  if (show_pvalues) {
    pval_mat <- filtered_data %>%
      select(pathway, transcriptome_pval, proteome_pval, metabolome_pval) %>%
      column_to_rownames("pathway") %>%
      as.matrix()
    pval_mat_t <- t(pval_mat)
    rownames(pval_mat_t) <- c("Transcriptome", "Proteome", "Metabolome")
  }
  col_annotation <- NULL
  ann_colors <- NULL
  if (show_source_annotation && length(selected_sources) > 1) {
    col_annotation <- data.frame(
      Source = filtered_data$source,
      row.names = filtered_data$pathway
    )
    n_sources <- length(selected_sources)
    if (n_sources == 1) {
      source_colors <- c("#8B5CF6")
    } else if (n_sources == 2) {
      source_colors <- c("#8B5CF6", "#EC4899")
    } else {
      source_colors <- c("#8B5CF6", "#EC4899", "#10B981", "#F59E0B", "#EF4444", "#3B82F6")
      source_colors <- source_colors[1:n_sources]
    }
    names(source_colors) <- selected_sources
    ann_colors <- list(Source = source_colors)
  }
  can_cluster <- TRUE
  if (ncol(mat_t) < 2) {
    can_cluster <- FALSE
    cat("Warning: Less than 2 pathways, disabling clustering\n")
  } else {
    col_issues <- apply(mat_t, 2, function(x) {
      all(is.na(x)) || (length(unique(x[!is.na(x)])) <= 1)
    })
    if (any(col_issues)) {
      cat("Warning: Some pathways have all NA or identical values, disabling clustering\n")
      can_cluster <- FALSE
    }
    if (any(is.infinite(mat_t), na.rm = TRUE)) {
      cat("Warning: Infinite values detected, replacing with max finite value\n")
      finite_vals <- mat_t[is.finite(mat_t)]
      if (length(finite_vals) > 0) {
        max_finite <- max(finite_vals, na.rm = TRUE)
        mat_t[is.infinite(mat_t)] <- max_finite
      }
    }
  }
  create_custom_colors <- function(mat) {
    non_na_values <- mat[!is.na(mat) & is.finite(mat)]
    if (length(non_na_values) == 0) {
      return(list(breaks = c(0, 1), colors = c("#CCCCCC")))
    }
    min_val <- max(0, min(non_na_values, na.rm = TRUE))
    max_val <- max(non_na_values, na.rm = TRUE)
    
    if (min_val == max_val) {
      breaks <- c(min_val - 0.1, min_val + 0.1)
      colors <- c("white")
    } else {
      # 非线性色阶映射：在低值区间更密集，高值区间更稀疏
      # 这样可以更好地显示大部分数据的差异
      
      # 计算分位数来创建非线性断点
      quantiles <- quantile(non_na_values, probs = seq(0, 1, 0.05), na.rm = TRUE)
      
      # 确保包含最小值和最大值
      breaks <- unique(c(min_val, quantiles, max_val))
      
      # 如果断点太少，添加更多中间点
      if (length(breaks) < 50) {
        # 在低值区间（前50%）添加更多断点
        low_breaks <- seq(min_val, quantile(non_na_values, 0.5, na.rm = TRUE), length.out = 50)
        high_breaks <- seq(quantile(non_na_values, 0.5, na.rm = TRUE), max_val, length.out = 30)
        breaks <- unique(c(low_breaks, high_breaks))
      }
      
      # 确保断点数量合理
      if (length(breaks) > 100) {
        breaks <- breaks[seq(1, length(breaks), length.out = 100)]
      }
      
      # 创建颜色渐变
      colors <- colorRampPalette(c("white", "#F5F0F3", "#E6D5DC", "#D7BAC5", "#C89FAE", "#B98497", "#AA6980", "#9B4E69", "#8C3352", "#8f2953"))(length(breaks) - 1)
    }
    return(list(breaks = breaks, colors = colors))
  }
  color_info <- create_custom_colors(mat_t)
  title_sources <- paste(selected_sources, collapse = " & ")
  p <- pheatmap(
    mat_t,
    color = color_info$colors,
    breaks = color_info$breaks,
    na_col = "white",
    annotation_col = col_annotation,
    annotation_colors = ann_colors,
    cluster_rows = FALSE,
    cluster_cols = can_cluster,
    clustering_distance_cols = if(can_cluster) "euclidean" else NULL,
    clustering_method = if(can_cluster) "complete" else NULL,
    main = paste(title_sources, "Pathways - Multi-omics Enrichment"),
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    angle_col = "315",
    border_color = "white",
    cellwidth = 20,
    cellheight = 25,
    display_numbers = show_values,
    number_format = if(show_pvalues) "%.2e" else "%.1f",
    fontsize_number = fontsize_number,
    cell_fun = if(show_pvalues && show_values) {
      function(j, i, x, y, width, height, fill) {
        if (!is.null(pval_mat_t) && !is.na(pval_mat_t[i, j])) {
          grid.text(sprintf("%.2e", pval_mat_t[i, j]), x, y, 
                    gp = gpar(fontsize = fontsize_number, fontface = fontface))
        }
      }
    } else NULL,
    treeheight_row = 0,
    treeheight_col = if(can_cluster) 50 else 0
  )
  return(p)
}

generate_pathway_heatmaps <- function(data, 
                                      selected_sources = NULL, 
                                      show_source_annotation = TRUE,
                                      fontsize = 11,
                                      fontsize_row = 12,
                                      fontsize_col = 9,
                                      fontsize_number = 7,
                                      fontface = "plain",
                                      show_values = TRUE,
                                      show_pvalues = FALSE,
                                      save_output = FALSE,
                                      output_format = "pdf",
                                      plot_width = 10,
                                      plot_height = 6,
                                      output_dir = "pathway_heatmaps") {
  data_processed <- prepare_heatmap_data(data)
  available_sources <- unique(data_processed$source)
  cat("Available sources:", paste(available_sources, collapse = ", "), "\n")
  if (is.null(selected_sources)) {
    selected_sources <- available_sources
    cat("No sources specified, using all available sources\n")
  } else {
    missing_sources <- setdiff(selected_sources, available_sources)
    if (length(missing_sources) > 0) {
      warning("Selected sources not found in data: ", paste(missing_sources, collapse = ", "))
    }
    selected_sources <- intersect(selected_sources, available_sources)
    cat("Using selected sources:", paste(selected_sources, collapse = ", "), "\n")
  }
  cat("Creating heatmap for selected sources...\n")
  heatmap_plot <- create_selected_sources_heatmap(data_processed, selected_sources, show_source_annotation,
                                                  fontsize, fontsize_row, fontsize_col, fontsize_number, fontface,
                                                  show_values, show_pvalues)
  if (is.null(heatmap_plot)) {
    cat("No heatmap generated - check your source selection\n")
    return(NULL)
  }
  if (save_output) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    sources_name <- paste(tolower(selected_sources), collapse = "_")
    filename_base <- file.path(output_dir, paste0(sources_name, "_pathways_heatmap"))
    if (output_format == "png") {
      filename <- paste0(filename_base, ".png")
      png(filename, width = plot_width * 100, height = plot_height * 100, res = 300)
    } else if (output_format == "svg") {
      filename <- paste0(filename_base, ".svg")
      svg(filename, width = plot_width, height = plot_height)
    } else {
      filename <- paste0(filename_base, ".pdf")
      pdf(filename, width = plot_width, height = plot_height)
    }
    grid::grid.newpage()
    grid::grid.draw(heatmap_plot$gtable)
    dev.off()
    cat(paste("Heatmap saved:", filename, "\n"))
    if (length(selected_sources) > 1) {
      for (source in selected_sources) {
        individual_plot <- create_source_separated_heatmap(data_processed, source,
                                                           fontsize, fontsize_row, fontsize_col, fontsize_number, fontface,
                                                           show_values, show_pvalues)
        if (!is.null(individual_plot)) {
          individual_filename_base <- file.path(output_dir, paste0(tolower(source), "_only_heatmap"))
          if (output_format == "png") {
            individual_filename <- paste0(individual_filename_base, ".png")
            png(individual_filename, width = plot_width * 100, height = plot_height * 100, res = 300)
          } else if (output_format == "svg") {
            individual_filename <- paste0(individual_filename_base, ".svg")
            svg(individual_filename, width = plot_width, height = plot_height)
          } else {
            individual_filename <- paste0(individual_filename_base, ".pdf")
            pdf(individual_filename, width = plot_width, height = plot_height)
          }
          grid::grid.newpage()
          grid::grid.draw(individual_plot$gtable)
          dev.off()
          cat(paste("Individual heatmap saved:", individual_filename, "\n"))
        }
      }
    }
  }
  cat("Heatmap generation completed!\n")
  return(heatmap_plot)
} 