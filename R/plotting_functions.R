
#' boxplot: Plot Response by Group with Statistical Test and Consistent Coloring
#'
#' This function creates a boxplot with jitter and compares response values
#' across groups, automatically choosing between t-test and Wilcoxon test,
#' and ensuring that the control group is always blue and the other orange.
#'
#' @param df A data.frame or data.table with the data
#' @param x Column name (character) for grouping variable
#' @param y Column name (character) for response variable
#' @param filename File name to save the plot (PNG format)
#' @param title Optional custom plot title
#' @param control_group Optional: specify control group name for color consistency
#' @export
boxplot <- function(df, x, y, filename, title = NULL, control_group = NULL) {
  df <- na.omit(df[, .(Group = get(x), Response = get(y))])
  choose_test_method <- function(df) {
    if (length(unique(df$Group)) != 2) return("wilcox.test")
    group_vals <- split(df$Response, df$Group)
    if (all(sapply(group_vals, length) >= 3)) {
      p_normality <- sapply(group_vals, function(g) shapiro.test(g)$p.value)
      if (all(p_normality > 0.05)) return("t.test")
    }
    return("wilcox.test")
  }
  method <- choose_test_method(df)
  group_levels <- unique(df$Group)
  if (is.null(control_group)) control_group <- group_levels[1]
  other_group <- setdiff(group_levels, control_group)
  color_map <- setNames(c("#1f77b4", "#ff7f0e"), c(control_group, other_group))
  p <- ggplot(df, aes(x = Group, y = Response, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, shape = 21, size = 2, color = "black") +
    stat_compare_means(method = method) +
    labs(
      title = ifelse(is.null(title), paste(y, "by", x), title),
      x = x,
      y = y
    ) +
    scale_fill_manual(values = color_map) +
    theme_bw()
  ggsave(filename, plot = p, width = 5, height = 4, dpi = 600)
}

#' correlation: Custom Correlation Plot with Pearson R and p-value
#'
#' This function generates a scatter plot with linear regression fit line
#' and displays correlation coefficient, p-value, and sample size.
#'
#' @param df A data.frame or data.table
#' @param xvar Column name for x-axis variable
#' @param yvar Column name for y-axis variable
#' @param outfile_prefix Output file prefix for saving PNG and PDF
#' @export
correlation <- function(df, xvar, yvar, outfile_prefix) {
  r_val <- cor(df[[xvar]], df[[yvar]], method = "pearson")
  p_val <- cor.test(df[[xvar]], df[[yvar]], method = "pearson")$p.value
  n_val <- nrow(df)
  p <- ggplot(df, aes_string(x = xvar, y = yvar)) +
    geom_point(alpha = 0.7, color = "black", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid", size = 1) +
    labs(x = xvar, y = yvar) +
    geom_text(label = paste0("n = ", n_val, "
p-value = ", signif(p_val, 3),
                             "
R = ", signif(r_val, 3)),
              x = min(df[[xvar]], na.rm = TRUE),
              y = max(df[[yvar]], na.rm = TRUE),
              hjust = 0, vjust = 1, color = "black", size = 5, inherit.aes = FALSE) +
    theme_bw()
  ggsave(paste0(outfile_prefix, ".png"), plot = p, width = 6, height = 6, dpi = 600)
  ggsave(paste0(outfile_prefix, ".pdf"), plot = p, width = 6, height = 6)
  return(p)
}

#' survival: Kaplan-Meier Survival Plot by Group
#'
#' @param fit A survfit object
#' @param data The data frame used to create the fit
#' @param title Plot title
#' @param legend_title Legend title
#' @param legend_labels Vector of legend labels
#' @param palette Vector of colors
#' @param filename File name (without extension) to save plot
#' @export
survival <- function(fit, data, title = "Kaplan-Meier Curve",
                     legend_title = "Group", legend_labels = c("Group1", "Group2"),
                     palette = c("#FFD700", "grey40"), filename = "KM_plot") {
  p <- ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    title = title,
    legend.title = legend_title,
    legend.labs = legend_labels,
    palette = palette,
    xlab = "Days",
    ylab = "Overall Survival Probability",
    surv.median.line = "hv",
    risk.table.height = 0.25,
    ggtheme = theme_minimal(base_size = 14)
  )
  ggsave(paste0(filename, ".png"), plot = print(p), width = 6, height = 6, dpi = 600)
  ggsave(paste0(filename, ".pdf"), plot = print(p), width = 6, height = 6)
}


#' longboxplot: Boxplot for Multiple Genes/Conditions with Custom Strip and Annotation
#'
#' This function creates a faceted-style boxplot comparing expression across genes or conditions,
#' with shaded background strips and annotated p-values.
#'
#' @param expr_merged Data frame with columns: GeneSymbol, Expression, EGFR_Status
#' @param strip_df Data frame with columns: xmin, xmax (for shaded strips)
#' @param pvals Data frame with columns: GeneSymbol, y, label, color
#' @param title_text Title for the plot
#' @param filename Output file name prefix (without extension)
#' @export
longboxplot <- function(expr_merged, strip_df, pvals, title_text, filename) {
  p <- ggplot(expr_merged, aes(x = GeneSymbol, y = Expression, fill = EGFR_Status)) +
    geom_rect(data = strip_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey95", inherit.aes = FALSE) +
    geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.6,
                 alpha = 0.9, color = "gray30") +
    geom_text(data = pvals, aes(x = GeneSymbol, y = y, label = label, color = color),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
    scale_color_identity() +
    scale_fill_manual(values = c("Wildtype" = "#FFFFFF", "Mutated" = "#F40009")) +
    labs(title = title_text, x = "", y = "Expression Ratio") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 14),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = alpha('white', 0.8), color = "gray80"),
      legend.title = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )

  ggsave(paste0(filename, ".png"), plot = p, width = 7.5, height = 5.5, dpi = 600)
  ggsave(paste0(filename, ".pdf"), plot = p, width = 7.5, height = 5.5)
}
