#' @include nexodiff-package.R
#' @include utils.r
NULL

#' Helper function for ExprData$plot_sum_per_type_per_sample
#'
#' @param data data frame with columns: batch, group, sample, type, total
#' @param design private$design see ExprData
#' @param same_scale see ExprData$plot_sum_per_type_per_sample
#' @param horizontal_bar see ExprData$plot_sum_per_type_per_sample
#' @param log2_expr see ExprData$plot_sum_per_type_per_sample
#' @param exclude_type see ExprData$plot_sum_per_type_per_sample
#' @return ggplot2 graph
plot_sum_per_type_helper <- function(
  data, design, same_scale = TRUE,
  horizontal_bar = TRUE, log2_expr = FALSE, exclude_type = c()
) {
  batch2label <- design$get_b_labels()
  group2label <- design$get_g_labels()

  data <- dplyr::filter(data, ! .data$type %in% exclude_type)

  y_label <- "Sum of expression"
  if (log2_expr) {
    y_label <- paste(y_label, "(log2-transformed)")
  }
  if (! same_scale) {
    f_scales <- "free"
  } else if (horizontal_bar) {
    f_scales <- "free_y"
  } else {
    f_scales <- "free_x"
  }
  data$batch <- factor(data$batch, levels = names(batch2label))
  data$group <- factor(data$group, levels = names(group2label))
  g <- ggplot2::ggplot(data,
    ggplot2::aes(.data$sample, .data$total, fill = .data$type)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggh4x::facet_nested_wrap(
      formula("~ batch + group"),
      nrow = length(unique(data$batch)),
      scales = f_scales
      ,
      labeller = ggplot2::labeller(
        batch = ggplot2::as_labeller(batch2label),
        group = ggplot2::as_labeller(group2label)
      )
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_scientific(digits = 2)
    ) +
    ggplot2::scale_fill_brewer(palette = "BrBG") +
    THEME_NEXOMIS +
    ggplot2::labs(
      fill = "RNA type",
      x = "Samples",
      y = y_label
    )
  if (horizontal_bar) {
    g <- g + ggplot2::coord_flip()
  }
  g
}

#' Helper function for ExprData$plot_dist_per_sample
#'
#' @param data data frame with columns: batch, group, sample, value
#' @param design private$design see ExprData
#' @param geoms see ExprData$plot_dist_per_sample
#' @param mean_fun see ExprData$plot_dist_per_sample
#' @param log2_expr see ExprData$plot_dist_per_sample
#' @return ggplot2 graph
plot_dist_per_sample_helper <- function(
  data, design, geoms = c("boxplot"), mean_fun = NULL, log2_expr = TRUE
) {

  if (length(setdiff(geoms, c("boxplot", "violin", "histo"))) > 0) {
    logging::logerror("unrecognized value for geoms")
    stop()
  }

  y_label <- "Expression"
  if (log2_expr) {
    y_label <- paste(y_label, "(log2-transformed)")
  }

  batch2label <- design$get_b_labels()
  group2label <- design$get_g_labels()

  data$batch <- factor(data$batch, levels = names(batch2label))
  data$group <- factor(data$group, levels = names(group2label))

  g <- ggplot2::ggplot(data,
    ggplot2::aes(.data$sample, .data$value)
  )
  if ("violin" %in% geoms) {
    g <- g +
      ggplot2::geom_violin(trim = FALSE)
  }
  if ("boxplot" %in% geoms) {
    g <- g +
      ggplot2::geom_boxplot(width = 0.1)
  }
  if ("histo" %in% geoms) {
    g <- g +
      ggmulti::geom_histogram_(bins = 30, alpha = 0.5)
  }
  if (! is.null(mean_fun)) {
    g <- g + ggplot2::stat_summary(
      fun = set_mean_function(mean_fun),
      geom = "point", shape = 21, size = 3,
      color = "black", fill = "red"
    )
  }

  g <- g +
    ggh4x::facet_nested_wrap(
      formula("~ batch + group"),
      nrow = length(unique(data$batch)),
      scales = "free_x",
      labeller = ggplot2::labeller(
        batch = ggplot2::as_labeller(batch2label),
        group = ggplot2::as_labeller(group2label)
      )
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_scientific(digits = 2)
    ) +
    THEME_NEXOMIS +
    ggplot2::labs(
      x = "Samples",
      y = y_label
    )
  g
}

