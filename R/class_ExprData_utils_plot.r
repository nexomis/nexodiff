#' @include nexodiff-package.R
#' @include utils.r
NULL

#' Helper function for ExprData$plot_dist_per_sample
#'
#' @param data data frame with columns: batch, group, sample, value
#' @param design private$design see ExprData
#' @param geoms see ExprData$plot_dist_per_sample
#' @param mean_fun see ExprData$plot_dist_per_sample
#' @param fn_str string representation of y-axis transform
#' @return ggplot2 graph
#' @keywords internal
plot_dist_per_sample_helper <- function(
  data, design, geoms = c("boxplot"), mean_fun = NULL, fn_str
) {

  assert_that(all(geoms %in% c("boxplot", "violin", "histo")))

  y_label <- paste("Expression", fn_str)

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

