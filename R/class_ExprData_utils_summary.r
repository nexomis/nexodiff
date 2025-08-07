#' @include nexodiff-package.R
#' @include utils.r
NULL

#' Helper function for ExprData$show_etags_summary
#'
#' @param annotation private$annotation (see ExprData)
#' @param at_gene_level private$at_gene_level (see ExprData)
#' @param data input data
#' @param tr_fn transformation function to apply to data
#' @param sum_fn aggregation function to apply after transformation
#' @return A list of count tables per feature
summarize_etags <- function(
  annotation, at_gene_level, data, tr_fn, sum_fn
) {
  # Apply transformation function to data
  assert_that(is.matrix(data))
  assert_that(is.function(sum_fn))
  if (!is.null(tr_fn)) {
    data <- tr_fn(data)
  }

  results <- list()
  etag_id <- if (at_gene_level) "tgid" else "txid"
  data <- as.data.frame(data)
  for (metric_id in c("tax_id", "tax_name", "type")) {
    data[, metric_id] <-
      annotation$generate_translate_dict(
        etag_id, metric_id
      )[row.names(data)]
    results[[metric_id]] <-
      aggregate(formula(paste(". ~", metric_id)), data, sum_fn)
    row.names(results[[metric_id]]) <- results[[metric_id]][, metric_id]
    results[[metric_id]][, metric_id] <- NULL
    data[, metric_id] <- NULL
  }
  results
}

#' Helper function for ExprData$sum_per_type_per_sample
#'
#' @param results see summarize_etags
#' @param design private$design
#' @return list of plots
plot_summarized_etags <- function(results, design) {
  data_design <-
    design$get_pairwise_design()[, c("batch", "group", "sample")]
  batch2label <- design$get_b_labels()
  group2label <- design$get_g_labels()
  tlabels <- c("Taxonomy ID", "Taxonomy", "RNA Type")
  names(tlabels) <- c("tax_id", "tax_name", "type")
  purrr::map(
    names(results),
    function(in_type) {
      data <- dplyr::mutate(
        results[[in_type]], type = row.names(results[[in_type]])
      ) %>%
        tidyr::pivot_longer(
          tidyselect::all_of(names(results[[in_type]])), names_to = "sample"
        ) %>%
        dplyr::right_join(data_design, by = "sample")
      ggplot2::ggplot(data,
        ggplot2::aes(.data$sample, .data$value, fill = .data$type)
      ) +
        ggplot2::geom_bar(stat = "identity") +
        ggh4x::facet_nested_wrap(
          formula("~ batch + group"),
          nrow = length(unique(batch2label)),
          scales = "free_x"
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
          fill = tlabels[in_type],
          x = "Samples",
          y = "Expression Values"
        )
    }
  ) %>%
    purrr::set_names(names(results))
}
