#' @include nexodiff-package.R
#' @include utils.r
NULL

#' Helper function for ExprData$show_etags_summary
#'
#' @param selected_ids private$selected_ids (see ExprData)
#' @param annotation private$annotation (see ExprData)
#' @param at_gene_level private$at_gene_level (see ExprData)
#' @param type see ExprData$show_etags_summary
#' @param data input data
#' @return A list of count tables per feature
summarize_etags <- function(
  selected_ids, annotation, at_gene_level, type = "etags",
  data
) {

  results <- list()
  etag_id <- if (at_gene_level) "tgid" else "txid"
  etags <- selected_ids

  if (type == "etags") {
    for (metric_id in c("tax_id", "tax_name", "type")) {
      results[[metric_id]] <- table(
        annotation$generate_translate_dict(etag_id, metric_id)[etags]
      )
    }
  } else {
    if (is.null(data)) {
      logging::logerror("Unknown type for summary")
      stop()
    }
    for (metric_id in c("tax_id", "tax_name", "type")) {
      data[, metric_id] <-
        annotation$generate_translate_dict(
          etag_id, metric_id
        )[row.names(data)]
      results[[metric_id]] <-
        aggregate(formula(paste(". ~", metric_id)), data, sum)
      row.names(results[[metric_id]]) <- results[[metric_id]][, metric_id]
      results[[metric_id]][, metric_id] <- NULL
      data[, metric_id] <- NULL
    }
  }
  results
}

#' Helper function for ExprData$sum_per_type_per_sample
#'
#' @param annotation private$annotation (see ExprData)
#' @param main_etag private$main_etag (see ExprData)
#' @param design private$design (see ExprData)
#' @param get_raw_data_callback A function to get the raw data
#' @param compute_norm_callback A function to compute the normalized data
#' @param intra_norm see ExprData$sum_per_type_per_sample
#' @param log2_expr see ExprData$sum_per_type_per_sample
#' @return data frame with counts with columns type, sample, total
sum_per_type_per_sample_helper <- function(
  annotation, main_etag, design, get_raw_data_callback,
  compute_norm_callback, intra_norm = FALSE, log2_expr = FALSE
) {
  to_type <- annotation$generate_translate_dict(
    main_etag, "type"
  )
  format_data <- function(data) {
    data$type <- to_type[row.names(data)]
    temp <- tidyr::pivot_longer(
      aggregate(. ~ type, data, sum),
      cols = -type, names_to = "sample", values_to = "total"
    )
    subset(temp, total > 10)
  }

  if (intra_norm) {
    data <- purrr::map_dfr(
      design$list_batches(),
      function(batch) {
        format_data(as.data.frame(compute_norm_callback(
          in_batch = batch, intra_norm = intra_norm, inter_norm = FALSE
        )))
      }
    )
  } else {
    data <- format_data(as.data.frame(get_raw_data_callback()))
  }
  if (log2_expr) {
    data$total <- log2(data$total + 2)
  }
  data
}
