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
#' @param data_matrix expression data matrix (raw or normalized)
#' @param log2_expr see ExprData$sum_per_type_per_sample
#' @return data frame with counts with columns type, sample, total
summarize_per_type <- function(
  annotation, main_etag, data_matrix, log2_expr = FALSE
) {
  to_type <- annotation$generate_translate_dict(main_etag, "type")

  data <- as.data.frame(data_matrix)
  data$type <- to_type[row.names(data)]
  temp <- tidyr::pivot_longer(
    aggregate(. ~ type, data, sum),
    cols = -type, names_to = "sample", values_to = "total"
  )
  result <- subset(temp, total > 10)

  if (log2_expr) {
    result$total <- log2(result$total + 2)
  }
  result
}
