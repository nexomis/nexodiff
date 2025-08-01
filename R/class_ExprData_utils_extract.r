#' @include nexodiff-package.R
#' @include utils.r
NULL

#' Helper function for ExprData$extract_pairwise_data_with_design
#'
#' @param design private$design (see ExprData)
#' @param raw_matrix raw count matrix
#' @param len_matrix length matrix
#' @param intra_norm_fact_matrix intra normalization factor matrix
#' @param inter_norm_fact inter normalization factors
#' @param in_batch see ExprData$extract_pairwise_data_with_design
#' @param in_group see ExprData$extract_pairwise_data_with_design
#' @param include_ctrl see ExprData$extract_pairwise_data_with_design
#' @return list with selected attributes (see ExprData$extract_pairwise_data_with_design)
extract_pairwise_data <- function(
  design, raw_matrix, len_matrix, intra_norm_fact_matrix, inter_norm_fact,
  in_batch, in_group, include_ctrl = TRUE
) {
  group_per_batches <-
    design$list_groups_per_batches(include_ctrl = TRUE)

  if (
    (length(in_batch) != 1) ||
      (! in_batch %in% names(group_per_batches))
  ) {
    logging::logerror(
      "in_batch is not recognized in design"
    )
    stop()
  }
  if (
    (length(in_group) != 1) ||
      (! in_group %in% group_per_batches[in_batch][[1]])
  ) {
    logging::logerror(
      "in_group is not recognized in design"
    )
    stop()
  }

  results <- list()

  # test if control group is the same... Sometimes we use this to get the
  # for a control group also.
  results$test_samples <-
    design$extract_sample_names(in_batch, in_group)

  # find batch control
  if (include_ctrl) {
    results$ctrl_group <-
      design$find_control_group_per_batches()[in_batch]
    results$ctrl_samples <-
      design$extract_sample_names(in_batch, results$ctrl_group)
    all_groups <- c(in_group, results$ctrl_group)
    results$all_samples <-
      unique(c(results$test_samples, results$ctrl_samples))
  } else {
    all_groups <- in_group
    results$all_samples <- results$test_samples
  }

  # norm is not computed because if the table is big, it's too long

  results$raw <- raw_matrix[, results$all_samples]
  results$len <- len_matrix[, results$all_samples]

  results$intra_norm_fact <-
    intra_norm_fact_matrix[, results$all_samples]

  if (!is.null(inter_norm_fact)) {
    results$inter_norm_fact <- inter_norm_fact[results$all_samples]
  } else {
    results$inter_norm_fact <- rep(1, length(results$all_samples))
    names(results$inter_norm_fact) <- results$all_samples
  }

  results$design_table <- design$get_pairwise_design(
    in_batch, all_groups
  )

  row.names(results$design_table) <- results$design_table$sample
  results$design_table <- results$design_table[results$all_samples, ]
  results
}
