root_path <- system.file("extdata/test_data1", package = "nexodiff")

design <- nexodiff::PairwiseDesignWithAnnotation$new(
    pairwise_design_file = paste(root_path, "design.csv", sep = "/"),
    annotation_file = paste(root_path, "annotation.txt", sep = "/"),
    idmapping_file = paste(root_path, "id_mapping.tab", sep = "/"),
    src_dir = root_path
)

expr_data_tx <- nexodiff::ExprDataTranscript$new(design)
expr_data_tx$filter_and_set_selected_ids("human", filtered_var = "tax_name")
expr_data_tx$filter_and_set_selected_ids("mRNA", filtered_var = "type")
expr_data <- nexodiff::ExprDataGene$new(expr_data_tx)
expr_data$compute_and_set_intra_norm_fact("fpk")
expr_data$compute_and_set_inter_norm_fact(
  method = "median", norm_scale = "design")
deseq <- nexodiff::PairwiseDESeq2$new(expr_data, ncpus = 8)
deg <- deseq$cross_args_and_generate_lists()


library(magrittr)
data_summary <- deseq$generate_summary() %>%
  dplyr::filter(.data$lfc_abs_lim == 1 & .data$min_signif == 0.01)

n_deg <- dplyr::filter(data_summary, .data$type == "deregulated")$n_deg
n_up <- dplyr::filter(data_summary, .data$type == "upregulated")$n_deg
n_down <- dplyr::filter(data_summary, .data$type == "downregulated")$n_deg

n_samples <- as.integer(deseq$get_results()[, c("n_test", "n_ctrl")])

deseq_paired <- nexodiff::PairwiseDESeq2$new(
  expr_data,
  add_vars = c("paired_id1"),
  ncpus = 8,
  min_base_mean = 100
)

deg_paired <- deseq_paired$cross_args_and_generate_lists()
data_summary_paired <- deseq_paired$generate_summary() %>%
  dplyr::filter(.data$lfc_abs_lim == 1 & .data$min_signif == 0.01)
n_deg_paired <-
  dplyr::filter(data_summary_paired, .data$type == "deregulated")$n_deg
n_up_paired <-
  dplyr::filter(data_summary_paired, .data$type == "upregulated")$n_deg
n_down_paired <-
  dplyr::filter(data_summary_paired, .data$type == "downregulated")$n_deg
n_samples_paired <-
  as.integer(deseq_paired$get_results()[, c("n_test", "n_ctrl")])

deseq_paired_only <- nexodiff::PairwiseDESeq2$new(
  expr_data,
  add_vars = c("paired_id1"),
  ncpus = 8,
  only_paired = TRUE,
  min_base_mean = 100
)

deg_paired_only <- deseq_paired_only$cross_args_and_generate_lists()
data_summary_paired_only <- deseq_paired_only$generate_summary() %>%
  dplyr::filter(.data$lfc_abs_lim == 1 & .data$min_signif == 0.01)
n_deg_paired_only <-
  dplyr::filter(data_summary_paired_only, .data$type == "deregulated")$n_deg
n_up_paired_only <-
  dplyr::filter(data_summary_paired_only, .data$type == "upregulated")$n_deg
n_down_paired_only <-
  dplyr::filter(data_summary_paired_only, .data$type == "downregulated")$n_deg
n_samples_paired_only <-
  as.integer(deseq_paired_only$get_results()[, c("n_test", "n_ctrl")])


deseq_bis <- nexodiff::PairwiseDESeq2$new(expr_data, ncpus = 8)

testthat::test_that("unit tests on pairwiseComp Object", {
  testthat::expect_identical(
    deseq,
    deseq_bis
  )
  testthat::expect_identical(
    n_deg,
    500,
    tolerance = 25
  )
  testthat::expect_identical(
    n_up,
    250,
    tolerance = 25
  )
  testthat::expect_identical(
    n_down,
    250,
    tolerance = 25
  )
  testthat::expect_identical(
    n_samples,
    c(10L, 10L)
  )
  testthat::expect_identical(
    n_deg_paired,
    500,
    tolerance = 25
  )
  testthat::expect_identical(
    n_up_paired,
    250,
    tolerance = 25
  )
  testthat::expect_identical(
    n_down_paired,
    250,
    tolerance = 25
  )
  testthat::expect_identical(
    n_samples_paired,
    c(10L, 10L)
  )
  testthat::expect_identical(
    n_deg_paired_only,
    500,
    tolerance = 25
  )
  testthat::expect_identical(
    n_up_paired_only,
    250,
    tolerance = 25
  )
  testthat::expect_identical(
    n_down_paired_only,
    250,
    tolerance = 25
  )
  testthat::expect_identical(
    n_samples_paired_only,
    c(8L, 8L)
  )
})




