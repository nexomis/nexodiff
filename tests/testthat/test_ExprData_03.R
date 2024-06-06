test03 <- nexodiff::make_test_data_from_xlsx("test03")

test <- test03

design <- nexodiff::PairwiseDesignWithAnnotation$new(
    pairwise_design_file = test$design,
    annotation_file = test$annotation,
    idmapping_file = test$id_mapping,
    src_dir = test$src_dir
)

# Load the expression data based on the

expr_data_tx <- nexodiff::ExprDataTranscript$new(design)
expr_data <- nexodiff::ExprDataGene$new(expr_data_tx)
expr_data$compute_and_set_intra_norm_fact("tpm")

g1 <- expr_data$plot_sum_per_type_per_sample(
  intra_norm = FALSE,
  exclude_type = c(),
  same_scale = TRUE,
  horizonal_bar = TRUE,
  log2_expr = FALSE
)

g2 <- expr_data$plot_sum_per_type_per_sample(
  intra_norm = TRUE,
  exclude_type = c("mRNA"),
  same_scale = FALSE,
  horizonal_bar = FALSE,
  log2_expr = TRUE
)

g3 <- expr_data$plot_sum_per_type_per_sample(
  intra_norm = FALSE,
  exclude_type = c(),
  same_scale = TRUE,
  horizonal_bar = FALSE,
  log2_expr = FALSE
)

pdf(NULL)

testthat::test_that("plot_type_sum_per_sample tests", {
  testthat::expect_true(
    ggplot2::is.ggplot(g1)
  )
  testthat::expect_true(
    ggplot2::is.ggplot(g2)
  )
  testthat::expect_true(
    ggplot2::is.ggplot(g3)
  )
  testthat::expect_error(
    g1,
    NA
  )
  testthat::expect_error(
    g2,
    NA
  )
  testthat::expect_error(
    g3,
    NA
  )
})

g1 <- expr_data$plot_dist_per_sample(intra_norm = TRUE, log2_expr = TRUE)
g2 <- expr_data$plot_dist_per_sample(intra_norm = FALSE, log2_expr = FALSE)

pdf(NULL)


testthat::test_that("plot_dist_per_sample tests", {
  testthat::expect_true(
    ggplot2::is.ggplot(g1)
  )
  testthat::expect_true(
    ggplot2::is.ggplot(g2)
  )
  testthat::expect_error(
    g1,
    NA
  )
  testthat::expect_error(
    g2,
    NA
  )
})

dev.off()

res <- expr_data$extract_pairwise_data_with_design("b1", "A")

res_no_ctrl_a <- expr_data$extract_pairwise_data_with_design("b1", "A",
  include_ctrl = FALSE)

res_no_ctrl_b <- expr_data$extract_pairwise_data_with_design("b1", "B",
  include_ctrl = FALSE)

testthat::test_that("get_pairwise_data_with_design tests", {
  testthat::expect_identical(
    colnames(res$raw),
    c("run1_A1", "run1_A2", "run1_B1", "run1_B2")
  )
  testthat::expect_identical(
    res$test_samples,
    c("run1_A1", "run1_A2")
  )
  testthat::expect_identical(
    res$ctrl_samples,
    c("run1_B1", "run1_B2")
  )
  testthat::expect_identical(
    colnames(res_no_ctrl_a$raw),
    c("run1_A1", "run1_A2")
  )
  testthat::expect_identical(
    colnames(res_no_ctrl_b$raw),
    c("run1_B1", "run1_B2")
  )
  testthat::expect_identical(
    res_no_ctrl_b$test_samples,
    c("run1_B1", "run1_B2")
  )
  testthat::expect_identical(
    colnames(res_no_ctrl_b$intra_norm_fact),
    c("run1_B1", "run1_B2")
  )
  testthat::expect_identical(
    names(res_no_ctrl_b$inter_norm_fact),
    c("run1_B1", "run1_B2")
  )
  testthat::expect_error(
    expr_data$extract_pairwise_data_with_design("b3", "A")
  )
  testthat::expect_error(
    expr_data$extract_pairwise_data_with_design("b1", "Y")
  )
})

expr_data_tx$reset()

expr_data <- nexodiff::ExprDataGene$new(expr_data_tx)
expr_data$compute_and_set_intra_norm_fact("tpm")
res <- expr_data$extract_pairwise_data_with_design("b1", "A")

expr_data <- nexodiff::ExprDataGene$new(expr_data_tx)
expr_data$compute_and_set_intra_norm_fact("tpm")
