test04 <- nexodiff::make_test_data_from_xlsx("test04")

test <- test04

design <- nexodiff::PairwiseDesignWithAnnotation$new(
    pairwise_design_file = test$design,
    annotation_file = test$annotation,
    idmapping_file = test$id_mapping,
    src_dir = test$src_dir
)

# Load the expression data based on the

expr_data_tx <- nexodiff::ExprDataTranscript$new(design)
expr_data_tx$filter_and_set_selected_ids("human", filtered_var = "tax_name")
expr_data_tx$filter_and_set_selected_ids("mRNA", filtered_var = "type")


expr_data_tx$compute_and_set_intra_norm_fact("fpk")

expr_data_tx$compute_and_set_inter_norm_fact(
  method = "median", norm_scale = "design")
exp_norm_design <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_fpk_mrn_bysp_rfall_desgn",
  colNames = TRUE
)
colnames(exp_norm_design) <- paste("run1", colnames(exp_norm_design), sep = "_")
exp_norm_design <- as.matrix(exp_norm_design)
norm_design <- expr_data_tx$compute_norm(inter_norm = TRUE)
exp_ratio_design <- norm_design
exp_ratio_design[] <- 1

expr_data_tx$compute_and_set_inter_norm_fact(
  method = "median", norm_scale = "batch")
exp_norm_batch <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_fpk_mrn_bysp_rfall_batch",
  colNames = TRUE
)
colnames(exp_norm_batch) <- paste("run1", colnames(exp_norm_batch), sep = "_")
exp_norm_batch <- as.matrix(exp_norm_batch)
norm_design_b1 <- expr_data_tx$compute_norm(in_batch = "b1", inter_norm = TRUE)
norm_design_b2 <- expr_data_tx$compute_norm(in_batch = "b2", inter_norm = TRUE)
b1_cols <- colnames(norm_design_b1)
b2_cols <- colnames(norm_design_b2)

expr_data_tx$compute_and_set_inter_norm_fact(
  method = "median", norm_scale = "group")

exp_norm_b1A <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_fpk_mrn_bysp_rfall_grp_b1A",
  colNames = TRUE
)
colnames(exp_norm_b1A) <- paste("run1", colnames(exp_norm_b1A), sep = "_")
exp_norm_b1A <- as.matrix(exp_norm_b1A)
b1A_cols <- colnames(exp_norm_b1A)
norm_design_b1A <- expr_data_tx$compute_norm(in_batch = "b1", in_group = "A",
  inter_norm = TRUE, include_ctrl = TRUE)
exp_norm_b2C <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_fpk_mrn_bysp_rfall_grp_b2C",
  colNames = TRUE
)
colnames(exp_norm_b2C) <- paste("run1", colnames(exp_norm_b2C), sep = "_")
exp_norm_b2C <- as.matrix(exp_norm_b2C)
b2C_cols <- colnames(exp_norm_b2C)
norm_design_b2C <- expr_data_tx$compute_norm(in_batch = "b2", in_group = "C",
  inter_norm = TRUE, include_ctrl = TRUE)
exp_norm_b2C <- as.matrix(exp_norm_b2C)
exp_norm_b2E <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_fpk_mrn_bysp_rfall_grp_b2E",
  colNames = TRUE
)
colnames(exp_norm_b2E) <- paste("run1", colnames(exp_norm_b2E), sep = "_")
exp_norm_b2E <- as.matrix(exp_norm_b2E)
b2E_cols <- colnames(exp_norm_b2E)
exp_norm_b2E <- as.matrix(exp_norm_b2E)
norm_design_b2E <- expr_data_tx$compute_norm(in_batch = "b2", in_group = "E",
  inter_norm = TRUE, include_ctrl = TRUE)

testthat::test_that("inter norm scale tests", {
  testthat::expect_identical(
    exp_ratio_design,
    norm_design / exp_norm_design,
    tolerance = 1e-5
  )
  testthat::expect_identical(
    exp_ratio_design[, b2_cols],
    norm_design_b2[, b2_cols] / exp_norm_batch[, b2_cols],
    tolerance = 1e-5
  )
  testthat::expect_identical(
    exp_ratio_design[, b2_cols],
    norm_design_b2[, b2_cols] / exp_norm_batch[, b2_cols],
    tolerance = 1e-5
  )
  testthat::expect_identical(
    exp_ratio_design[, b1A_cols],
    norm_design_b1A[, b1A_cols] / exp_norm_b1A[, b1A_cols],
    tolerance = 1e-5
  )
  testthat::expect_identical(
    exp_ratio_design[, b2C_cols],
    norm_design_b2C[, b2C_cols] / exp_norm_b2C[, b2C_cols],
    tolerance = 1e-5
  )
  testthat::expect_identical(
    exp_ratio_design[, b2E_cols],
    norm_design_b2E[, b2E_cols] / exp_norm_b2E[, b2E_cols],
    tolerance = 1e-5
  )
})

testthat::test_that("inter error handling tests", {
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(
      method = "tmm")
  )
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(
      method = "tmm", m_trim_mean = "median")
  )
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(
      method = "tmm", m_trim_prop = 0)
  )
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(
      method = "tmm", m_trim_prop = 0.6)
  )
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(
      norm_by = "group"
    )
  )
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(
      norm_ref = "YO"
    )
  )
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(
      norm_scale = "design",
      norm_ref = "ctrl"
    )
  )
})
