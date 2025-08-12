# Setup block (run once for all tests in the file)
tmp_dir <- tempdir(TRUE)
test <- nexodiff::make_test_data("test03", tmp_dir)

# Expected ratios from the simulation, which normalization should recover
exp_ratios <- readr::read_csv(
  file.path(test$config, "exp_ratios.csv"),
  show_col_types = FALSE
)
exp_ratios <- as.matrix(exp_ratios)

# --- Initialize objects for testing ---
annot <- nexodiff::Annotation$new(
  annotation = test$annotation,
  format_gff = FALSE,
  idmapping = test$id_mapping
)

design <- nexodiff::PairwiseDesign$new(
  pairwise_design_file = test$design,
  src_dir = test$src_dir
)

# Create a base ExprData object to be used across tests
expr_data_tx <- nexodiff::ExprDataTranscript$new(
  design,
  annot
)
# Filter to keep only relevant data for the test
expr_data_tx$filter_and_set_selected_ids(
  "human",
  filtered_var = "tax_name", filter_type = "keep"
)
expr_data_tx$filter_and_set_selected_ids(
  "mRNA",
  filtered_var = "type", filter_type = "keep"
)
# Use TPM for intra-normalization to account for library size differences
# simulated in the test data.
expr_data_tx$compute_and_set_intra_norm_fact("tpm")


# --- Tests ---

testthat::test_that("Inter-norm with norm_scale = 'design' works", {
  # Test with ref_type = "all" (balanced design)
  expr_data_tx$compute_and_set_inter_norm_fact(
    ref_type = "all", norm_scale = "design", norm_by = "sample"
  )
  norm_mat_all <- expr_data_tx$compute_norm(
    inter_norm = TRUE, rescale_inter_norm = FALSE
  )

  # Ratios relative to a stable control sample (A3). The biological ratio for
  # A3 is 1 for all genes, so this should recover the expected ratios.
  ratios_all <- sweep(norm_mat_all, 1, norm_mat_all[, "A3"], "/")

  testthat::expect_equal(
    ratios_all[!is.nan(ratios_all)],
    exp_ratios[!is.nan(ratios_all)], tolerance = 1e-5
  )

  # Test with ref_type = "ctrl"
  expr_data_tx$compute_and_set_inter_norm_fact(
    ref_type = "ctrl", norm_scale = "design", norm_by = "sample"
  )
  norm_mat_ctrl <- expr_data_tx$compute_norm(
    inter_norm = TRUE, rescale_inter_norm = FALSE
  )
  ratios_ctrl <- sweep(norm_mat_ctrl, 1, norm_mat_ctrl[, "A3"], "/")
  testthat::expect_equal(
    ratios_ctrl[!is.nan(ratios_ctrl)],
    exp_ratios[!is.nan(ratios_ctrl)], tolerance = 1e-5
  )
})

testthat::test_that("Inter-norm with norm_scale = 'batch' works", {
  expr_data_tx$compute_and_set_inter_norm_fact(
    ref_type = "ctrl", norm_scale = "batch", norm_by = "sample"
  )

  # Test batch 1 (reference is group A)
  norm_mat_b1 <- expr_data_tx$compute_norm(
    in_batch = "b1", inter_norm = TRUE, rescale_inter_norm = FALSE
  )
  b1_cols <- colnames(norm_mat_b1)
  ratios_b1 <- sweep(norm_mat_b1, 1, norm_mat_b1[, "A3"], "/")
  testthat::expect_equal(
    ratios_b1[!is.nan(ratios_b1)], exp_ratios[, b1_cols][!is.nan(ratios_b1)],
    tolerance = 1e-5
  )

  # Test batch 2 (reference is group D)
  norm_mat_b2 <- expr_data_tx$compute_norm(
    in_batch = "b2", inter_norm = TRUE, rescale_inter_norm = FALSE
  )
  b2_cols <- colnames(norm_mat_b2)
  ratios_b2 <- sweep(norm_mat_b2, 1, norm_mat_b2[, "D1"], "/")
  testthat::expect_equal(
    ratios_b2[!is.nan(ratios_b2)], exp_ratios[, b2_cols][!is.nan(ratios_b2)],
    tolerance = 1e-5
  )
})

testthat::test_that("Inter-norm with norm_scale = 'group' works", {
  expr_data_tx$compute_and_set_inter_norm_fact(
    ref_type = "ctrl", norm_scale = "group", norm_by = "sample"
  )

  # Test group B in batch b1 (control is A)
  norm_mat_b1B <- expr_data_tx$compute_norm(
    in_batch = "b1", in_group = "B", inter_norm = TRUE,
    include_ctrl = TRUE, rescale_inter_norm = FALSE
  )
  b1B_cols <- colnames(norm_mat_b1B)
  ratios_b1B <- sweep(norm_mat_b1B, 1, norm_mat_b1B[, "A3"], "/")
  testthat::expect_equal(
    ratios_b1B[!is.nan(ratios_b1B)],
    exp_ratios[, b1B_cols][!is.nan(ratios_b1B)],
    tolerance = 1e-5
  )

  # Test group E in batch b2 (control is D)
  norm_mat_b2E <- expr_data_tx$compute_norm(
    in_batch = "b2", in_group = "E", inter_norm = TRUE,
    include_ctrl = TRUE, rescale_inter_norm = FALSE
  )
  b2E_cols <- colnames(norm_mat_b2E)
  ratios_b2E <- sweep(norm_mat_b2E, 1, norm_mat_b2E[, "D1"], "/")
  testthat::expect_equal(
    ratios_b2E[!is.nan(ratios_b2E)],
    exp_ratios[, b2E_cols][!is.nan(ratios_b2E)],
    tolerance = 1e-5
  )
})

testthat::test_that("Rescaling of inter-normalization factors works", {
  # When rescale=TRUE (default), the geometric mean of factors should be 1
  expr_data_tx$compute_and_set_inter_norm_fact(
    ref_type = "all", norm_scale = "design"
  )

  # For design scale
  factors_design_rescaled <- expr_data_tx$get_inter_norm_fact(rescale = TRUE)
  testthat::expect_equal(exp(mean(log(factors_design_rescaled))), 1)

  factors_design_raw <- expr_data_tx$get_inter_norm_fact(rescale = FALSE)
  testthat::expect_true(abs(exp(mean(log(factors_design_raw))) - 1) > 1e-9)

  # For batch scale
  expr_data_tx$compute_and_set_inter_norm_fact(
    ref_type = "ctrl", norm_scale = "batch"
  )

  factors_b1_rescaled <- expr_data_tx$get_inter_norm_fact(
    in_batch = "b1", rescale = TRUE
  )
  testthat::expect_equal(exp(mean(log(factors_b1_rescaled))), 1)

  factors_b1_raw <- expr_data_tx$get_inter_norm_fact(
    in_batch = "b1", rescale = FALSE
  )
  testthat::expect_true(abs(exp(mean(log(factors_b1_raw))) - 1) > 1e-9)
})

testthat::test_that("Inter-normalization error handling works", {
  # Invalid arguments for compute_and_set_inter_norm_fact
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(norm_scale = "invalid_scale")
  )
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(norm_by = "invalid_by")
  )
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(ref_type = "invalid_ref")
  )
  testthat::expect_error(
    expr_data_tx$compute_and_set_inter_norm_fact(m_trim_prop = 0.6)
  )

  # Calling compute_norm with inter_norm=T but intra_norm=F
  testthat::expect_error(
    expr_data_tx$compute_norm(inter_norm = TRUE, intra_norm = FALSE),
    "`intra_norm` cannot be FALSE if `inter_norm` is TRUE"
  )

  # Calling get_inter_norm_fact before computing it
  expr_data_tx_fresh <- expr_data_tx$clone(deep = TRUE)
  expr_data_tx_fresh$reset() # reset removes inter_norm_fact
  testthat::expect_error(
    expr_data_tx_fresh$get_inter_norm_fact()
  )
})