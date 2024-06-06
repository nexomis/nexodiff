

test <- nexodiff::make_test_data_from_xlsx("test01")

design <- nexodiff::PairwiseDesignWithAnnotation$new(
    pairwise_design_file = test$design,
    annotation_file = test$annotation,
    idmapping_file = test$id_mapping,
    src_dir = test$src_dir
)

# Load the expression data based on the

expr_data_tx <- nexodiff::ExprDataTranscript$new(design)

expr_data_tx_fixed <- nexodiff::ExprDataTranscript$new(
  design,
  with_fixed_length = TRUE
)

expected_raw <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_raw",
  colNames = TRUE
)

expected_len <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_len",
  colNames = TRUE
)

raw <- expr_data_tx$get_raw()
expected_raw <- as.matrix(expected_raw)
mode(raw) <- mode(expected_raw)
dimnames(expected_raw) <- dimnames(raw)

len <- expr_data_tx$get_len()
expected_len <- as.matrix(expected_len)
mode(len) <- mode(expected_len)
dimnames(expected_len) <- dimnames(len)

len_fixed <- expr_data_tx_fixed$get_len()
expected_len_fixed <- expected_len
expected_len_fixed[] <- 1
mode(len_fixed) <- mode(expected_len_fixed)
dimnames(expected_len_fixed) <- dimnames(len_fixed)

testthat::test_that("loading tx data from test", {
  testthat::expect_identical(
    raw,
    expected_raw
  )
  testthat::expect_identical(
    len,
    expected_len
  )
  testthat::expect_identical(
    len_fixed,
    expected_len_fixed
  )
})

expected_tpm <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_tpm",
  colNames = TRUE
)

expected_fpkm <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_fpkm",
  colNames = TRUE
)

expected_fpm <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_fpm",
  colNames = TRUE
)

expected_fpk <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "tx_fpk",
  colNames = TRUE
)

expr_data_tx$compute_and_set_intra_norm_fact("tpm")
tpm <- expr_data_tx$compute_norm()
expected_tpm <- as.matrix(expected_tpm)
mode(tpm) <- mode(expected_tpm)
dimnames(expected_tpm) <- dimnames(tpm)

expr_data_tx$reset()
norm_reset <- expr_data_tx$compute_norm()
mode(norm_reset) <- mode(expected_raw)

expr_data_tx$compute_and_set_intra_norm_fact("fpkm")
fpkm <- expr_data_tx$compute_norm()
expected_fpkm <- as.matrix(expected_fpkm)
mode(fpkm) <- mode(expected_fpkm)
dimnames(expected_fpkm) <- dimnames(fpkm)

expr_data_tx$compute_and_set_intra_norm_fact("fpm")
fpm <- expr_data_tx$compute_norm()
expected_fpm <- as.matrix(expected_fpm)
mode(fpm) <- mode(expected_fpm)
dimnames(expected_fpm) <- dimnames(fpm)

expr_data_tx$compute_and_set_intra_norm_fact("fpk")
fpk <- expr_data_tx$compute_norm()
expected_fpk <- as.matrix(expected_fpk)
mode(fpk) <- mode(expected_fpk)
dimnames(expected_fpk) <- dimnames(fpk)

expr_data_tx$compute_and_set_intra_norm_fact("none")
nn <- expr_data_tx$compute_norm()
mode(nn) <- mode(expected_raw)

testthat::test_that("normalization at tx level", {
  testthat::expect_identical(
    nn,
    expected_raw
  )
  testthat::expect_identical(
    norm_reset,
    expected_raw
  )
  testthat::expect_identical(
    tpm,
    expected_tpm,
    tolerance = 1e-5
  )
  testthat::expect_identical(
    fpkm,
    expected_fpkm,
    tolerance = 1e-5
  )
  testthat::expect_identical(
    fpm,
    expected_fpm,
    tolerance = 1e-5
  )
  testthat::expect_identical(
    fpk,
    expected_fpk,
    tolerance = 1e-5
  )
  testthat::expect_error(
    expr_data_tx$compute_and_set_intra_norm_fact("unknown method")
  )
})

expr_data_tx$reset()

expr_data <- nexodiff::ExprDataGene$new(expr_data_tx)

expected_graw <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "g_raw",
  colNames = TRUE
)

expected_glen <- openxlsx::read.xlsx(
  xlsxFile = test$config,
  sheet = "g_len",
  colNames = TRUE
)

graw <- expr_data$get_raw()
expected_graw <- as.matrix(expected_graw)
mode(graw) <- mode(expected_graw)
dimnames(expected_graw) <- dimnames(graw)

glen <- expr_data$get_len()
expected_glen <- as.matrix(expected_glen)
mode(glen) <- mode(expected_glen)
dimnames(expected_glen) <- dimnames(glen)

testthat::test_that("summarisation", {
  testthat::expect_identical(
    graw,
    expected_graw
  )
  testthat::expect_identical(
    glen,
    expected_glen,
    tolerance = 1e-5
  )
})
