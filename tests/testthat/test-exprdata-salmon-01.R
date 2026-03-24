
tmp_dir <- tempdir(TRUE)
test <- nexodiff::make_test_data("test01", tmp_dir, quant_source = "salmon")

annot <- nexodiff::Annotation$new(
  annotation = test$annotation,
  format_gff = FALSE,
  idmapping = test$id_mapping
)

design <- nexodiff::PairwiseDesign$new(
  pairwise_design_file = test$design,
  src_dir = test$src_dir,
  quant_source = "salmon"
)

# Load the expression data based on the

expr_data_tx <- nexodiff::ExprDataTranscript$new(
  design,
  annot
)

expr_data_tx_fixed <- nexodiff::ExprDataTranscript$new(
  design,
  annot,
  with_fixed_length = TRUE
)

expected_raw <- readr::read_csv(
  file.path(test$config, "tx_raw.csv"),
  show_col_types = FALSE
)

expected_len <- readr::read_csv(
  file.path(test$config, "tx_len.csv"),
  show_col_types = FALSE
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

testthat::test_that("loading tx data from salmon test", {
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
