tmp_dir <- tempdir(TRUE)
test <- nexodiff::make_test_data("test02", tmp_dir)

annot <- nexodiff::Annotation$new(
  annotation = test$annotation,
  format_gff = FALSE,
  idmapping = test$id_mapping
)

design <- nexodiff::PairwiseDesign$new(
  pairwise_design_file = test$design,
  src_dir = test$src_dir
)

# Load the expression data based on the
expr_data_tx <- nexodiff::ExprDataTranscript$new(
  design,
  annot
)

# Test show_etags_summary with default parameters (raw counts)
summary_default <- expr_data_tx$show_etags_summary()

# Test show_etags_summary with transformation and sum functions
summary_transformed <- expr_data_tx$show_etags_summary(
  tr_fn = function(x) log2(x + 2),
  sum_fn = sum
)

# Test show_etags_summary with count function
summary_count <- expr_data_tx$show_etags_summary(
  tr_fn = NULL,
  sum_fn = function(x) sum(as.integer(x > 0))
)

summary_count_transformed <- expr_data_tx$show_etags_summary(
  tr_fn = function(x) log2(x + 2),
  sum_fn = function(x) sum(as.integer(x > 1))
)

# Test filtering by species (tax_id)
# Filter to keep only human data (tax_id = 9606)
expr_data_tx$filter_and_set_selected_ids(c("9606"), "tax_id", "keep")
summary_human <- expr_data_tx$show_etags_summary()

# Reset and filter to exclude human data
expr_data_tx$reset()
expr_data_tx$filter_and_set_selected_ids(c("9606"), "tax_id", "excl")
summary_non_human <- expr_data_tx$show_etags_summary()

# Reset and filter by RNA type
expr_data_tx$reset()
expr_data_tx$filter_and_set_selected_ids(c("mRNA"), "type", "keep")
summary_mrna <- expr_data_tx$show_etags_summary()

# Reset and filter by multiple RNA types
expr_data_tx$reset()
expr_data_tx$filter_and_set_selected_ids(c("mRNA", "rRNA"), "type", "keep")
summary_all_rna <- expr_data_tx$show_etags_summary()

# Test with normalized data
expr_data_tx$reset()
expr_data_tx$compute_and_set_intra_norm_fact("tpm")
summary_norm <- expr_data_tx$show_etags_summary(intra_norm = TRUE)

testthat::test_that("show_etags_summary works with different parameters", {
  # Check that we get results for all metrics
  testthat::expect_true(
    all(c("tax_id", "tax_name", "type") %in% names(summary_default))
  )

  # Check that results are matrices or data frames
  testthat::expect_true(is.data.frame(summary_default$tax_id))
  testthat::expect_true(is.data.frame(summary_default$tax_name))
  testthat::expect_true(is.data.frame(summary_default$type))

  # Check that transformed results are different from default
  testthat::expect_false(
    identical(summary_default$tax_id, summary_transformed$tax_id)
  )

  # Check that length results are different from sum results
  testthat::expect_false(
    identical(summary_default$tax_id, summary_count$tax_id)
  )

  # Check filtering results
  testthat::expect_true(
    nrow(summary_human$tax_id) < nrow(summary_default$tax_id)
  )
  testthat::expect_true(
    nrow(summary_non_human$tax_id) < nrow(summary_default$tax_id)
  )
  testthat::expect_true(
    nrow(summary_mrna$type) < nrow(summary_default$type)
  )

  # Check that human filter only contains human data
  testthat::expect_true(all(rownames(summary_human$tax_id) == "9606"))

  # Check that mRNA filter only contains mRNA data
  testthat::expect_true(all(rownames(summary_mrna$type) == "mRNA"))
})

# Test gene-level expression data
expr_data_tx$reset()
expr_data <- nexodiff::ExprDataGene$new(expr_data_tx)

summary_gene <- expr_data$show_etags_summary()
summary_gene_transformed <- expr_data$show_etags_summary(
  tr_fn = function(x) log2(x + 1),
  sum_fn = sum
)

testthat::test_that("gene-level show_etags_summary works", {
  # Check that we get results for all metrics
  testthat::expect_true(
    all(c("tax_id", "tax_name", "type") %in% names(summary_gene))
  )

  # Check that results are matrices or data frames
  testthat::expect_true(is.data.frame(summary_gene$tax_id))
  testthat::expect_true(is.data.frame(summary_gene$tax_name))
  testthat::expect_true(is.data.frame(summary_gene$type))
})

testthat::test_that("summary values are ok", {
  # Check tax_name dimension names
  testthat::expect_equal(
    dimnames(summary_default$tax_name)[[1]],
    c("human", "synthetic_construct")
  )

  # Check tax_id values
  testthat::expect_equal(
    summary_default$tax_id[1, "run1_A1"],
    700
  )
  testthat::expect_equal(
    summary_default$tax_id[1, "run1_B1"],
    700
  )
  testthat::expect_equal(
    summary_default$tax_id[2, "run1_A1"],
    1200
  )
  testthat::expect_equal(
    summary_default$tax_id[2, "run1_B1"],
    1200
  )

  # Check type values
  testthat::expect_equal(
    summary_default$type[1, "run1_A1"],
    1250
  )
  testthat::expect_equal(
    summary_default$type[1, "run1_B1"],
    1250
  )
  testthat::expect_equal(
    summary_default$type[2, "run1_A1"],
    650
  )
  testthat::expect_equal(
    summary_default$type[2, "run1_B1"],
    650
  )

  # Check filtered results row names
  testthat::expect_equal(
    dimnames(summary_mrna$type)[[1]],
    "mRNA"
  )
  testthat::expect_equal(
    dimnames(summary_all_rna$type)[[1]],
    c("mRNA", "rRNA")
  )
  testthat::expect_equal(
    dimnames(summary_non_human$tax_name)[[1]],
    "synthetic_construct"
  )
  testthat::expect_equal(
    dimnames(summary_human$tax_name)[[1]],
    "human"
  )
})
