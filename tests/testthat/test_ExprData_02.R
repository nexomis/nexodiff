
test02 <- nexodiff::make_test_data_from_xlsx("test02")

test <- test02

design <- nexodiff::PairwiseDesignWithAnnotation$new(
    pairwise_design_file = test$design,
    annotation_file = test$annotation,
    idmapping_file = test$id_mapping,
    src_dir = test$src_dir
)

# Load the expression data based on the

expr_data_tx <- nexodiff::ExprDataTranscript$new(design)


testthat::test_that("testing etags type counts", {
  testthat::expect_identical(
    as.matrix(expr_data_tx$show_etags_summary()$tax_id),
    as.matrix(c(
    "32630" = 4L,
    "9606" = 15L
    ))
  )
  testthat::expect_identical(
    as.matrix(expr_data_tx$show_etags_summary()$tax_name),
    as.matrix(c(
    "human" = 15L,
    "synthetic_construct" = 4L
    ))
  )
  testthat::expect_identical(
    as.matrix(expr_data_tx$show_etags_summary()$type),
    as.matrix(c(
    "mRNA" = 14L,
    "rRNA" = 5L
    ))
  )
})

expr_data_tx$filter_and_set_selected_ids("rRNA", filtered_var = "type",
  filter_type = "excl")

testthat::test_that("testing etags type counts after filtering 1/3", {
  testthat::expect_identical(
    as.matrix(expr_data_tx$show_etags_summary()$type),
    as.matrix(c(
    "mRNA" = 14L
    ))
  )
})

expr_data_tx$reset()
expr_data_tx$filter_and_set_selected_ids("human", filtered_var = "tax_name")

testthat::test_that("testing etags type counts after filtering 2/3", {
  testthat::expect_identical(
    as.matrix(expr_data_tx$show_etags_summary()$tax_name),
    as.matrix(c(
    "human" = 15L
    ))
  )
})

expr_data_tx$filter_and_set_selected_ids("mRNA", filtered_var = "type")

testthat::test_that("testing etags type counts after filtering 3/3", {
  testthat::expect_identical(
    as.matrix(expr_data_tx$show_etags_summary()$tax_name),
    as.matrix(c(
    "human" = 10L
    ))
  )
    testthat::expect_identical(
    as.matrix(expr_data_tx$show_etags_summary()$type),
    as.matrix(c(
    "mRNA" = 10L
    ))
  )
})

expr_data_tx$reset()
expr_data <- nexodiff::ExprDataGene$new(expr_data_tx)
expr_data$compute_and_set_intra_norm_fact("tpm")
expr_data$compute_and_set_intra_norm_fact("none")
raw_values <- as.numeric(unlist(expr_data$show_etags_summary(type = "raw")))
names(raw_values) <- NULL

norm_values <- as.numeric(
    unlist(expr_data$show_etags_summary(type = "norm", in_batch = "b1"))
)
names(norm_values) <- NULL

expected <- c(700, 1200, 700, 1200, 1200, 700, 1200, 700, 1250, 650, 1250, 650)

testthat::test_that("misc tests", {
  testthat::expect_identical(
    raw_values,
    expected
  )
  testthat::expect_identical(
    norm_values,
    expected
  )
  testthat::expect_error(
    expr_data$show_etags_summary(type = "unknown")
  )
  testthat::expect_error(
    expr_data$filter_and_set_selected_ids(c("a", "b"),
    filtered_var  = "unknown")
  )
  testthat::expect_error(
    expr_data$filter_and_set_selected_ids(c("a", "b"), filter_type  = "unknown")
  )
  testthat::expect_identical(
    expr_data$get_intra_norm_fact_method(),
    "none"
  )
  testthat::expect_identical(
    expr_data$is_at_gene_level(),
    TRUE
  )
  testthat::expect_identical(
    expr_data_tx$is_at_gene_level(),
    FALSE
  )
})

expr_data$reset()
expr_data$compute_and_set_intra_norm_fact("tpm")
raw_sum_per_sample <- expr_data$sum_per_type_per_sample()
lograw_sum_per_sample <- expr_data$sum_per_type_per_sample(log2_expr = TRUE)
norm_sum_per_sample <- expr_data$sum_per_type_per_sample(intra_norm = TRUE)

testthat::test_that("get_type_sum_per_sample tests", {
  testthat::expect_identical(
    raw_sum_per_sample$total,
    c(1250, 1250, 650, 650)
  )
  testthat::expect_identical(
    log2(raw_sum_per_sample$total + 2),
    lograw_sum_per_sample$total
  )
  testthat::expect_identical(
    sum(norm_sum_per_sample$total),
    2e+06,
    tolerance = 1e-5
  )
})