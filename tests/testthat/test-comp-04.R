# Setup block (run once for all tests in the file)
tmp_dir <- tempdir(TRUE)
# The withr::defer block ensures that the temporary directory is cleaned up
# even if the tests fail.
withr::defer(unlink(tmp_dir, recursive = TRUE))

test <- nexodiff::make_test_data("test04", tmp_dir)

annot <- nexodiff::Annotation$new(
  annotation = test$annotation,
  format_gff = FALSE,
  idmapping = test$id_mapping
)

design <- nexodiff::PairwiseDesign$new(
  pairwise_design_file = test$design,
  src_dir = test$src_dir
)

# Load the expression data
expr_data <- nexodiff::ExprDataTranscript$new(
  design,
  annot
)

# Filter and normalize
expr_data$filter_and_set_selected_ids(
  "human",
  filtered_var = "tax_name", filter_type = "keep"
)
expr_data$filter_and_set_selected_ids(
  "mRNA",
  filtered_var = "type", filter_type = "keep"
)
expr_data$compute_and_set_intra_norm_fact("tpm")
expr_data$compute_and_set_inter_norm_fact(ref_type = "ctrl")

# Initialize the DESeq2 comparison object
ds <- nexodiff::PairwiseDESeq2$new(expr_data)


testthat::test_that("PairwiseDESeq2 plots are generated correctly", {
  p_ma <- ds$plot_ma()
  vdiffr::expect_doppelganger("PairwiseDESeq2 MA Plot", p_ma)

  p_vulcano <- ds$plot_vulcano() + ggplot2::facet_wrap(batch ~ group)
  vdiffr::expect_doppelganger("PairwiseDESeq2 Faceted Vulcano Plot", p_vulcano)

  hm <- ds$plot_heatmap()
  vdiffr::expect_doppelganger("PairwiseDESeq2 Heatmap", hm$main)
})

testthat::test_that("Top gene selection and plotting works correctly", {
  # 1. Select top 10 deregulated genes
  top10_genes <- ds$generate_a_list(
    in_batch = "b1",
    in_group = "B",
    id = "symbol",
    type = "deregulated",
    top_x = 10
  )
  testthat::expect_length(top10_genes, 10)
  testthat::expect_type(top10_genes, "character")

  # 2. Test MA plot with selected genes highlighted
  p_ma_selected <- ds$plot_ma(
    in_batches = "b1",
    select_ids = top10_genes,
    tag_id_select = "symbol"
  )
  vdiffr::expect_doppelganger("MA Plot with Top 10 Selected", p_ma_selected)

  # 3. Test Heatmap with selected genes
  hm_selected <- ds$plot_heatmap(
    select_ids = top10_genes,
    tag_id_select = "symbol",
    in_batches = "b1"
  )
  vdiffr::expect_doppelganger("Heatmap of Top 10 Genes", hm_selected$main)

  # 4. Test LFC per group plots for selected genes
  p_lfc_bar <- ds$plot_lfc_per_group_facet_tags(
    select_ids = top10_genes,
    tag_id_select = "symbol",
    in_batches = "b1",
    geoms = c("bar", "errorbar")
  )
  vdiffr::expect_doppelganger("LFC Bar Plot for Top 10", p_lfc_bar)

})

testthat::test_that("DESeq2 results are consistent with simulated truth", {
  # 1. Get results from the DESeq2 analysis
  df <- ds$filter_and_get_results(
    in_batch = "b1", in_group = "B", add_ids = c("gid", "symbol")
  )
  df <- dplyr::arrange(df, as.integer(.data$gid))

  # 2. Load the ground truth data for this simulation
  tdf <- read.csv(file.path(test$config, "true_betas.csv"))

  # 3. Combine estimated results with simulated ground truth
  comp <- data.frame(
    base = df$baseMean,
    est = df$log2FoldChange,
    sim = tdf$true_beta_b1
  ) %>%
    dplyr::filter(
      !(is.na(est) | is.na(sim) | is.nan(est) | is.nan(sim))
    )

  # 4. Run simulations to establish baseline correlation and error metrics
  #    This helps determine if our results are better than random chance.
  x <- 10 # Number of repetitions
  beta_list <- lapply(1:x, function(i) {
    dds_demo <- DESeq2::makeExampleDESeqDataSet(
      n = 10000, m = 6,
      sizeFactors = runif(6, min = 0.8, max = 1.2),
      betaSD = 1
    )
    SummarizedExperiment::rowData(dds_demo)$trueBeta
  })

  combinations <- utils::combn(seq_along(beta_list), 2)
  sim_results <- apply(combinations, 2, function(pair) {
    beta1 <- beta_list[[pair[1]]]
    beta2 <- beta_list[[pair[2]]]
    c(
      correlation = stats::cor(beta1, beta2, method = "pearson"),
      mean_abs_err = mean(abs(beta1 - beta2))
    )
  })

  # 5. Perform the validation checks
  # We exclude the lowest 25% of expressed genes (q25) as they are often noisy
  comp_filtered <- dplyr::filter(comp, base >= quantile(comp$base, 0.25))

  # Calculate metrics from our test data
  est_correlation <- cor(
    comp_filtered$sim, comp_filtered$est,
    method = "pearson"
  )
  est_mae <- mean(abs(comp_filtered$sim - comp_filtered$est))

  # The correlation should be higher than the max correlation from random sims
  testthat::expect_true(
    est_correlation > max(sim_results["correlation", ]),
    label = paste(
      "Estimated correlation", round(est_correlation, 3),
      "is not greater than the simulation baseline max",
      round(max(sim_results["correlation", ]), 3)
    )
  )

  # The Mean Absolute Error should be lower than the min MAE from random sims
  testthat::expect_true(
    est_mae < min(sim_results["mean_abs_err", ]),
    label = paste(
      "Estimated MAE", round(est_mae, 3),
      "is not less than the simulation baseline min",
      round(min(sim_results["mean_abs_err", ]), 3)
    )
  )
})
