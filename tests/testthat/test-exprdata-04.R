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

# Load the expression data based on the
expr_data <- nexodiff::ExprDataTranscript$new(
  design,
  annot
)

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


testthat::test_that("ExprData plot_dist_per_sample works correctly", {
  p_dist_boxplot <- expr_data$plot_dist_per_sample(
    intra_norm = TRUE,
    inter_norm = TRUE,
    geoms = c("violin", "histo", "boxplot")
  )
  vdiffr::expect_doppelganger(
    "Sample Distribution (Vln, Histo, Box)",
    p_dist_boxplot
  )

  p_dist_violin <- expr_data$plot_dist_per_sample(
    intra_norm = TRUE,
    inter_norm = TRUE,
    geoms = c("violin")
  )
  vdiffr::expect_doppelganger(
    "Sample Distribution (Violin only)",
    p_dist_violin
  )
})

testthat::test_that("ExprData complex plots work correctly", {
  # PCA plots
  p_pca_design <- expr_data$plot_prcomp(
    plot_scale = "design",
    pca_plot_dims = c(1, 2, 3),
    intra_norm = TRUE,
    inter_norm = TRUE,
    point_size = 5
  )
  vdiffr::expect_doppelganger("PCA Plot (design scale)", p_pca_design)

  p_pca_batch <- expr_data$plot_prcomp(
    plot_scale = "batch",
    intra_norm = TRUE,
    inter_norm = TRUE
  )
  vdiffr::expect_doppelganger("PCA Plot (batch scale)", p_pca_batch)

  p_pca_group <- expr_data$plot_prcomp(
    in_batch = "b1",
    pca_plot_dims = c(1, 2, 3),
    plot_scale = "group",
    intra_norm = TRUE,
    inter_norm = TRUE
  )
  vdiffr::expect_doppelganger("PCA Plot (group scale)", p_pca_group)

  # Correlation plots
  p_corr_design <- expr_data$plot_corr(
    plot_scale = "group",
    intra_norm = TRUE,
    inter_norm = TRUE
  )
  vdiffr::expect_doppelganger("Correlation Plot (group scale)", p_corr_design)

  p_corr_batch <- expr_data$plot_corr(
    plot_scale = "batch",
    intra_norm = TRUE,
    inter_norm = TRUE
  )
  vdiffr::expect_doppelganger("Correlation Plot (batch scale)", p_corr_batch)

  # Hierarchical clustering plots
  p_hclust_design <- expr_data$plot_hclust(
    plot_scale = "design",
    intra_norm = TRUE,
    inter_norm = TRUE,
    clust_bar_var = c("batch", "group")
  )
  vdiffr::expect_doppelganger("Hclust Plot (design scale)", p_hclust_design)

  p_hclust_batch <- expr_data$plot_hclust(
    plot_scale = "batch",
    dist_method = "manhattan",
    hclust_method = "average",
    intra_norm = TRUE,
    inter_norm = TRUE
  )
  vdiffr::expect_doppelganger("Hclust Plot (batch scale)", p_hclust_batch)
})