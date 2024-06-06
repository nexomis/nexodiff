
root_path <- system.file("extdata/sim_data", package = "nexodiff")

design <- nexodiff::PairwiseDesignWithAnnotation$new(
    pairwise_design_file = paste(root_path, "samples.yml", sep = "/"),
    annotation_file = paste(root_path, "annotation.txt", sep = "/"),
    idmapping_file = paste(root_path, "id_mapping.tab", sep = "/"),
    src_dir = root_path
)

design_fake <- nexodiff::PairwiseDesignWithAnnotation$new(
    pairwise_design_file =
      paste(root_path, "samples_with_fake_path.yml", sep = "/"),
    annotation_file = paste(root_path, "annotation.txt", sep = "/"),
    idmapping_file = paste(root_path, "id_mapping.tab", sep = "/"),
    src_dir = root_path
)

testthat::test_that("test loading of design with fake path", {
  testthat::expect_false(
    is.null(design_fake)
  )
})

expr_data_tx <- nexodiff::ExprDataTranscript$new(design)

expr_data_tx$filter_and_set_selected_ids("human", filtered_var = "tax_name")

expr_data_tx$filter_and_set_selected_ids("mRNA", filtered_var = "type")

expr_data <- nexodiff::ExprDataGene$new(expr_data_tx)

expr_data$compute_and_set_intra_norm_fact("fpk")

expr_data$compute_and_set_inter_norm_fact(
  method = "median", norm_scale = "design")

g1 <- expr_data$plot_dist_per_sample()

g2 <- expr_data$plot_dist_per_sample(geoms = c("boxplot", "violin"))

g3 <- expr_data$plot_dist_per_sample(
  geoms = c("boxplot", "violin", "histo"),
  mean_fun = "nz.geometric")

g4 <- expr_data$plot_prcomp(plot_scale = "group")
g5 <- expr_data$plot_prcomp(plot_scale = "batch")
g6 <- expr_data$plot_corr(in_batch = "batch1", plot_scale = "batch")
g7 <- expr_data$plot_prcomp(in_batch = "batch1", tr_fn = (function(x) (x^2)),
  plot_scale = "design")
g8 <- expr_data$plot_complex(plot_type = "hclust",
  in_group = c("A", "B"), plot_scale = "design")
g9 <- expr_data$plot_prcomp(in_group = c("A", "B"), plot_scale = "batch")
g10 <- expr_data$plot_prcomp(in_group = c("A", "B"), plot_scale = "group")

pdf(NULL)

testthat::test_that("dist plots tests", {
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
  testthat::expect_error(
    g4,
    NA
  )
  testthat::expect_error(
    g5,
    NA
  )
  testthat::expect_error(
    g6,
    NA
  )
  testthat::expect_error(
    g7,
    NA
  )
  testthat::expect_error(
    g8,
    NA
  )
  testthat::expect_error(
    g9,
    NA
  )
  testthat::expect_error(
    g10,
    NA
  )
})

dev.off()