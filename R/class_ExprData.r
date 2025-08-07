#' @include utils.r
#' @include class_Annotation.r
#' @include class_PairwiseDesign.r
#' @include class_ExprData_utils_filter.r
#' @include class_ExprData_utils_summary.r
#' @include class_ExprData_utils_norm.r
#' @include class_ExprData_utils_plot.r
#' @include class_ExprData_utils_extract.r
#' @include class_ExprData_utils_plot_complex.r

NULL

#' An R6 class representing an RNA-Seq expression dataset.
#'
#' @description This class represents a RNA sequencing dataset containing the
#' following information:
#'
#' * expression values (raw count) of expressed tags (genes or transcripts) for
#' all samples.
#' * effective lengths of expressed tags (genes or transcripts) for all samples.
#' * normalization factors for expressed tags (genes or transcripts) for all
#' samples.
#' * design and annotation information; see \link{PairwiseDesignWith}
#' for details.
#'
#' By convention, the dataset is represented as a matrix:
#' * A column represents the expression values for one sample.
#' * A row represents the expression values for one expression tag.
#' Column names are the sample names, and row names are the expression tag
#' names.
#'
#' @param in_batch (optional) A character vector specifying the batch(es) to
#' include.
#' @param in_group (optional) A character vector specifying the group(s) to
#' include.
#' @param in_run (optional) A character vector specifying the run(s) to include.
#' @param intra_norm A boolean indicating whether to apply intra-sample
#' normalization.
#' @param log2_expr A boolean indicating whether to use the log2(x+2)
#' transformation for expression data.
#' @param group_variable The variable used to group, color, or fill samples.
#' @param show_all_text A boolean indicating whether to show all text on graphs.
#' @param inter_norm A boolean indicating whether to apply inter-sample
#' normalization.
#' @param include_ctrl A boolean indicating whether to include the control
#' group. This argument is only valid when specifying only one argument for
#' 'in_group' and 'in_batch'.
#' @param nrow integer, number of raw in the arranged graph layout
#' @param plot_type plot type
#' - "prcomp" PCA analysis
#' - "corr" correlation analysis
#' - "hclust"
#' @param plot_scale At which scale are the samples correlated. There are 3
#' possible option:
#' - "design": all samples are compared together
#' - "batch": all samples are compared together within batch
#' - "group": all samples are compared together within group including ctrl
#' Note that with inter_scale=TRUE there might be cases where it is
#' not compatible with the scale inf inter norm.
#' (see compute_and_set_inter_norm)
#' @param tr_fn data transformation function
#' @param ggplot_mod ggplot modifier that will be added to the graph that
#' are later arranged.
#' @param prcomp_args args to be used with prcomp. It is used for the
#' the following value of plot_type:
#' - prcomp
#' - hclust when dim_reduce is not null
#' @param prcomp_autoplot_args list of argument that are used with autoplot
#' to plot the prcomp object (see ggfortify autoplot.pca_common). Note that
#' the design is given as data therefore the variable defined in the design
#' can be used.
#' @param include_ctrl_at_group_scale whether to include the controls at
#' group scale
#' @param dist_method distance method, see philentropy::distance function 
#' for available methods
#' @param hclust_method hierarchical clustering method (see hclust function)
#' @param dim_reduce reduction before clustering (not yet implemented)
#' @param clust_bar_var list of variable to include as legend as a color bar
#' @param height_main height value for main graph when combined vertically
#' with others. Others heights will be set to 1. Note that there can be more
#' than 1 others graphs. Used in the following plot_type:
#' - hclust
#' @param width_main width value for main graph when combined horizontically
#' with others. Others widths will be set to 1. Note that there can be more
#' than 1 others graphs.
#' - hclust
#' @param tags vector of tag ids to use for the plot
#' @param tag_type name of the tag type to use:
#' - uniprot
#' - symbol
#' - type
#' - tax_id (taxonomy)
#' - tax_name (taxonomy)
#' if null, then the default id will be used (typed gene id at gene level
#' and transcript id at transcript level)
#' @param df_design_filter data frame to filter samples based on their name
#' and eventually their batch and/or group depending on the plot_scale.
#' - if plot_scale = "design", expected column is "sample_name"
#' - if plot_scale = "batch", expected columns are "sample_name" & "batch"
#' - if plot_scale = "group", expected columns are "sample_name" & "batch" &
#' "code"
#' Note that unexpected columns will be ignored.
#' Default is null meaning that there is no filtering.

ExprData <- R6::R6Class("ExprData", # nolint
  public = list(

    #' @title Show expression summary per tax_id, tax_name and rna type for
    #' selected expressed tags
    #' @description
    #' This function computes a summary per tax_id, tax_name and rna type for
    #' selected expressed tags, and returns a list of tables per feature.
    #' @param tr_fn A function to transform the expression values before
    #' summarizing. Default is identity function (no transformation).
    #' @param sum_fn A function to aggregate the transformed values. Default
    #' is "sum", but can be "sum_fn = function(x) sum(as.integer(x > 0))"
    #' to count non zero genes for example.
    #' @param intra_norm A boolean indicating whether to apply intra-sample
    #' normalization. Only used when type="norm".
    #' @param inter_norm A boolean indicating whether to apply inter-sample
    #' normalization. Only used when type="norm".
    #' @return A list of tables per feature
    #' @examples
    #' show_etags_summary()
    #' show_etags_summary(intra_norm=TRUE)
    show_etags_summary = function(
      in_batch = NULL, tr_fn = NULL, sum_fn = sum,
      intra_norm = FALSE, inter_norm = FALSE
    ) {

      data <- self$compute_norm(
        in_batch = in_batch, intra_norm = intra_norm, inter_norm = inter_norm
      )

      summarize_etags(
        private$annotation, private$at_gene_level, data, tr_fn, sum_fn
      )
    },

    #' @title Plot results from show_etags_summary.
    #' @param tr_fn A function to transform the expression values before
    #' summarizing. Default is identity function (no transformation).
    #' @param sum_fn A function to aggregate the transformed values. Default
    #' is "sum", but can be "sum_fn = function(x) sum(as.integer(x > 0))"
    #' to count non zero genes for example.
    #' @param intra_norm A boolean indicating whether to apply intra-sample
    #' normalization. Only used when type="norm".
    #' @param inter_norm A boolean indicating whether to apply inter-sample
    #' normalization. Only used when type="norm".
    #' @return A list of tables per feature
    #' @examples
    #' show_etags_summary()
    #' show_etags_summary(intra_norm=TRUE)
    plot_etags_summary = function(
      in_batch = NULL, tr_fn = NULL, sum_fn = sum,
      intra_norm = FALSE, inter_norm = FALSE
    ) {

      results <- self$show_etags_summary(in_batch, tr_fn, sum_fn, intra_norm,
                                         inter_norm)
      plot_summarized_etags(results, private$design)
    },

    #' @description
    #' Select expressed tags based on filtering on taxon or rna type.
    #' The results will be an intersect with the previous selection.
    #' You can reset the object if it's not desired.
    #' @param values A list of values used for filtering.
    #'   e.g c("mRNA", "transcript")
    #' @param filtered_var The variable being filtered. Possible values are
    #' - "tax_id" for filtering based on taxon id
    #' - "tax_name" for filtering based on taxon name
    #' - "type" for filtering based on rna type (default)
    #' @param filter_type The filtering method. Possible values are:
    #' - "keep" to keep only an expression tag if its value for the variable
    #' designated with `filtered_var` is in `values`
    #' - "excl" to keep only an expression tag if its value for the variable
    #' designated with `filtered_var` is not in `values`.
    filter_and_set_selected_ids = function(
      values, filtered_var = "type", filter_type = "keep"
    ) {
      private$selected_ids <- filter_etags(
        private$selected_ids, private$annotation, private$main_etag,
        values, filtered_var, filter_type
      )
    },

    #' @description
    #' Compute and set inter normalization factors
    #' @param norm_scale Scale to apply the normalization:
    #' - "design" for all samples in one go
    #' - "batch" for applying norm_factors per batch
    #' - "group" for applying norm factors per batch per group.
    #' The scale will have an impact on the reference selection if ref type is
    #' "ctrl" or "all". In consequences the same sample can be associated with
    #' different scaling factor if it appears in different batch/group.
    #' Therefore data normalized at the "batch" or "group" level can not be
    #' analyzed at the design level. (and group scale not at the batch level
    #' either). Indeed at group/batch levels the inter_norm can only be applied
    #' at group/batch level.
    #' @param norm_by level by which M and A vector are computed followed by.
    #' the scaling factor for normalization.
    #' @param ref_type Type of reference samples to use for normalization.
    #' - "all" for using all samples as reference
    #' - "ctrl" for using only control samples as reference
    #' - "specified" for using specific samples `ref_samples` as reference
    #' @param ref_samples Samples to use as reference for normalization.
    #' @param ref_mean Method for the mean of gene expression between samples of
    #' the reference.
    #' - "median"
    #' - "geometric"
    #' - "nz.geometric" geometric without zero
    #' - "mod.geometric" modified geometric with epsilon = 1e-05
    #' (see https://arxiv.org/abs/1806.06403)
    #' - "arithmetic"
    #' @param norm_mean Method for the mean of ratios between the target and the
    #' reference. See `ref_mean` for possible values.
    #' @param tgt_mean Method for the mean of gene expression between target
    #' samples prior to the calculation of M and A vector.
    compute_and_set_inter_norm_fact <- function( # nolint: object_length_linter.
      norm_scale = "group", norm_by = "sample",
      ref_type = "all", ref_samples = NULL,
      ref_mean = "mod.geometric", norm_mean = "median",
      tgt_mean = "mod.geometric", a_mean = "mod.geometric", a_trim_value = 1,
      m_trim_prop = 0
    ) {

      private$inter_norm_fact <- compute_inter_norm(
        private$raw, private$intra_norm_fact,
        ref_type, ref_samples, private$design, norm_scale, norm_by,
        ref_mean, norm_mean, tgt_mean, a_mean,
        a_trim_value, m_trim_prop
      )

      private$inter_norm_fact_opts <- list(
        norm_scale = norm_scale,
        norm_by = norm_by,
        ref_type = ref_type,
        ref_samples = ref_samples,
        ref_mean = ref_mean,
        norm_mean = norm_mean,
        m_trim_prop = m_trim_prop,
        tgt_mean = tgt_mean,
        a_mean = a_mean,
        a_trim_value = a_trim_value
      )

    },

    #' @description
    #' Set intra normalization factors.
    #' This function replaces any NaN values in the intra normalization factors
    #' with 1. NaN values may occur when the denominator in the calculations is
    #' 0. Replacing NaN values with 1 ensures that the normalization factors do
    #' not have missing or undefined values, which can cause issues in further
    #' analyses.
    #' @param method Normalization method:
    #'   * "none" for no normalization
    #'   * "fpm" for Fragment per million mapped
    #'   * "fpk" for Fragment per kilobase of transcript
    #'   * "tpm" for Transcript per million
    #'   * "fpkm" for Fragment per kilobase of transcript per million mapped
    compute_and_set_intra_norm_fact = function(method) {
      private$intra_norm_fact_method <- method
      private$intra_norm_fact <- compute_intra_norm_factor(
        method, private$selected_ids, self$get_raw(), self$get_len()
      )
    },

    #' @description
    #' Reset object, this affects:
    #' * expression tags selection (as initially)
    #' * samples selection (as initially)
    #' * normalization (No norm)
    #' @return interger matrix with raw counts for expressed tags
    reset = function() {
      private$selected_ids <- row.names(private$raw)
      self$compute_and_set_intra_norm_fact("none")
      private$inter_norm_fact <- NULL
      private$inter_norm_fact_opts <- NULL
    },

    #' @description
    #' Get raw count matrix
    #' @return interger matrix with raw counts for expressed tags
    get_raw = function() {
      private$raw[private$selected_ids, ]
    },

    #' @description
    #' Filter and get raw count matrix
    #' @return interger matrix with raw counts for expressed tags
    filter_and_get_raw = function(in_batch = NULL, in_group = NULL) {
      private$raw[private$selected_ids,
        private$design$extract_sample_names(in_batch, in_group)]
    },

    #' @description
    #' Get effective length matrix
    #' @return numeric matrix with effective length  for expressed tags
    get_len = function() {
      private$len[private$selected_ids, ]
    },

    #' @description
    #' filter and get effective length matrix
    #' @return numeric matrix with effective length  for expressed tags
    filter_and_get_len = function(in_batch = NULL, in_group = NULL) {
      private$len[private$selected_ids,
        private$design$extract_sample_names(in_batch, in_group)]
    },

    #' @description
    #' Get normalization factors matrix
    #' @return Numeric matrix with normalization factors for expressed tags
    compute_norm_fact = function(in_batch = NULL, in_group = NULL,
      inter_norm = FALSE, intra_norm = TRUE) {
      compute_norm_fact_helper(
        private$selected_ids, private$design, private$intra_norm_fact,
        private$inter_norm_fact, private$inter_norm_fact_opts,
        in_batch, in_group, inter_norm, intra_norm
      )
    },

    #' @description
    #' Get intra normalization method
    #' @return method name
    get_intra_norm_fact_method = function() {
      private$intra_norm_fact_method
    },

    #' @description
    #' Get inter normalization methods
    #' @return lists with all normalizations options
    get_inter_norm_fact_opts = function() {
      private$inter_norm_fact_opts
    },

    #' @description
    #' Compute the normalized expression data from raw and norm matrices.
    #'
    #' **IMPORTANT**: Default is no normalization. Please use the function below
    #' beforehand.
    #' - "compute_and_set_intra_norm_fact"
    #' - "compute_and_set_inter_norm_fact"
    #' Note that intra normalization is always applied, however inter norm is
    #' optional. Depending on the `norm_scale` used for inter normalization,
    #' the `in_batch` argument might be mandatory.
    #' @return numeric matrix with norm counts for expressed tags
    compute_norm = function(
      in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = FALSE
    ) {
      if (intra_norm) {
        self$filter_and_get_raw(in_batch, in_group) *
          self$compute_norm_fact(in_batch, in_group, inter_norm = inter_norm)
      } else {
        self$filter_and_get_raw(in_batch, in_group)
      }
    },

    #' @description
    #' Are the expression tags lengths fixed ?
    #' @return numeric matrix with raw counts for expressed tags
    is_with_fixed_length = function() {
      private$with_fixed_length
    },

    #' @description
    #' Get the main expression tag id. The tag id corresponds to the unique id
    #' linked with a expression point. It can be tx_id (for transcript id) or
    #' "tgid" for "typed gene id". Why typed gene ID ? because some genes have
    #' transcripts from different types. As it would not be biologically
    #' relevant to aggregate different typse during the sum, a new id is create
    #' that is gene and type specific.
    #' @return chararacter representing the main espression tag id
    get_main_etag = function() {
      private$main_etag
    },

    #' @description
    #' Whether the expression data is at gene level or not.
    #' @return T or F
    is_at_gene_level = function() {
      private$at_gene_level
    },

    #' @description
    #' Extract the PairwiseDesign instance that were used to
    #' build the expression data object.
    #' @return `PairwiseDesign` object
    get_design = function() {
      private$design
    },


    #' @description
    #' Extract the Annotation instance that were used to
    #' build the expression data object.
    #' @return `Annotation` object
    get_annotation = function() {
      private$annotation
    },

    #' @description
    #' Plot the distribution of count or normalized values per RNA type for each
    #' sample.
    #' @param geoms Character vector, types of plot to draw, can be a
    #' combination of:
    #' - "boxplot": draw a box plot (fast)
    #' - "violin": draw a violin plot (slow)
    #' - "histo": draw an hiostogram (slow)
    #' @param tr_fn function to transform expression values before plot.
    #' @param mean_fun (Optional) Function to plot the mean per sample.
    #' Can be one of the following:
    #' * "median"
    #' * "geometric"
    #' * "nz.geometric" geometric without zero
    #' * "mod.geometric" modified geometric with epsilon = 1e-05
    #' (see https://arxiv.org/abs/1806.06403)
    #' * "arithmetic"
    #' @return ggplot2 graph
    plot_dist_per_sample = function(
      intra_norm = TRUE, inter_norm = TRUE, geoms = c("boxplot"),
      tr_fn = (function(x) log2(x + 2) - 1), mean_fun = NULL
    ) {
      data_design_batch_group <- dplyr::distinct(
        private$design$get_pairwise_design()[, c("batch", "group")]
      )


      data <- purrr::map2_dfr(
        data_design_batch_group$batch,
        data_design_batch_group$group,
        function(in_batch, in_group) {
          as.tibble(self$compute_norm(in_batch, in_group, inter_norm)) %>%
            dplyr::mutate(batch = in_batch, batch = in_batch) %>%
            tidyr::pivot_longer(
              ! tidyselect::all_of(c("batch", "group")),
              names_to = "sample", values_to = "value"
            ) %>%
            dplyr::mutate(value = tr_fn(value))
        }
      )

      plot_dist_per_sample_helper(
        data, private$design, geoms, mean_fun, log2_expr
      )
    },


    #' @description
    #' Get formatted data for a group with its matching control within a batch
    #' @return list with selected attributes:
    #' - "raw" : matrix of raw values for test and control samples specified as
    #' column names
    #' - "len" : matrix of tag lenght values for test and control samples
    #' specified as column names
    #' - "intra_norm_fact" : matrix intra normalization factors
    #' - "inter_norm_fact" : vector of inter normalization factors
    #' - "norm_fact" : vector of inter normalization factors
    #' - "test_samples" : vector with test sample names
    #' - "ctrl_samples" : vector with control sample names
    #' - "design_table" : design table only for those samples
    #' - "ctrl_group" : the control group
    extract_pairwise_data_with_design = function(in_batch, in_group) {
      group_per_batches <-
        private$design$list_groups_per_batches(include_ctrl = TRUE)

      assert_that(length(in_batch) == 1 && length(in_group) == 1)
      assert_that(in_batch %in% names(group_per_batches))
      assert_that(in_group %in% group_per_batches[in_batch][[1]])

      results <- list()

      # test if control group is the same... Sometimes we use this to get the
      # for a control group also.
      results$test_samples <-
        private$design$extract_sample_names(in_batch, in_group)

      # find batch control
      results$ctrl_group <-
        private$design$find_control_group_per_batches()[in_batch]
      results$ctrl_samples <-
        private$design$extract_sample_names(in_batch, results$ctrl_group)
      all_groups <- c(in_group, results$ctrl_group)
      results$all_samples <-
        unique(c(results$test_samples, results$ctrl_samples))

      # norm is not computed because if the table is big, it's too long

      results$raw <- self$filter_and_get_raw(in_batch, all_groups)[,
        results$all_samples
      ]
      results$len <- self$filter_and_get_len(in_batch, all_groups)[,
        results$all_samples
      ]

      results$intra_norm_fact <- self$compute_norm_fact(
        in_batch, all_groups, inter_norm = FALSE, intra_norm = TRUE
      )[, results$all_samples]

      results$inter_norm_fact <- self$compute_norm_fact(
        in_batch, all_groups, inter_norm = TRUE, intra_norm = FALSE
      )[results$all_samples]

      results$design_table <-
        private$design$get_pairwise_design(in_batch, all_groups)

      row.names(results$design_table) <- results$design_table$sample
      results$design_table <- results$design_table[results$all_samples, ]
      results
    },

    #' @description
    #' Clean version suffixes from expression tag IDs.
    #'
    #' This method applies the `remove_id_version_suffix` utility to the
    #' row names of internal data matrices (`raw`, `len`, `intra_norm_fact`)
    #' and to `selected_ids`.
    #'
    #' @param suffix_pattern (character) The regex pattern for the version
    #'   suffix. Default is `"-\\d+$"`. If `NULL`, `FALSE`, or an empty string,
    #'   no cleaning is performed.
    clean_id_versions = function(suffix_pattern = "\\.\\d+$") {
      rownames(private$raw) <-
        remove_id_version_suffix(rownames(private$raw), suffix_pattern)
      rownames(private$len) <-
        remove_id_version_suffix(rownames(private$len), suffix_pattern)
      rownames(private$intra_norm_fact) <-
        remove_id_version_suffix(
          rownames(private$intra_norm_fact), suffix_pattern
        )
      private$selected_ids <-
        remove_id_version_suffix(private$selected_ids, suffix_pattern)
    },

    #' @description
    #' Function to draw PCA plots
    #' @param in_batch (optional) A character vector specifying the batch(es) to include.
    #' @param in_group (optional) A character vector specifying the group(s) to include.
    #' @param intra_norm A boolean indicating whether to apply intra-sample normalization.
    #' @param inter_norm A boolean indicating whether to apply inter-sample normalization.
    #' @param tr_fn data transformation function
    #' @param plot_scale At which scale are the samples correlated.
    #' @param ggplot_mod ggplot modifier that will be added to the graph.
    #' @param color_palette name of the palette of color to use
    #' @param prcomp_args args to be used with prcomp.
    #' @param prcomp_autoplot_args list of argument that are used with autoplot.
    #' @param include_ctrl_at_group_scale whether to include the controls at group scale
    #' @param tags vector of tag ids to use for the plot
    #' @param tag_type name of the tag type to use
    #' @param df_design_filter data frame to filter samples
    #' @return ggplot2 graph
    plot_prcomp = function(
      in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL, color_palette = "BrBg",
      prcomp_args = list(), prcomp_autoplot_args = list(),
      include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL, df_design_filter = NULL
    ) {
      private$plot_complex(
        "prcomp", in_batch, in_group, intra_norm, inter_norm, tr_fn,
        plot_scale, ggplot_mod, color_palette,
        prcomp_args, prcomp_autoplot_args,
        dist_method = NULL, hclust_method = NULL, dim_reduce = NULL,
        clust_bar_var = NULL, height_main = NULL, width_main = NULL,
        include_ctrl_at_group_scale,
        tags, tag_type, df_design_filter
      )
    },

    #' @description
    #' Function to draw correlation plots
    #' @param in_batch (optional) A character vector specifying the batch(es) to include.
    #' @param in_group (optional) A character vector specifying the group(s) to include.
    #' @param intra_norm A boolean indicating whether to apply intra-sample normalization.
    #' @param inter_norm A boolean indicating whether to apply inter-sample normalization.
    #' @param tr_fn data transformation function
    #' @param plot_scale At which scale are the samples correlated.
    #' @param ggplot_mod ggplot modifier that will be added to the graph.
    #' @param include_ctrl_at_group_scale whether to include the controls at group scale
    #' @param tags vector of tag ids to use for the plot
    #' @param tag_type name of the tag type to use
    #' @param df_design_filter data frame to filter samples
    #' @return ggplot2 graph
    plot_corr = function(
      in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL,
      include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL, df_design_filter = NULL
    ) {
      private$plot_complex(
        "corr", in_batch, in_group, intra_norm, inter_norm, tr_fn,
        plot_scale, ggplot_mod, color_palette = NULL,
        prcomp_args = NULL, prcomp_autoplot_args = NULL,
        dist_method = NULL, hclust_method = NULL, dim_reduce = NULL,
        clust_bar_var = NULL, height_main = NULL, width_main = NULL,
        include_ctrl_at_group_scale,
        tags, tag_type, df_design_filter
      )
    },

    #' @description
    #' Function to draw hierarchical clustering plots
    #' @param in_batch (optional) A character vector specifying the batch(es) to include.
    #' @param in_group (optional) A character vector specifying the group(s) to include.
    #' @param intra_norm A boolean indicating whether to apply intra-sample normalization.
    #' @param inter_norm A boolean indicating whether to apply inter-sample normalization.
    #' @param tr_fn data transformation function
    #' @param plot_scale At which scale are the samples correlated.
    #' @param ggplot_mod ggplot modifier that will be added to the graph.
    #' @param color_palette name of the palette of color to use
    #' @param dist_method distance method
    #' @param hclust_method hierarchical clustering method
    #' @param dim_reduce reduction before clustering
    #' @param clust_bar_var list of variable to include as legend as a color bar
    #' @param height_main height value for main graph
    #' @param width_main width value for main graph
    #' @param prcomp_args args to be used with prcomp.
    #' @param include_ctrl_at_group_scale whether to include the controls at group scale
    #' @param tags vector of tag ids to use for the plot
    #' @param tag_type name of the tag type to use
    #' @param df_design_filter data frame to filter samples
    #' @return ggplot2 graph
    plot_hclust = function(
      in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL, color_palette = "BrBg",
      dist_method = "euclidean", hclust_method = "ward.D2", dim_reduce = NULL,
      clust_bar_var = c(), height_main = 10, width_main = 4,
      prcomp_args = list(),
      include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL, df_design_filter = NULL
    ) {
      private$plot_complex(
        "hclust", in_batch, in_group, intra_norm, inter_norm, tr_fn,
        plot_scale, ggplot_mod, color_palette,
        prcomp_args = prcomp_args,
        prcomp_autoplot_args = NULL,
        dist_method, hclust_method, dim_reduce,
        clust_bar_var, height_main, width_main,
        include_ctrl_at_group_scale,
        tags, tag_type, df_design_filter
      )
    }
  ),
  private = list(
    #A boolean indicating whether the expression data is at gene level rather
    #than transcripts-level.
    at_gene_level = NULL,
    #A integer matrix with raw counts per genes/tx per sample. gene or tx id
    #as rownames and sample as colnames.
    raw = NULL,
    #A integer matrix with (effective) length of genes/tx per sample. gene or
    # tx id as rownames and sample as colnames.
    len = NULL,
    #A list per batch with a numeric matrix with normalization factors per
    #genes/tx per sample. gene or tx id as rownames and sample as colnames.
    intra_norm_fact = NULL,
    inter_norm_fact = NULL,
    #A PairwiseDesign object that describes the samples
    design = NULL,
    annotation = NULL,
    # Normalization method for norm factors
    intra_norm_fact_method = NULL,
    inter_norm_fact_opts = NULL,
    # etags IDs selected based on filters
    selected_ids = NULL,
    # sample names selected based on filters
    selected_samples = NULL,
    # whether the expression tag lengths are fixed or not
    with_fixed_length = NULL,
    # main expression tag which is the row-nams of data frame raw, len and
    # norm_fact
    main_etag = NULL,
    plot_complex = function(
      plot_type, in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL, color_palette = "BrBg",
      prcomp_args = list(), prcomp_autoplot_args = list(),
      dist_method = "euclidean", hclust_method = "ward.D2", dim_reduce = NULL,
      clust_bar_var = c(), height_main = 10, width_main = 4,
      include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL, df_design_filter = NULL
    ) {
      # Check arguments
      plot_complex_check(
        plot_scale,inter_norm, private$inter_norm_fact_opts, in_batch
      )

      # Prepare tag selection
      selected_tag_ids <- prepare_tag_selection(
        private$main_etag, private$annotation,
        private$selected_ids, tags, tag_type
      )

      # Prepare data based on plot scale
      data_prep <- prepare_complex_plot_data(
        selected_tag_ids, private$design, NULL, in_batch, in_group,
        df_design_filter, plot_scale, include_ctrl_at_group_scale, tr_fn
      )

      data_design <- data_prep$data_design
      df_design_filter <- data_prep$df_design_filter

      if (plot_scale == "group") {
        graphs <- purrr::map2(
          data_design$batch,
          data_design$group,
          function(x, y) {
            in_data <- self$compute_norm(x, y, inter_norm)
            if (include_ctrl_at_group_scale) {
              ctrl_group <- private$design$find_control_group_per_batches()[x]
              ctrl_data <- self$compute_norm(x, ctrl_group, inter_norm)
              in_data <- cbind(in_data, ctrl_data)
            }
            selected_samples <- intersect(
              unique(dplyr::filter(df_design_filter,
                .data$batch == .env$x & .data$group == .env$y)$sample),
              colnames(in_data)
            )
            in_data <- in_data[selected_tag_ids, selected_samples]
            in_title <- paste(
              private$design$get_b_labels()[x],
              private$design$get_g_labels()[y],
              sep = ": "
            )

            # Dispatch to the appropriate plotting function based on plot_type
            if (plot_type == "prcomp") {
              create_prcomp_plot(
                in_data, in_title,
                private$design$get_pairwise_design() %>%
                  dplyr::filter(.data$sample %in% colnames(in_data)) %>%
                  dplyr::distinct(),
                prcomp_args, prcomp_autoplot_args, ggplot_mod, color_palette
              )
            } else if (plot_type == "corr") {
              create_corr_plot(in_data, in_title, ggplot_mod)
            } else if (plot_type == "hclust") {
              create_hclust_plot(
                in_data, in_title,
                private$design$get_pairwise_design() %>%
                  dplyr::filter(.data$sample %in% colnames(in_data)) %>%
                  dplyr::distinct(),
                dist_method, hclust_method, dim_reduce, clust_bar_var,
                height_main, width_main, color_palette, prcomp_args, ggplot_mod
              )
            }
          }
        )
      } else if (plot_scale == "batch") {
        graphs <- purrr::map(
          unique(data_design$batch),
          function(x) {
            if (intra_norm) {
              data <- self$compute_norm(x, in_group, inter_norm)
            } else {
              data <- self$filter_and_get_raw(x, in_group)
            }
            selected_samples <- intersect(
              unique(dplyr::filter(df_design_filter,
                .data$batch == .env$x)$sample),
              colnames(data)
            )
            data <- data[selected_tag_ids, selected_samples]
            in_title <- private$design$get_b_labels()[x]

            # Adjust autoplot args for batch scale
            prcomp_autoplot_args_mod <- prcomp_autoplot_args
            n_args <- names(prcomp_autoplot_args_mod)
            if (! "colour" %in% n_args) {
              prcomp_autoplot_args_mod$colour <- "g_label"
            }

            # Dispatch to the appropriate plotting function based on plot_type
            if (plot_type == "prcomp") {
              create_prcomp_plot(
                data, in_title,
                private$design$get_pairwise_design() %>%
                  dplyr::filter(.data$sample %in% colnames(data)) %>%
                  dplyr::distinct(),
                prcomp_args, prcomp_autoplot_args_mod, ggplot_mod, color_palette
              )
            } else if (plot_type == "corr") {
              create_corr_plot(data, in_title, ggplot_mod)
            } else if (plot_type == "hclust") {
              create_hclust_plot(
                data, in_title,
                private$design$get_pairwise_design() %>%
                  dplyr::filter(.data$sample %in% colnames(data)) %>%
                  dplyr::distinct(),
                dist_method, hclust_method, dim_reduce, clust_bar_var,
                height_main, width_main, color_palette, prcomp_args, ggplot_mod
              )
            }
          }
        )
      } else if (plot_scale == "design") {
        if (intra_norm) {
          data <- self$compute_norm(in_batch, in_group, inter_norm)
        } else {
          data <- self$filter_and_get_raw(in_batch, in_group)
        }
        selected_samples <- intersect(
          unique(df_design_filter$sample),
          colnames(data)
        )
        data <- data[selected_tag_ids, selected_samples]
        in_title <- "Design"

        # Adjust autoplot args for design scale
        prcomp_autoplot_args_mod <- prcomp_autoplot_args
        n_args <- names(prcomp_autoplot_args_mod)
        if (! "colour" %in% n_args) {
          prcomp_autoplot_args_mod$colour <- "g_label"
        }
        if (! "shape" %in% n_args) {
          prcomp_autoplot_args_mod$shape <- "b_label"
        }

        # Dispatch to the appropriate plotting function based on plot_type
        if (plot_type == "prcomp") {
          graphs <- list(create_prcomp_plot(
            data, in_title,
            private$design$get_pairwise_design() %>%
              dplyr::filter(.data$sample %in% colnames(data)) %>%
              dplyr::distinct(),
            prcomp_args, prcomp_autoplot_args_mod, ggplot_mod, color_palette
          ))
        } else if (plot_type == "corr") {
          graphs <- list(create_corr_plot(data, in_title, ggplot_mod))
        } else if (plot_type == "hclust") {
          graphs <- list(create_hclust_plot(
            data, in_title,
            private$design$get_pairwise_design() %>%
              dplyr::filter(.data$sample %in% colnames(data)) %>%
              dplyr::distinct(),
            dist_method, hclust_method, dim_reduce, clust_bar_var,
            height_main, width_main, color_palette, prcomp_args, ggplot_mod
          ))
        }
      }

      arrange_complex_plots(graphs)
    }
  )
)


#' Class representing an expression dataset at gene-level
#'
#' @description
#' Extension of expression dataset at gene-level.
#'
#' @export

ExprDataGene <- R6::R6Class("ExprDataGene", # nolint
  inherit = ExprData,
  public = list(

    #' @description
    #' Initialize `ExprDataGene` object.
    #'
    #' The `ExprDataGene` object is initiated from a `ExprDatatranscript` object
    #'
    #' Raw count are summarized per genes.
    #'
    #' Lengths are computed using a weighted mean. Weights are raw count divided
    #' by length
    #'
    #' @param expr_data_transcript Expression data at transcript level.
    #' (`ExprDataTranscript` object)
    #' @return A new `ExprDataGene` object.
    initialize = function(expr_data_transcript) {

      private$design <- expr_data_transcript$get_design()
      private$annotation <- expr_data_transcript$get_annotation()
      private$with_fixed_length <- expr_data_transcript$is_with_fixed_length()
      private$intra_norm_fact_method <- "none"
      private$inter_norm_fact_opts <- NULL
      private$inter_norm_fact <- NULL
      private$main_etag <- "tgid"
      private$at_gene_level <- TRUE
      # Summarize raw counts

      txid2tgid <- private$annotation$generate_translate_dict(
        expr_data_transcript$get_main_etag(), private$main_etag
      )

      private$raw <- as.data.frame(expr_data_transcript$get_raw())
      private$raw$gene <-
        txid2tgid[row.names(private$raw)]
      private$raw <- as.data.frame(private$raw %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(dplyr::across(tidyselect::everything(), sum))
      )
      row.names(private$raw) <- private$raw$gene
      private$raw$gene <- NULL
      private$raw <- as.matrix(private$raw)

      # Initialize selected etags and norm_fact
      private$selected_ids <- row.names(private$raw)

      private$intra_norm_fact <- matrix(rep(1, length(private$raw)),
        nrow = nrow(private$raw),
        ncol = ncol(private$raw),
        dimnames = list(
          row.names(private$raw),
          colnames(private$raw)
        )
      )

      # Initialize lengths
      if (expr_data_transcript$is_with_fixed_length()) {
        private$len <- private$intra_norm_fact # all lengths are 1
      } else {
        # Ponderation of lengths with normalized values (raw/length)
        norm <- expr_data_transcript$get_raw() / expr_data_transcript$get_len()
        private$len <- as.data.frame(expr_data_transcript$get_len() * norm)

        # Summarize normalized values for later depondaration
        norm <- as.data.frame(norm)
        norm$gene <- txid2tgid[row.names(norm)]
        g_norm <- as.data.frame(norm %>%
          dplyr::group_by(gene) %>%
          dplyr::summarise(dplyr::across(tidyselect::everything(), sum))
        )
        row.names(g_norm) <- g_norm$gene
        g_norm$gene <- NULL
        g_norm <- as.matrix(g_norm)
        # WARNING G_NORM CAN BE ZERO !!! DIVISION BY ZERO

        # Summarize ponderated lengths
        private$len$gene <-
          txid2tgid[row.names(private$len)]
        private$len <- as.data.frame(private$len %>%
          dplyr::group_by(gene) %>%
          dplyr::summarise(dplyr::across(tidyselect::everything(), sum))
        )
        row.names(private$len) <- private$len$gene
        private$len$gene <- NULL
        # deponderate
        private$len <- as.matrix(private$len) / g_norm
        # replace length of 0 by 1 to avoid dividing by 0 later
        private$len[private$len == 0] <- 1
        # remove nan values
        private$len[is.nan(private$len)] <- 1
      }
    }
  )
)


#' Class representing an expression dataset at transcript-level
#'
#' @description
#' Extension of expression dataset at transcript-level
#' @export

ExprDataTranscript <- R6::R6Class("ExprDataTranscript", # nolint
  inherit = ExprData,
  public = list(
    #' @description
    #' Initialize `ExprDataTranscript` object
    #' @param design A pairwise design object \link{PairwiseDesign}
    #' @param annotation An annotation object \link{Annotation}
    #' @param with_fixed_length A boolean indicating whether the NGS library is
    #' of a fixed length like 3 prime or not
    #' @param log_level Level of logging (see logging package). Default = WARN
    #' @param format Format of files. Default=kallisto
    #' @param suffix_pattern See clean_txid_versions method
    #' Default = txid
    #' @return A new `ExprDataTranscript` object.
    initialize = function(
      design, annotation, with_fixed_length = FALSE, log_level = "WARN",
      format = "kallisto", suffix_pattern = "\\.\\d+$"
    ) {
      logging::basicConfig(log_level)

      private$design <- design$clone()
      private$annotation <- annotation$clone()
      private$at_gene_level <- FALSE
      private$intra_norm_fact_method <- "none"
      private$inter_norm_fact_opts <- NULL
      private$inter_norm_fact <- NULL
      private$with_fixed_length <- with_fixed_length
      private$main_etag <- "txid"
      # todo try to detect if fixed length ...

      # Import data

      logging::loginfo("Import transcript abundance files (.h5)")

      files <- private$design$build_file_paths()

      if (!all(file.exists(files))) {
        logging::logerror("The following files are missing")
        logging::logerror(files[!file.exists(files)])
        stop()
      }

      if (format == "kallisto") {

        private$raw <- NULL
        private$selected_ids <- NULL

        for (sample_name in names(files)) {

          fpath <- files[[sample_name]]
          ids <- as.character(rhdf5::h5read(fpath, "aux/ids"))

          # Check for duplicates in ids
          assert_that(!any(duplicated(ids)))

          if (is.null(private$selected_ids)) {
            private$selected_ids <- ids
            private$raw <- matrix(1, nrow = length(ids),
                                  ncol = length(files))
            rownames(private$raw) <- ids
            colnames(private$raw) <- names(files)
            private$len <- private$raw
            slug_ids <- digest::digest(sort(ids))
          } else {
            assert_that(digest::digest(sort(ids)) == slug_ids)
          }

          # Read counts and verify dimensions before assignment
          counts <- rhdf5::h5read(fpath, "est_counts")
          assert_that(length(counts) == length(ids))

          private$raw[ids, sample_name] <- counts

          if (! with_fixed_length) {
            lengths <- rhdf5::h5read(fpath, "aux/eff_lengths")
            assert_that(length(lengths) == length(ids))
            private$len[ids, sample_name] <- lengths
          }
        }

      } else {
        logging::logerror("wrong value for `format` argument")
        stop()
      }

      private$intra_norm_fact <- matrix(rep(1, length(private$raw)),
        nrow = nrow(private$raw),
        ncol = ncol(private$raw),
        dimnames = list(
          row.names(private$raw),
          colnames(private$raw)
        )
      )

      private$len[is.nan(private$len)] <- 1

      self$clean_id_versions(suffix_pattern = suffix_pattern)
    }
  )
)
