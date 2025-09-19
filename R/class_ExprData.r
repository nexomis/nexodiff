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

#' @title An R6 class representing an RNA-Seq expression dataset.
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

    #' Show expression summary per tax_id, tax_name and rna type for
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

    #' Plot results from show_etags_summary.
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

    #' Compute and set inter normalization factors
    #' @description
    #' Compute and set inter normalization factors
    #' @description
    #' This function computes scaling factors using a ratio-based
    #' method to correct for systematic technical biases between samples, such
    #' as differences in sequencing depth or library composition. By default, it
    #' performs Median Ratio Normalization (MRN), as `norm_mean` is `"median"`
    #' and M-value trimming is disabled (`m_trim_prop = 0`). It can be
    #' configured to perform TMM-style normalization by setting `m_trim_prop`
    #' and norm_mean="geometric"` or `norm_mean="mod.geometric"`. Note that the
    #' weigthed mean for TMM is not available yet.
    #'
    #' The normalization process is as follows:
    #' 1.  **Reference Selection**: A reference expression profile is created
    #' from samples chosen via `ref_type` and `ref_samples`. Gene expression
    #' values are averaged using the `ref_mean` method.
    #' 2.  **Target Processing**: For each target sample or group (`norm_by`), a
    #' target expression profile is created using `tgt_mean`.
    #' 3.  **M and A Vector Calculation**:
    #'     - The **A-vector** (average expression) is used to filter out
    #'     low-expression genes, which can be noisy.
    #'     - The **M-vector** (expression ratio: target/reference) is used to
    #'     derive the normalization factor from its central tendency.
    #' 4.  **Trimming Process**: A sequential trimming process refines the set
    #' of genes used for normalization:
    #'     - **A-vector trimming**: Genes with average expression below
    #'     `a_trim_value` are removed.
    #'     - **M-vector trimming**: A proportion (`m_trim_prop`) of genes with
    #'     the most extreme expression ratios are removed. This is disabled by
    #'     default (`m_trim_prop = 0`).
    #'     - **Extreme ratio trimming**: If `trim_extreme` is `TRUE`, genes with
    #'     infinite or zero ratios are removed. This is advised when `norm_mean`
    #'     is not `"median"` or when `m_trim_prop` is too small to remove these
    #'     values.
    #' 5.  **Factor Calculation**: The final normalization factor is computed by
    #' applying `norm_mean` to the ratio (M-vector) of the remaining genes.
    #' 6.  **Normalization Scope (`norm_scale`)**: This parameter determines the
    #' scope and reference for normalization.
    #'     - **`design`**: Normalization is performed across all samples at once
    #'     This is the generally recommended approach as it makes all samples
    #'     comparable across the entire experiment.
    #'     - **`batch`** or **`group`**: Normalization is performed
    #'     independently within each context. This can mitigate biases if a
    #'     particular batch or group is compositionally very different. However,
    #'     it renders data comparable only *within* that same context
    #'     (e.g., within a batch), not across the entire design.
    #'
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
    #' @param ref_mean Method for averaging gene expression across reference
    #'   samples. It is recommended to use the same method for `ref_mean`,
    #'   `tgt_mean`, and `a_mean`.
    #'   - **`"arithmetic"`**: Recommended when using library size-scaled
    #'   intra-normalization (e.g., "tpm", "fpkm", "fpm"). This preserves the
    #'   sum of proportions.
    #'   - **`"geometric"`**: Recommended when intra-normalization does not
    #'   account for library size (e.g., "none", "fpk"). The library size effect
    #'   is then absorbed into the inter-sample normalization factor. This
    #'   approach, with "none" intra-normalization, mirrors the default
    #'   behavior in DESeq2. Note that this can be problematic if normalization
    #'   is performed per group (`norm_by = "group"`).
    #'   - Other options: `"median"`, `"nz.geometric"`, `"mod.geometric"`.
    #' @param norm_mean Method for the mean of ratios between the target and the
    #' reference. See `ref_mean` for possible values.
    #' @param tgt_mean Method for averaging gene expression across target
    #'   samples. See `ref_mean` for recommended usage and available options.
    #' @param a_mean Method for averaging gene expression between target and
    #'   reference samples to calculate the A-vector. See `ref_mean` for
    #'   recommended usage and available options.
    #' @param a_trim_value Value for trimming the A vector.
    #' @param m_trim_prop Proportion of values to trim from the M vector.
    #' @param trim_extreme Whether to trim extreme values from the M vector,
    #' before calculating the mean. If the mean function `norm_mean` is not
    #' "median", extreme values will causes issues with the calculation of the
    #' mean.
    compute_and_set_inter_norm_fact = function( # nolint: object_length_linter.
      norm_scale = "design", norm_by = "sample",
      ref_type = "all", ref_samples = NULL,
      ref_mean = "arithmetic", norm_mean = "median",
      tgt_mean = "arithmetic", a_mean = "arithmetic", a_trim_value = 1,
      m_trim_prop = 0, trim_extreme = FALSE
    ) {

      private$inter_norm_fact <- compute_inter_norm(
        self$get_raw(), private$intra_norm_fact,
        ref_type, ref_samples, private$design, norm_scale, norm_by,
        ref_mean, norm_mean, tgt_mean, a_mean,
        a_trim_value, m_trim_prop, trim_extreme
      )

      private$inter_norm_fact_opts <- list(
        norm_scale = norm_scale,
        norm_by = norm_by,
        ref_type = ref_type,
        ref_samples = ref_samples,
        ref_mean = ref_mean,
        norm_mean = norm_mean,
        tgt_mean = tgt_mean,
        a_mean = a_mean,
        a_trim_value = a_trim_value,
        m_trim_prop = m_trim_prop,
        trim_extreme = trim_extreme
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
    #' @param include_ctrl when in_group is specified, will include batch ctrl
    #' @return interger matrix with raw counts for expressed tags
    filter_and_get_raw = function(in_batch = NULL, in_group = NULL,
      include_ctrl = FALSE
    ) {
      ids <- private$design$extract_sample_names(in_batch, in_group)
      if (include_ctrl) {
        assert_that(! is.null(in_group))
        ctrl_group <- private$design$find_control_group_per_batches()[in_batch]
        ids <- unique(c(
          private$design$extract_sample_names(in_batch, ctrl_group), ids
        ))
      }
      private$raw[
        private$selected_ids,
        ids
      ]
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
    #' @param include_ctrl if TRUE and if `in_batch` and `in_group` are
    #' specified, control samples from the specified batch will be included.
    #' @param rescale_inter_norm Rescale inter-sample normalization factors so
    #' that their geometric mean is 1. Default is TRUE.
    #' @return Numeric matrix with normalization factors for expressed tags
    compute_norm_fact = function(
      in_batch = NULL, in_group = NULL, inter_norm = FALSE, intra_norm = TRUE,
      include_ctrl = FALSE, rescale_inter_norm = TRUE
    ) {
      compute_norm_fact_helper(
        private$selected_ids, private$design, private$intra_norm_fact,
        private$inter_norm_fact, private$inter_norm_fact_opts,
        in_batch, in_group, inter_norm, intra_norm, include_ctrl,
        rescale_inter_norm
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
    #' @param include_ctrl include control samples in the normalization
    #' @param rescale_inter_norm Rescale inter-sample normalization factors so
    #' that their geometric mean is 1. Default is TRUE.
    #' @return numeric matrix with norm counts for expressed tags
    compute_norm = function(
      in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = FALSE, include_ctrl = FALSE,
      rescale_inter_norm = TRUE
    ) {
      if (intra_norm) {
        intra <- self$filter_and_get_raw(in_batch, in_group, include_ctrl)
        inter <-
          self$compute_norm_fact(
            in_batch, in_group,
            inter_norm = inter_norm, intra_norm = intra_norm,
            include_ctrl = include_ctrl,
            rescale_inter_norm = rescale_inter_norm
          )
        assert_that(all(colnames(inter) %in% colnames(intra)))
        intra * inter[, colnames(intra)]
      } else {
        assert_that(! inter_norm,
          msg = "`intra_norm` cannot be FALSE if `inter_norm` is TRUE"
        )
        self$filter_and_get_raw(in_batch, in_group, include_ctrl)
      }
    },

    #' @description
    #' Are the expression tags lengths fixed ?
    #' @return numeric matrix with raw counts for expressed tags
    is_with_fixed_length = function() {
      private$with_fixed_length
    },

    #' @description
    #' Get the intra normalization factor
    #' @return numeric matrix with intra normalization factor
    get_intra_norm_fact = function(in_batch = NULL, in_group = NULL) {
      if (is.null(in_batch) && is.null(in_group)) {
        private$intra_norm_fact[private$selected_ids, ]
      } else {
        assert_that(! is.null(in_batch))
        private$intra_norm_fact[
          private$selected_ids,
          private$design$extract_sample_names(in_batch, in_group)
        ]
      }
    },
    #' @description
    #' Get the inter normalization factor
    #' @param rescale Rescale the factors so that their geometric mean is 1.
    get_inter_norm_fact = function(in_batch = NULL, in_group = NULL,
                                   rescale = TRUE) {
      assert_that(!is.null(private$inter_norm_fact))
      norm_scale <- private$inter_norm_fact_opts$norm_scale
      assert_that(! (norm_scale != "design" && is.null(in_batch)))
      assert_that(! (norm_scale == "group" && is.null(in_group)))
      inter_norm_factors <- switch(private$inter_norm_fact_opts$norm_scale,
        design = private$inter_norm_fact,
        batch = private$inter_norm_fact[[in_batch]],
        group = private$inter_norm_fact[[in_batch]][[in_group]]
      )
      sample_names <- private$design$extract_sample_names(in_batch, in_group)
      inter_norm_factors <- inter_norm_factors[sample_names]
      if (rescale && ! is.null(inter_norm_factors)) {
        inter_norm_factors <- inter_norm_factors /
          exp(mean(log(inter_norm_factors)))
      }
      inter_norm_factors
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
    #' @param ... Additional arguments passed to `$compute_norm`
    #' @return ggplot2 graph
    plot_dist_per_sample = function(
      geoms = c("histo", "boxplot"),
      tr_fn = (function(x) log2(x + 2) - 1),
      mean_fun = NULL,
      ...
    ) {
      sdesign <- private$design$get_simple_design()

      data <- purrr::map2_dfr(
        sdesign$batch,
        sdesign$group,
        ~ tibble::as_tibble(self$compute_norm(.x, .y, ...)) %>%
          dplyr::mutate(batch = .x, group = .y) %>%
          tidyr::pivot_longer(
            ! tidyselect::all_of(c("batch", "group")),
            names_to = "sample", values_to = "value"
          ) %>%
          dplyr::mutate(value = tr_fn(value))
      )

      plot_dist_per_sample_helper(
        data, private$design, geoms, mean_fun, deparse(body(tr_fn))
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
    #' @param rescale_inter_norm Rescale inter-sample normalization factors so
    #' that their geometric mean is 1. Default is TRUE.
    extract_pairwise_data_with_design = function(
      in_batch, in_group, rescale_inter_norm = TRUE
    ) {
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

      sf <- self$get_inter_norm_fact(
        in_batch, all_groups, rescale_inter_norm
      )

      # mandatory for rescaling to be OK
      assert_that(all(names(sf) %in% results$all_samples))

      results$inter_norm_fact <- sf[results$all_samples]

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
    #' @param tr_fn data transformation function
    #' @param plot_scale At which scale are the samples correlated.
    #' @param ggplot_mod ggplot modifier that will be added to the graph.
    #' @param color_palette name of the palette of color to use
    #' @param prcomp_args args to be used with prcomp.
    #' @param prcomp_autoplot_args list of argument that are used with autoplot.
    #' @param include_ctrl_at_group_scale whether to include the controls at
    #' group scale
    #' @param tags vector of tag ids to use for the plot
    #' @param tag_type name of the tag type to use
    #' @param pca_plot_dims integer id the the principal component to plot
    #' @param mshape Which variable to use for shape mapping
    #' @param mcolor Which variable to use for color mapping
    #' @param point_size point size for geom plot
    #' @return ggplot2 graph
    plot_prcomp = function(
      in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL, color_palette = "Paired",
      prcomp_args = list(), include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL, pca_plot_dims = c(1L, 2L),
      mshape = setNames("b_label", "Batch"),
      mcolor = setNames("g_label", "Group"), point_size = 5
    ) {
      private$plot_complex(
        "prcomp", in_batch, in_group, intra_norm, inter_norm, tr_fn,
        plot_scale, ggplot_mod, color_palette,
        prcomp_args,
        dist_method = NULL, hclust_method = NULL, dim_reduce = NULL,
        clust_bar_var = NULL, height_main = NULL, width_main = NULL,
        include_ctrl_at_group_scale,
        tags, tag_type, pca_plot_dims, mshape, mcolor, point_size
      )
    },

    #' @description
    #' Function to draw correlation plots
    #' @param tr_fn data transformation function
    #' @param plot_scale At which scale are the samples correlated.
    #' @param ggplot_mod ggplot modifier that will be added to the graph.
    #' @param include_ctrl_at_group_scale whether to include the controls at
    #' group scale
    #' @param tags vector of tag ids to use for the plot
    #' @param tag_type name of the tag type to use
    #' @return ggplot2 graph
    plot_corr = function(
      in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL,
      include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL
    ) {
      private$plot_complex(
        "corr", in_batch, in_group, intra_norm, inter_norm, tr_fn,
        plot_scale, ggplot_mod, color_palette = NULL,
        prcomp_args = NULL,
        dist_method = NULL, hclust_method = NULL, dim_reduce = NULL,
        clust_bar_var = NULL, height_main = NULL, width_main = NULL,
        include_ctrl_at_group_scale,
        tags, tag_type,
        pca_plot_dims = NULL, mshape = NULL, mcolor = NULL, point_size = NULL
      )
    },

    #' @description
    #' Function to draw hierarchical clustering plots
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
    #' @param include_ctrl_at_group_scale whether to include the controls at
    #' group scale
    #' @param tags vector of tag ids to use for the plot
    #' @param tag_type name of the tag type to use
    #' @return ggplot2 graph
    plot_hclust = function(
      in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL, color_palette = "BrBG",
      dist_method = "euclidean", hclust_method = "ward.D2", dim_reduce = NULL,
      clust_bar_var = c(), height_main = 10, width_main = 4,
      prcomp_args = list(),
      include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL
    ) {
      private$plot_complex(
        "hclust", in_batch, in_group, intra_norm, inter_norm, tr_fn,
        plot_scale, ggplot_mod, color_palette,
        prcomp_args = prcomp_args,
        dist_method, hclust_method, dim_reduce,
        clust_bar_var, height_main, width_main,
        include_ctrl_at_group_scale,
        tags, tag_type,
        pca_plot_dims = NULL, mshape = NULL, mcolor = NULL, point_size = NULL
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
    # main expression tag which is the row-nams-of data frame raw, len and
    # norm_fact
    main_etag = NULL,
    plot_complex = function(
      plot_type, in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL, color_palette = "Paired",
      prcomp_args = list(),
      dist_method = "euclidean", hclust_method = "ward.D2", dim_reduce = NULL,
      clust_bar_var = c(), height_main = 10, width_main = 4,
      include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL,
      pca_plot_dims, mshape = "batch", mcolor = "group", point_size = 5
    ) {

      # Prepare tag selection
      selected_tag_ids <- prepare_tag_selection(
        private$main_etag, private$annotation,
        private$selected_ids, tags, tag_type
      )

      sdesign <- private$design$get_simple_design(include_ctrl = TRUE)

      data <- switch(plot_scale,
        design = list(list(
          norm = self$compute_norm(
            intra_norm = inter_norm, inter_norm = inter_norm
          )[selected_tag_ids, ],
          design = private$design$get_pairwise_design(),
          title = "Design"
        )),
        batch = purrr::map(
          private$design$list_batches(),
          ~ list(
            norm = self$compute_norm(
              in_batch = .x, intra_norm = inter_norm, inter_norm = inter_norm
            )[selected_tag_ids, ],
            design = private$design$get_pairwise_design(in_batch = .x),
            title = private$design$get_b_labels()[.x]
          )
        ),
        group = purrr::map2(
          sdesign$batch,
          sdesign$group,
          ~ list(
            norm = self$compute_norm(
              in_batch = .x, in_group = .y,
              include_ctrl = include_ctrl_at_group_scale,
              intra_norm = inter_norm, inter_norm = inter_norm,
            )[selected_tag_ids, ],
            design = private$design$get_pairwise_design(in_batch = .x, in_group = .y),
            title = paste(private$design$get_b_labels()[.x], private$design$get_g_labels()[.y])
          )
        )
      )
      names(data) <- purrr::map_chr(data, ~ .x$title)

      graphs <- setNames(purrr::map(
        data,
        ~ switch(plot_type,
          prcomp = create_prcomp_plot(
            .x$norm, .x$title, .x$design, pca_plot_dims, mshape,
            mcolor, prcomp_args, ggplot_mod, color_palette, point_size
          ),
          corr = create_corr_plot(.x$norm, .x$title, ggplot_mod),
          hclust = create_hclust_plot(
            .x$norm, .x$title, .x$design, dist_method, hclust_method,
            dim_reduce, clust_bar_var, height_main, width_main, color_palette,
            prcomp_args, ggplot_mod
          )
        )
      ), names(data))

      if (plot_scale == "design") {
        return(graphs[[1]])
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
      private$raw <- as.data.frame(
        private$raw %>%
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
        g_norm <- as.data.frame(
          norm %>%
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
        private$len <- as.data.frame(
          private$len %>%
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


