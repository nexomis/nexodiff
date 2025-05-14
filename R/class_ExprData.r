#' @include utils.r
#' @include class_Annotation.r
#' @include class_PairwiseDesign.r

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
#' @param dist_method distand method, see dist function
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
  public <- list(

    #' @title Show raw count summary per tax_id, tax_name and rna type for
    #' selected expressed tags
    #' @description
    #' This function computes a raw count summary per tax_id, tax_name and rna
    #' type for selected expressed tags, and returns a list of count tables per
    #' feature.
    #' @param type A string specifying the type of summary to compute. Possible
    #' values are:
    #' - "etags" for the number of expression tags
    #' - "raw" for the sum of raw reads counts per samples
    #' - "norm" for the sum of normalized reads counts per samples.
    #' @return A list of count tables per feature
    #' @examples
    #' show_etags_summary("raw")
    show_etags_summary = function(type = "etags",
      in_batch = NULL) {
      results <- list()

      etag_id <- if (private$at_gene_level) "tgid" else "txid"
      etags <- private$selected_ids

      if (type == "etags") {
        for (metric_id in c("tax_id", "tax_name", "type")) {
          results[[metric_id]] <- table(
            private$annotation$generate_translate_dict(etag_id, metric_id)[etags]
          )
        }
        results$tax_id <- table(
          private$annotation$generate_translate_dict(etag_id, "tax_id")[etags])
        results$tax_name <- table(
          private$annotation$generate_translate_dict(etag_id, "tax_name")[etags])
        results$type <- table(
          private$annotation$generate_translate_dict(etag_id, "type")[etags])
      } else {
        if (type == "raw") {
          data <- as.data.frame(self$filter_and_get_raw(in_batch))
        } else if (type == "norm") {
          data <- as.data.frame(self$compute_norm(in_batch))
        } else {
          logging::logerror("Unknown type for summary")
          stop()
        }
        for (metric_id in c("tax_id", "tax_name", "type")) {
          data[, metric_id] <-
            private$annotation$generate_translate_dict(
              etag_id, metric_id)[row.names(data)]
          results[[metric_id]] <-
            aggregate(formula(paste(". ~", metric_id)), data, sum)
          row.names(results[[metric_id]]) <- results[[metric_id]][, metric_id]
          results[[metric_id]][, metric_id] <- NULL
          data[, metric_id] <- NULL
        }
      }
    results
    },

    #' @description
    #' Select expressed tags based on filtering on taxon or rna type.
    #' The results will be an intersect with the previous selection.
    #' You can reset the object if it's not desired.
    #' @param values A list of values used for filtering.
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
      values, filtered_var = "type", filter_type = "keep") {
      etag_id <- private$main_etag

      if (! filtered_var %in% c("tax_id", "tax_name", "type")) {
        logging::logerror("Wrong value for filtered_var argument")
        stop()
      }

      etags <- row.names(private$raw)
      filter_dict <- private$annotation$generate_translate_dict(
        etag_id, filtered_var)

      if (filter_type == "keep") {
        keep_ids <- names(filter_dict)[filter_dict %in% values]
      } else if (filter_type == "excl") {
        keep_ids <- names(filter_dict)[!filter_dict %in% values]
      } else {
        logging::logerror("Wrong value for filter_type argument")
        stop()
      }

      private$selected_ids <- intersect(keep_ids, private$selected_ids)
    },

    #' @description
    #' Compute and set inter normalization factors with extensible
    #' parametrization. Please be aware that this function is designed to handle
    #'  various normalization methods, including median and TMM. However, at the
    #'  time of implementation, only the median method has been thoroughly
    #' tested. Other methods such as TMM might not produce accurate results and
    #' should be used with caution.
    #' @param method Normalization method to use:
    #'   * "none" apply no normalization
    #'   * "median" apply median ratio normalization
    #'   * "tmm" apply trimmed mean of median normalization
    #' TODO: add litterature reference for each
    #' @param norm_scale Scale to apply the normalization:
    #' - "design" for all samples in one go
    #' - "batch" for applying norm_factors per batch
    #' - "group" for applying norm factors per batch per group.
    #' @param norm_by Where to apply the normalization factor:
    #' - "samples" for the sample level
    #' - "group" for the group level (not implemented).
    #' The normalized can be computed for each sample or can be computed based
    #' on the mean expression within a group.
    #' @param norm_ref Sample(s) to use as reference for normalization:
    #' - "all" for all samples
    #' - "ctrl" for control samples (not tested and not compatible with scale
    #' design)
    #' - character vector with sample names (not tested)
    #' - integer vector with sample rank in expression matrix (not tested).
    #' @param norm_ref_mean Method for the mean of gene expression between
    #' samples of the reference. Possible options are given for `m_trim_mean`.
    #' @param m_trim_prop Proportion of expression ratio (M-value) to trim at
    #' both tails (only for TMM).
    #' @param m_trim_mean Method for the mean of gene expression ratios in the
    #' TMM method to compute the normalization factor. Possible options are:
    #' - "median"
    #' - "geometric"
    #' - "nz.geometric" geometric without zero
    #' - "mod.geometric" modified geometric with epsilon = 1e-05
    #' (see https://arxiv.org/abs/1806.06403)
    #' - "arithmetic"
    #' @param a_trim_mean Method for the mean of gene expression between all
    #' samples. Possible options are given for `m_trim_mean`.
    #' @param a_trim_norm Whether to use (intra) normalized value to trim based
    #' on mean gene expression (only FALSE tested).
    #' @param ncpus Number of cpus to use for computation (default = 1)
    #' @param a_trim_value Threshold value to trim genes from the calculation of
    #' the normalization factors based on mean gene expression (A-value).
    #' @param replace_zero_by Replace all zero raw count by the given value.
    compute_and_set_inter_norm_fact = function(method = "median",
      norm_scale = "group", norm_by = "sample", norm_ref = "all",
      norm_ref_mean = "mod.geometric", m_trim_prop = 0.3,
      m_trim_mean = "mod.geometric", a_trim_value = 0.5, a_trim_norm = FALSE,
      a_trim_mean = "mod.geometric", replace_zero_by = 0, ncpus = 1) {

      if (method == "tmm" & m_trim_mean == "median") {
        logging::logwarning(paste("tmm method with median as mean function...",
          "Is that not simply the median method ?")
        )
      }

      if (method == "tmm" & m_trim_prop == 0) {
        logging::logerror("`m_trim_prop` cannot be 0 with `method` tmm")
        stop()
      }

      if (m_trim_prop <= 0 | m_trim_prop >= 0.5) {
        logging::logerror("`m_trim_prop must` be between 0 and 0.5 excluded")
        stop()
      }

      norm_ref_mean_fun <- set_mean_function(norm_ref_mean)
      a_trim_mean_fun <- set_mean_function(a_trim_mean)
      m_trim_mean_fun <- set_mean_function(m_trim_mean)

      if (method == "median") {
        calc_norm_fact <- function(norm, sample_ref) {
          # filter out when sample_ref = 0
          keep_row <- sample_ref != 0
          norm <- norm[keep_row, ] / sample_ref[keep_row]
          norm_fact <- 1 / matrixStats::colMedians(norm)
          names(norm_fact) <- colnames(norm)
          return(norm_fact)
        }
      } else if (method == "tmm") {
        logging::logerror("TMM method is not implemented yet")
        stop()
      }

      if (norm_by == "group") {
        logging::logerror("`norm_by` group is not yet implemented")
        stop()
      }

      get_refrence <- function(norm, ctrl_samples) {
        if (all(norm_ref == "all")) {
          norm_ref_samples <- colnames(norm)
        } else if (all(ref == "ctrl")) {
          norm_ref_samples <- ctrl_samples
        } else if (all(ref %in% colnames(norm)) |
          all(ref %in% seq_len(ncol(norm)))
        ) {
          norm_ref_samples <- ref
        } else {
          logging::logerror("`ref` argument not recognized")
          stop()
        }
        norm_ref <- unlist(lapply(seq(1, nrow(norm)), function(i) {
          norm_ref_mean_fun(norm[i, norm_ref_samples])
        }))
        return(norm_ref)
      }



      get_trimmed_mat <- function(raw, norm_fact) {

        trimmed_mat <- raw * norm_fact

        if (a_trim_norm) {
          base_mat <- trimmed_mat
        } else {
          base_mat <- raw
        }
        base_mat_mean <-  unlist(lapply(seq(1, nrow(base_mat)), function(i) {
          a_trim_mean_fun(base_mat[i, ])
        }))
        trimmed_mat[trimmed_mat == 0] <- replace_zero_by
        trimmed_mat[base_mat_mean > a_trim_value, ]
      }

      if (ncpus > parallelly::availableCores() - 1) {
        logging::logwarn("The number of core asked was more than available")
      }

      if (ncpus == 1) {
        future::plan(strategy = future::sequential)
      } else if (parallelly::supportsMulticore()) {
        future::plan(
          strategy = future::multicore,
          workers = min(ncpus, max(1, parallelly::availableCores() - 1))
        )
      } else {
        future::plan(
          strategy = future::multisession,
          workers = min(ncpus, max(1, parallelly::availableCores() - 1))
        )
      }

      if (norm_scale == "group") {
        map_base <- dplyr::distinct(private$design$get_pairwise_design(
          )[, c("batch", "group")])
        tmp <- furrr::future_map2(
          map_base$batch,
          map_base$group,
          function(batch, group) {
            res <- self$extract_pairwise_data_with_design(batch, group)
            # here we only use intra_norm !
            norm <- get_trimmed_mat(res$raw, res$intra_norm_fact)
            sample_ref <- get_refrence(norm, res$ctrl_samples)
            calc_norm_fact(norm, sample_ref)
          }
        )
        names(tmp) <- paste(map_base$batch, map_base$group, sep = "")
        private$inter_norm_fact <- list()

        for (batch in private$design$list_batches()) {
          private$inter_norm_fact[[batch]] <- list()
          for (group in private$design$list_groups_per_batches(
            include_ctrl = TRUE)[[batch]]) {
            private$inter_norm_fact[[batch]][[group]] <-
              tmp[[paste(batch, group, sep = "")]]
          }
        }

        names(private$inter_norm_fact) <- private$design$list_batches()
        for (batch in private$design$list_batches()) {
          names(private$inter_norm_fact[[batch]]) <-
            private$design$list_groups_per_batches(include_ctrl = TRUE)[[batch]]
        }
      } else if (norm_scale == "batch") {
        private$inter_norm_fact <- furrr::future_map(
          private$design$list_batches(),
          function(batch) {
            norm <- get_trimmed_mat(self$filter_and_get_raw(batch),
              self$compute_norm_fact(batch))
            sample_ref <- get_refrence(norm,
              private$design$find_control_group_per_batches()[batch])
            calc_norm_fact(norm, sample_ref)
          }
        )
        names(private$inter_norm_fact) <- private$design$list_batches()
      } else if (norm_scale == "design") {
        norm <- get_trimmed_mat(self$get_raw(),
          self$compute_norm_fact())

        if (norm_ref == "ctrl") {
          logging::logerror(
            "`norm_ref` ctrl is not compatible with `norm_scale` design"
          )
          stop()
        }
        sample_ref <- get_refrence(norm, NULL)
        private$inter_norm_fact <-  calc_norm_fact(norm, sample_ref)
      }

      future::plan(future::sequential)

      private$inter_norm_fact_opts <- list(
        method = method,
        norm_scale = norm_scale,
        norm_by = norm_by,
        norm_ref = norm_ref,
        norm_ref_mean = norm_ref_mean,
        m_trim_prop = m_trim_prop,
        m_trim_mean = m_trim_mean,
        a_trim_value = a_trim_value,
        a_trim_norm = a_trim_norm,
        a_trim_mean = a_trim_mean,
        replace_zero_by = replace_zero_by
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
      for (batch in private$design$list_batches()) {
        if (method == "tpm") {
          private$intra_norm_fact[private$selected_ids, ] <-
            sweep(1e+06 / self$get_len(),
              2,
              colSums(self$get_raw() /
                self$get_len()),
              FUN = "/"
            )
        } else if (method == "fpk") {
          private$intra_norm_fact[private$selected_ids, ] <-
            1e+03 / self$get_len()
        } else if (method == "fpm") {
          private$intra_norm_fact[private$selected_ids, ] <-
            sweep(
              matrix(1e+06,
                nrow = length(private$selected_ids),
                ncol = ncol(private$raw),
                dimnames = list(
                  private$selected_ids,
                  colnames(private$raw)
                )
              ),
              2,
              colSums(self$get_raw()),
              FUN = "/"
            )
        } else if (method == "fpkm") {
          private$intra_norm_fact[private$selected_ids, ] <-
            sweep(1e+09 / self$get_len(),
              2,
              colSums(self$get_raw()),
              FUN = "/"
            )
        } else if (method == "none") {
          private$intra_norm_fact[private$selected_ids, ][] <- 1
        }else {
          logging::logdebug("Unrecognized normalization method")
          stop()
        }

        # There will be NaN value where colSums(private$raw) or
        # colSums(private$raw/private$len) is 0
        # We are going to replace those NaN value by 0
        private$intra_norm_fact[is.nan(private$intra_norm_fact)] <- 1
      }
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
      inter_norm = FALSE, intra_norm = TRUE, include_ctrl = FALSE) {

      row_ids <- private$selected_ids

      all_group <- in_group

      if ((length(in_batch) != 1 | length(in_group) != 1) & include_ctrl) {
        logging::logerror(paste(
          "include_ctrl set to TRUE is possible only when specifying",
          "one and only one in_group and in_batch args"
        ))
        stop()
      }

      if (include_ctrl) {
        if ((! is.null(in_batch)) & (! is.null(in_group))) {
          ctrl_group <- private$design$find_control_group_per_batches(
            )[in_batch]
          all_group <- c(in_group, ctrl_group)
        }
      }

      col_ids <- private$design$extract_sample_names(in_batch, all_group)

      if (intra_norm) {
        norm_fact <- private$intra_norm_fact[row_ids, col_ids]
      } else {
        norm_fact <- NULL
      }

      norm_fact_inter <- NULL

      if (inter_norm) {
        if (! is.null(private$inter_norm_fact)) {
          if (private$inter_norm_fact_opts$norm_scale == "design") {
            norm_fact_inter <- private$inter_norm_fact[col_ids]
          } else {
            if (is.null(in_batch)) {
              logging::logerror(
                "`in_batch` cannot be null when norm scale is not design"
              )
              stop()
            }
            if (private$inter_norm_fact_opts$norm_scale == "batch") {
              norm_fact_inter <- private$inter_norm_fact[[in_batch]][col_ids]
            } else {
              if (is.null(in_group)) {
                logging::logerror(
                  "`in_group` cannot be null when norm scale is group"
                )
                stop()
              }
              norm_fact_inter <-
                private$inter_norm_fact[[in_batch]][[in_group]][col_ids]
            }
          }
        } else {
          norm_fact_inter <- rep(1, length(col_ids))
          names(norm_fact_inter) <- col_ids
        }
      }
      if (! is.null(norm_fact_inter)) {
        if (intra_norm) {
          norm_fact <- norm_fact %*% diag(norm_fact_inter[col_ids])
          colnames(norm_fact) <- col_ids
        } else {
          norm_fact <- norm_fact_inter[col_ids]
        }
      }
      if (is.null(norm_fact)) {
        logging::logerror(
          "`intra_norm` and `inter_norm` cannot be FALSE together"
        )
        stop()
      }
      norm_fact
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
    #' 
    #' - "compute_and_set_intra_norm_fact"
    #' - "compute_and_set_inter_norm_fact"
    #' Note that intra normalization is always applied, however inter norm is
    #' optional. Depending on the `norm_scale` used for inter normalization,
    #' the `in_batch` argument might be mandatory.
    #' @return numeric matrix with norm counts for expressed tags
    compute_norm = function(in_batch = NULL, in_group = NULL,
      inter_norm = FALSE, include_ctrl = FALSE) {
      if (include_ctrl) {
        ctrl_group <- private$design$find_control_group_per_batches()[in_batch]
      } else {
        ctrl_group <- c()
      }
      norm_fact <-
        self$compute_norm_fact(in_batch, in_group, inter_norm = inter_norm,
          include_ctrl = include_ctrl)
      norm <- self$filter_and_get_raw(in_batch, c(in_group, ctrl_group))
      return(norm * norm_fact[, colnames(norm)])
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
    #' Get sum of count or normalized value per rna type
    #' @return data frame with counts with colums
    #'  - type: RNA type
    #'  - sample: sample name
    #'  - total: Sum of count or normalized value
    sum_per_type_per_sample = function(intra_norm = FALSE, log2_expr = FALSE) {
      to_type <- private$annotation$generate_translate_dict(
        private$main_etag, "type")
      format_data <- function(data) {
        data$type <- to_type[row.names(data)]
          temp <- tidyr::pivot_longer(aggregate(. ~ type, data, sum),
            cols = -type, names_to = "sample", values_to = "total")
         return(subset(temp, total > 10))
      }

      if (intra_norm) {
        data <- purrr::map_dfr(
          private$design$list_batches(),
          function(batch) {
            format_data(as.data.frame(self$compute_norm(in_batch = batch)))
          }
        )
      } else {
        data <- format_data(as.data.frame(self$get_raw()))
      }
      if (log2_expr) {
        data$total <- log2(data$total + 2)
      }
      data
    },

    #' @description
    #' plot sum of count or intra-normalized values per rna type (Inter norm is
    #' not available for this)
    #' @param same_scale whether to represent the scale with the same
    #' @param exclude_type character vector with the rne type to exclude from
    #' the graph
    #' @param horizonal_bar whether the bar in the plot are horizontal rather
    #' than vertical
    #' @return ggplot2 graph
    plot_sum_per_type_per_sample = function(intra_norm = FALSE,
      exclude_type = c(), same_scale = TRUE, horizonal_bar = TRUE,
      log2_expr = FALSE) {
      data_design <-
        private$design$get_pairwise_design()[, c("batch", "group", "sample")]
      data <- dplyr::left_join(
        data_design,
        self$sum_per_type_per_sample(intra_norm = intra_norm, log2_expr),
        by = "sample") %>%
        dplyr::filter(! .data$type %in% exclude_type)
      batch2label <- private$design$get_b_labels()
      group2label <- private$design$get_g_labels()
      y_label <- "Sum of expression"
      if (log2_expr) {
        y_label <- paste(y_label, "(log2-transformed)")
      }
      if (! same_scale) {
        f_scales <- "free"
      } else if (horizonal_bar) {
        f_scales <- "free_y"
      } else {
        f_scales <- "free_x"
      }
      data$batch <- factor(data$batch, levels = names(batch2label))
      data$group <- factor(data$group, levels = names(group2label))
      g <- ggplot2::ggplot(data,
        ggplot2::aes(.data$sample, .data$total, fill = .data$type)
      ) +
        ggplot2::geom_bar(stat = "identity") +
        ggh4x::facet_nested_wrap(
          formula("~ batch + group"),
          nrow = length(unique(data$batch)),
          scales = f_scales
          ,
          labeller = ggplot2::labeller(
            batch = ggplot2::as_labeller(batch2label),
            group = ggplot2::as_labeller(group2label)
          )
        ) +
        ggplot2::scale_y_continuous(
          labels = scales::label_scientific(digits = 2)) +
        ggplot2::scale_fill_brewer(palette = "BrBG") +
        THEME_NEXOMIS +
        ggplot2::labs(
          fill = "RNA type",
          x = "Samples",
          y = y_label
        )
      if (horizonal_bar) {
        g <- g + ggplot2::coord_flip()
      }
      g
    },

    #' @description
    #' Plot the distribution of count or normalized values per RNA type for each
    #' sample.
    #' @param geoms Character vector, types of plot to draw, can be a
    #' combination of:
    #' - "boxplot": draw a box plot (fast)
    #' - "violin": draw a violin plot (slow)
    #' - "histo": draw an hiostogram (slow)
    #' @param mean_fun (Optional) Function to plot the mean per sample.
    #' Can be one of the following:
    #' * "median"
    #' * "geometric"
    #' * "nz.geometric" geometric without zero
    #' * "mod.geometric" modified geometric with epsilon = 1e-05
    #' (see https://arxiv.org/abs/1806.06403)
    #' * "arithmetic"
    #' @return ggplot2 graph
    plot_dist_per_sample = function(intra_norm = TRUE, inter_norm = TRUE,
      log2_expr = TRUE, geoms = c("boxplot"),
      mean_fun = NULL) {

      if (length(setdiff(geoms, c("boxplot", "violin", "histo"))) > 0) {
        logging::logerror("unrecognized value for geoms")
        stop()
      }

      data_design <-
        private$design$get_pairwise_design()[, c("batch", "group", "sample")]
      data_design_batch_group <- dplyr::distinct(
        private$design$get_pairwise_design()[, c("batch", "group")]
      )

      y_label <- "Expression"

      if (log2_expr) {
        y_label <- paste(y_label, "(log2-transformed)")
      }
      batch2label <- private$design$get_b_labels()
      group2label <- private$design$get_g_labels()
      get_data <- function(in_batch, in_group) {
        results <- self$extract_pairwise_data_with_design(in_batch, in_group,
          include_ctrl = FALSE)
        keep_samples <- results$test_samples
        data <- results$raw[, keep_samples]
        if (intra_norm) {
           data <- data *
            results$intra_norm_fact[, keep_samples]
        }
        if (inter_norm) {
          data <- data %*%
            diag(results$inter_norm_fact[keep_samples])
        }
        if (log2_expr) {
          data <- log2(data + 2)
        }
        data <- as.data.frame(data)
        data$batch <- in_batch
        data$group <- in_group
        data %>%
          tidyr::pivot_longer(
            ! tidyselect::all_of(c("batch", "group")),
            names_to = "sample",
            values_to = "value"
          )
      }

      data <- purrr::map2_dfr(
        data_design_batch_group$batch,
        data_design_batch_group$group,
        get_data
      )

      data$batch <- factor(data$batch, levels = names(batch2label))
      data$group <- factor(data$group, levels = names(group2label))

      g <- ggplot2::ggplot(data,
        ggplot2::aes(.data$sample, .data$value)
      )
      if ("violin" %in% geoms) {
        g <- g +
          ggplot2::geom_violin(trim = FALSE)
      }
      if ("boxplot" %in% geoms) {
        g <- g +
          ggplot2::geom_boxplot(width = 0.1)
      }
      if ("histo" %in% geoms) {
        g <- g +
          ggmulti::geom_histogram_(bins = 30, alpha = 0.5)
      }
      if (! is.null(mean_fun)) {
        g <- g + ggplot2::stat_summary(
          fun = set_mean_function(mean_fun),
          geom = "point", shape = 21, size = 3,
          color = "black", fill = "red"
        )
      }

      g <- g +
        ggh4x::facet_nested_wrap(
          formula("~ batch + group"),
          nrow = length(unique(data$batch)),
          scales = "free_x",
          labeller = ggplot2::labeller(
            batch = ggplot2::as_labeller(batch2label),
            group = ggplot2::as_labeller(group2label)
          )
        ) +
        ggplot2::scale_y_continuous(
          labels = scales::label_scientific(digits = 2)) +
        THEME_NEXOMIS +
        ggplot2::labs(
          x = "Samples",
          y = y_label
        )
      return(g)
    },


    #' @description
    #' Function to draw complex plots
    #' @param color_palettes name of the palette of color to use for plot in the
    #' order of variable added.
    #' @return ggplot2 graph
    plot_complex = function(plot_type, in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL, color_palettes = c("Set1",
      "Set2", "Set3", "Pastel1", "Pastel2", "Paired", "Dark2", "Accent"),
      prcomp_args = list(), prcomp_autoplot_args = list(),
      dist_method = "euclidean", hclust_method = "ward.D2", dim_reduce = NULL,
      clust_bar_var = c(), height_main = 10, width_main = 4,
      include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL, df_design_filter = NULL) {

      if (is.null(tag_type)) {
        dict_ids <- private$selected_ids
        names(dict_ids) <- dict_ids
      } else {
        dict_ids <- private$annotation$generate_translate_dict(
          private$main_etag,
          tag_type
        )
      }

      if (is.null(tags)) {
        selected_ids <- names(dict_ids)
      } else {
        selected_ids <- names(dict_ids)[dict_ids %in% tags]
      }

      selected_ids <- intersect(selected_ids, private$selected_ids)

      if (length(selected_ids) < 3) {
        logging::logerror("tags must have at least 3 valid choice")
        stop()
      }

      make_plot_complex_args <- list(
        prcomp_args = prcomp_args,
        prcomp_autoplot_args = prcomp_autoplot_args,
        ggplot_mod = ggplot_mod,
        plot_type = plot_type,
        dist_method = dist_method,
        hclust_method = hclust_method,
        dim_reduce = dim_reduce,
        clust_bar_var = clust_bar_var,
        height_main = height_main,
        width_main = width_main,
        color_palettes = color_palettes
      )

      tr_fn_df <- function(df) {
        rn <- row.names(df)
        cn <- colnames(df)
        df <- tr_fn(as.matrix(df))
        row.names(df) <- rn
        colnames(df) <- cn
        return(df)
      }

      data_design <-
        private$design$get_pairwise_design()[, c("batch", "group")] %>%
        dplyr::distinct()
      if (! is.null(in_batch)) {
        data_design <- dplyr::filter(data_design,
          .data$batch %in% in_batch)
      }
      if (! is.null(in_group)) {
        data_design <- dplyr::filter(data_design,
          .data$group %in% in_group)
      }

      if (!is.null(df_design_filter)) {
        if (plot_scale == "batch" && (! "batch" %in% names(df_design_filter))) {
          logging::logerror(
            "at batch scale, df_design_filters must precise batch for filters")
          stop()
        } else if (plot_scale == "group" &&
          (! all(c("batch", "group") %in% names(df_design_filter)))
        ) {
          logging::logerror(paste(
            "at group scale,",
            "df_design_filters must precise batch and group for filters"))
          stop()
        }
        inner_join_cols <- c()
        for (c_id in c("batch", "group")) {
          if (c_id %in% names(df_design_filter)) {
            inner_join_cols <- c(inner_join_cols, c_id)
          }
        }

        if (length(inner_join_cols) != 0) {
          data_design <- dplyr::inner_join(
            data_design,
            df_design_filter,
            by = inner_join_cols,
            multiple = "all",
            suffix = c("", ".filter")
          )[, c("batch", "group")] %>%
            dplyr::distinct()
        }
      } else {
        df_design_filter <- private$design$get_pairwise_design()
      }

      y_label <- "Expression"

      if (plot_scale == "group") {
        graphs <- purrr::map2(
          data_design$batch,
          data_design$group,
          function(x, y) {
            results <- self$extract_pairwise_data_with_design(
              x, y, include_ctrl = include_ctrl_at_group_scale)
            selected_samples <- intersect(
              unique(dplyr::filter(df_design_filter,
                .data$batch == .env$x & .data$group == .env$y)$sample),
              colnames(results$raw)
            )
            in_data <- results$raw[selected_ids, selected_samples]
            rn <- row.names(in_data)
            cn <- colnames(in_data)
            if (intra_norm) {
              in_data <- in_data *
                results$intra_norm[selected_ids, selected_samples]
            }
            if (inter_norm) {
              in_data <- in_data %*%
                diag(results$inter_norm_fact[selected_samples])
            }
            in_data <- tr_fn_df(in_data)
            row.names(in_data) <- rn
            colnames(in_data) <- cn
            in_title <- paste(
              private$design$get_b_labels()[x],
              private$design$get_g_labels()[y],
              sep = ": "
            )
            make_plot_complex_args$in_data <- in_data
            make_plot_complex_args$in_title <- in_title
            do.call(private$make_plot_complex, make_plot_complex_args)
          }
        )
      } else if (plot_scale == "batch") {
        if (inter_norm & private$inter_norm_fact_opts$norm_scale == "group") {
          logging::logerror(paste(
            "plot_scale=batch is not compatible with group-scale inter",
            "normalization")
          )
          stop()
        }
        graphs <- purrr::map(
          unique(data_design$batch),
          function(x) {
            if (intra_norm) {
              data <- self$compute_norm(x, in_group, inter_norm, FALSE)
            } else {
              data <- self$filter_and_get_raw(x, in_group)
            }
            selected_samples <- intersect(
              unique(dplyr::filter(df_design_filter,
                .data$batch == .env$x)$sample),
              colnames(data)
            )
            data <- data[selected_ids, selected_samples]
            data <- tr_fn_df(data)
            in_title <- private$design$get_b_labels()[x]
            make_plot_complex_args$in_data <- data
            make_plot_complex_args$in_title <- in_title
            n_args <- names(make_plot_complex_args$prcomp_autoplot_args)
            if (! "colour" %in% n_args) {
              make_plot_complex_args$prcomp_autoplot_args$colour <- "g_label"
            }
            do.call(private$make_plot_complex, make_plot_complex_args)
          }
        )
      } else if (plot_scale == "design") {
        if (inter_norm & private$inter_norm_fact_opts$norm_scale != "design") {
          logging::logerror(paste(
            "plot_scale=design is not compatible with group-scale or",
            "batch-scale in normalization")
          )
          stop()
        }
        if (intra_norm) {
          data <- self$compute_norm(in_batch, in_group, inter_norm, FALSE)
        } else {
          data <- self$filter_and_get_raw(in_batch, in_group)
        }
        selected_samples <- intersect(
          unique(df_design_filter$sample),
          colnames(data)
        )
        data <- data[selected_ids, selected_samples]
        data <- tr_fn_df(data)
        in_title <- "Design"
        make_plot_complex_args$in_data <- data
        make_plot_complex_args$in_title <- in_title
        n_args <- names(make_plot_complex_args$prcomp_autoplot_args)
        if (! "colour" %in% n_args) {
          make_plot_complex_args$prcomp_autoplot_args$colour <- "g_label"
        }
        if (! "shape" %in% n_args) {
          make_plot_complex_args$prcomp_autoplot_args$shape <- "b_label"
        }
        graphs <-
          list(do.call(private$make_plot_complex, make_plot_complex_args))
      } else {
        logging::logerror("plot_scale argument not valid")
      }

      p <- cowplot::plot_grid(
          plotlist = graphs,
          nrow = as.integer(sqrt(length(graphs)))
        )

      return(p)

    },

    #' @description
    #' plot principal components
    plot_prcomp = function(in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL,
      prcomp_args = list(), prcomp_autoplot_args = list(),
      include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL, df_design_filter = NULL) {

      self$plot_complex(
        "prcomp",
        in_batch = in_batch,
        in_group = in_group,
        intra_norm = intra_norm,
        inter_norm = inter_norm,
        tr_fn = tr_fn,
        plot_scale = plot_scale,
        ggplot_mod = ggplot_mod,
        prcomp_args = prcomp_args,
        prcomp_autoplot_args = prcomp_autoplot_args,
        include_ctrl_at_group_scale = include_ctrl_at_group_scale,
        tags = tags,
        tag_type = tag_type,
        df_design_filter = df_design_filter
      )
    },

    #' @description
    #' plot correlations
    plot_corr = function(in_batch = NULL, in_group = NULL,
      intra_norm = TRUE, inter_norm = TRUE, tr_fn = (function(x) log2(x + 2)),
      plot_scale = "group", ggplot_mod = NULL,
      include_ctrl_at_group_scale = FALSE,
      tags = NULL, tag_type = NULL, df_design_filter = NULL) {

      self$plot_complex(
        "corr",
        in_batch = in_batch,
        in_group = in_group,
        intra_norm = intra_norm,
        inter_norm = inter_norm,
        tr_fn = tr_fn,
        plot_scale = plot_scale,
        ggplot_mod = ggplot_mod,
        include_ctrl_at_group_scale = include_ctrl_at_group_scale,
        tags = tags,
        tag_type = tag_type,
        df_design_filter = df_design_filter
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
    extract_pairwise_data_with_design = function(in_batch, in_group,
      include_ctrl = TRUE) {
      group_per_batches <-
        private$design$list_groups_per_batches(include_ctrl = TRUE)

      if ((length(in_batch) != 1) |
        ! (in_batch %in% names(group_per_batches))
      ) {
        logging::logerror(
          "in_batch is not recognized in design")
        stop()
      }
      if ((length(in_group) != 1) |
        ! (in_group %in% group_per_batches[in_batch][[1]])
      ) {
        logging::logerror(
          "in_group is not recognized in design")
        stop()
      }

      results <- list()

      # test if control group is the same... Sometimes we use this to get the
      # for a control group also.
      results$test_samples <-
        private$design$extract_sample_names(in_batch, in_group)

      # find batch control
      if (include_ctrl) {
        results$ctrl_group <-
          private$design$find_control_group_per_batches()[in_batch]
        results$ctrl_samples <-
          private$design$extract_sample_names(in_batch, results$ctrl_group)
        all_groups <- c(in_group, results$ctrl_group)
        results$all_samples <-
        unique(c(results$test_samples, results$ctrl_samples))
      } else {
        all_groups <- in_group
        results$all_samples <- results$test_samples
      }

      # norm is not computed because if the table is big, it's too long

      results$raw <- self$filter_and_get_raw(in_batch, all_groups)[,
        results$all_samples]
      results$len <- self$filter_and_get_len(in_batch, all_groups)[,
        results$all_samples]

      results$intra_norm_fact <-
        self$compute_norm_fact(in_batch, in_group, inter_norm = FALSE,
        intra_norm = TRUE, include_ctrl = include_ctrl)[, results$all_samples]

      results$inter_norm_fact <-
        self$compute_norm_fact(in_batch, in_group, inter_norm = TRUE,
        intra_norm = FALSE, include_ctrl = include_ctrl)[results$all_samples]

      results$design_table <- private$design$get_pairwise_design(
        in_batch, all_groups)

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
    make_plot_complex = function(in_data, in_title, plot_type, ggplot_mod,
      prcomp_args, prcomp_autoplot_args, dist_method, hclust_method, dim_reduce,
      clust_bar_var, height_main, width_main, color_palettes) {

      var_base <- names(private$design$get_pairwise_design())
      data_grouped <- private$design$get_pairwise_design() %>%
        dplyr::filter(.data$sample %in% colnames(in_data)) %>%
        dplyr::distinct() %>%
        dplyr::arrange(.data$sample)

      sample_order <- unique(data_grouped[["sample"]])

      data_grouped <- data_grouped %>%
        dplyr::group_by(.data$sample)

      data_design <- as.data.frame(purrr::map_dfc(
        setdiff(var_base, "sample"),
        function(x) {
          df <- dplyr::summarise(
            data_grouped, var = paste0(.data[[x]], collapse = ":")
          ) %>%
            dplyr::arrange(.data$sample)
          df$sample <- NULL
          names(df) <- x
          return(df)
        }
      ))

      data_design$sample <- sample_order

      if (plot_type == "corr") {
        lowerfun <- function(data, mapping) {
          ggplot2::ggplot(data = data, mapping = mapping) +
            ggplot2::geom_point(size = 0.5) +
            ggplot2::geom_abline(intercept = 0, slope = 1, color = "red")
        }

        g <- GGally::ggpairs(
          in_data,
          lower = list(continuous = GGally::wrap(lowerfun)),
          title = in_title
        ) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank()
        ) + THEME_NEXOMIS
        if (! is.null(ggplot_mod)) {
          g <- g + ggplot_mod
        }
        return(GGally::ggmatrix_gtable(g))
      } else if (plot_type == "prcomp") {
        prcomp_args$x <- t(in_data)
        pca_res <- do.call(prcomp, prcomp_args)
        prcomp_autoplot_args$object <- pca_res
        prcomp_autoplot_args$data <- data_design
        g <- do.call(ggplot2::autoplot, prcomp_autoplot_args) +
          ggplot2::ggtitle(in_title) +
          THEME_NEXOMIS +
          ggplot2::theme(legend.position = "bottom")
        if (! is.null(ggplot_mod)) {
          g <- g + ggplot_mod
        }
        return(g)
      } else if (plot_type == "hclust") {
        if (is.null(dim_reduce)) {
          d <- dist(t(in_data), method = dist_method)
        } else {
          prcomp_args$x <- t(in_data)
          prcomp_args$rank. <- as.integer(dim_reduce)
          pca_res <- do.call(prcomp, prcomp_args)
          str(pca_res)
          d <- dist(pca_res$x, method = dist_method)
        }
        hc <- hclust(d, hclust_method)

        ordered_labels <- hc$labels[hc$order]

        data_design$sample <- factor(
          data_design$sample,
          levels = ordered_labels,
          ordered = TRUE
        )

        p1_dendro <- ggdendro::dendro_data(hc)
        list_legs <- list()
        list_plots <- list(
          ggdendro::ggdendrogram(hc) +
            ggplot2::coord_cartesian(
              xlim = c(-1, nrow(data_design) + 1),
              ylim = c(-1, max(p1_dendro$segments$y)),
              expand = FALSE) +
            ggplot2::ggtitle(in_title) +
            THEME_NEXOMIS +
            ggplot_mod
        )
        hs <- c(height_main)
        ws <- integer()

        data_design <- as.data.frame(
          dplyr::arrange(data_design, .data$sample)
        )

        if (!identical(
          as.character(data_design$sample),
          as.character(ordered_labels)
        )) {
          print(data_design)
          str(data_design)
          print(hc$labels)
          str(hc$labels)
          print(as.character(data_design$sample))
          print(as.character(hc$labels))
          stop(
            "Error: The labels in data_design and hc are not in the same order."
          )
        }

        for (var in clust_bar_var) {
          if (! var %in% names(data_design)) {
            logging::logerror("var inclust_bar_var arg is not recognized")
            stop()
          }
          hs <- c(hs, 1)
          data_design$voi <- data_design[[var]]
          nth <- length(list_legs) + 1
          var_plot <-
            ggplot2::ggplot(data_design,
              ggplot2::aes(.data$sample, y = 1, fill = .data$voi)
            ) +
              ggplot2::ylab("") +
              ggplot2::geom_tile(color = "black") + ggplot2::theme_minimal() +
              ggplot2::coord_cartesian(
                xlim = c(-1, nrow(data_design) + 1),
                expand = FALSE) +
              ggplot2::scale_fill_brewer(palette = color_palettes[nth])
          list_legs[[nth]] <-
            ggpubr::as_ggplot(ggpubr::get_legend(
              var_plot +
                ggplot2::theme(
                  plot.margin = ggplot2::margin(0, 0, 0, 0, "cm")
                ) +
                ggplot2::guides(
                  fill = ggplot2::guide_legend(title = var, ncol = 1))
            ))
          ws <- c(ws, 1 + length(unique(data_design[[var]])))
          list_plots[[length(list_plots) + 1]] <- var_plot +
              ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_text(
                  angle = 0, vjust = 0.5, hjust = 1),
                axis.ticks = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                legend.position = "none",
                line = ggplot2::element_blank())
        }
        if (length(list_legs) == 0) {
          return(egg::ggarrange(
            plots = list_plots,
            ncol = 1,
            heights = hs, draw = FALSE
          ))
        } else {
          return(egg::ggarrange(
            ggpubr::as_ggplot(egg::ggarrange(
              plots = list_plots,
              ncol = 1,
              heights = hs, draw = FALSE
            )),
            ggpubr::as_ggplot(egg::ggarrange(
              plots = list_legs,
              ncol = 1,
              heights = ws,
              draw = FALSE
            )),
            ncol = 2,
            widths = c(width_main, 1),
            draw = FALSE
          ))
        }
      }
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
        expr_data_transcript$get_main_etag(), private$main_etag)

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
    initialize = function(design, annotation, with_fixed_length = FALSE, log_level = "WARN",
      format = "kallisto", suffix_pattern = "\\.\\d+$") {
      logging::basicConfig(log_level)

      private$design <- design
      private$annotation <- annotation
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

        if (!requireNamespace("rhdf5", quietly = TRUE)) {
          logging::logerror(paste("reading kallisto results from hdf5 files",
          "requires Bioconductor package `rhdf5`"))
          stop()
        }

        private$raw <- NULL
        private$selected_ids <- NULL

        for (sample_name in names(files)) {

          fpath <- files[[sample_name]]
          ids <- as.character(rhdf5::h5read(fpath, "aux/ids"))

          if (is.null(private$selected_ids)) {

            private$selected_ids <- ids
            private$raw <- matrix(1, nrow = length(ids),
              ncol = length(files))
            rownames(private$raw) <- ids
            colnames(private$raw) <- names(files)
            private$len <- private$raw
            slug_ids <- digest::digest(sort(ids))


          } else {

            if (digest::digest(sort(ids)) != slug_ids) {
              logging::logerror(paste("The are transcript ids mismatch",
                "between samples"))
              stop()
            }

          }

          private$raw[ids, sample_name] <- rhdf5::h5read(fpath, "est_counts")

          if (! with_fixed_length) {
            private$len[ids, sample_name] <- rhdf5::h5read(fpath,
              "aux/eff_lengths")
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
