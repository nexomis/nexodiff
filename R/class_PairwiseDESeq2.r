#' @include utils.r
#' @include class_PairwiseDesignWithAnnotation.r
#' @include class_ExprData.r
#' @include class_PairwiseComp.r

NULL

#' @title DESeq2 Results from Pairwise Design
#' @description A class representing DESeq2 results from a pairwise design.
#'
#' @details
#' This class generates results from expression data associated with a pairwise
#' design using the DESeq2 package.
#'
#' This class inherits from the "PairwiseComp" class.
#' @export

PairwiseDESeq2 <- R6::R6Class("PairwiseDESeq2", # nolint
  inherit = PairwiseComp,
  public <- list(

    #' @description
    #' Initialize `PairwiseDESeq2` object.
    #'
    #' @param expr_data ExprData : expression data
    #' @param add_vars Character vactor with additional variable name in deisgn.
    #' Those given variable will be add to the analysis formula. However the
    #' fold-change in results will be the coefficient from the variable `group`
    #' @param fit_formula Character with the formula which must contains group.
    #' 
    #' @param only_paired Wether to apply a filtering to extract only the
    #' samples which are fully paired given the concatenated values in add_vars.
    #' @param auto_filtering whether or not to perform independent filtering.
    #' @param min_base_mean minimum value of base mean for filtering. Apply if
    #' not null and is priority before auto_filtering. default is null.
    #' @param ncpus number of cpus to use for computation (default = 1)
    #' @param seed seed for random computation in DESeq2
    #' @param force_deseq2_norm renorm with deseq2 in place of internorm
    #' @param shrinkage The parameter for shrinkage type, if null, no shrinkage
    #' is applied. see DESeq2 documentation for shrinkage types.
    #'
    #' @return A new `PairWiseCompDESeq2` object.
    initialize = function(expr_data, add_vars = c(), auto_filtering = TRUE,
      ncpus = 1, seed = 1234567, force_deseq2_norm = FALSE, shrinkage = NULL,
      min_base_mean = NULL, only_paired = FALSE, fit_formula = NULL
    ) {
      private$expr_data <- expr_data
      private$opts <- list()
      private$opts$add_vars <- add_vars
      if (length(add_vars) == 0) {
        private$opts$only_paired <- FALSE
      } else {
        private$opts$only_paired <- only_paired
      }
      private$opts$auto_filtering <- auto_filtering
      private$opts$force_deseq2_norm <- force_deseq2_norm
      private$opts$min_base_mean <- min_base_mean
      private$opts$shrinkage <- shrinkage
      private$opts$fit_formula <- fit_formula

      data_design <- expr_data$get_design()$get_pairwise_design()

      data_design <-
        subset(
          data_design,
          ! data_design$ctrl
        )[, c("batch", "group")] %>%
          dplyr::distinct()

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

      private$results <- furrr::future_map2_dfr(
        data_design$batch,
        data_design$group,
        private$make_one_pairwise_comp,
        .options = furrr::furrr_options(seed = seed)
      )

      future::plan(future::sequential)
    }

  ),
  private <- list(
    # function for one pairwise comparison with DESeq2
    make_one_pairwise_comp = function(in_batch, in_group) {

      add_vars <- private$opts$add_vars
      auto_filtering <- private$opts$auto_filtering
      shrinkage <- private$opts$shrinkage
      only_paired <- private$opts$only_paired

      results <- private$expr_data$extract_pairwise_data_with_design(
        in_batch = in_batch, in_group = in_group)

      groupvars <- c("group", add_vars)
      design_deseq2 <- results$design_table

      if (only_paired) {
        design_filter_table <- data.frame(
          sample = design_deseq2$sample,
          group = design_deseq2$group,
          replicate = ""
        )
      }

      design_deseq2[, "group"] <- as.factor(design_deseq2[, "group"])

      for (var in add_vars) {
        if (! var %in% colnames(design_deseq2)) {
          logging::logerror("`add_vars` not recognized in design")
          stop()
        }
        if (only_paired) {
          design_filter_table[, "replicate"] <- paste0(
            design_filter_table[, "replicate"],
            as.character(design_deseq2[, var])
          )
        }
        design_deseq2[, var] <- as.factor(design_deseq2[, var])
      }

      design_deseq2 <- dplyr::select(
        design_deseq2, dplyr::all_of(groupvars))

      count_deseq2 <- round(results$raw)
      mode(count_deseq2) <- "integer"

      if (only_paired) {
        ugroups <- unique(design_filter_table$group)
        assertthat::assert_that(length(ugroups) == 2)
        keep_replicates <- intersect(
          unique(subset(design_filter_table, group == ugroups[1])$replicate),
          unique(subset(design_filter_table, group == ugroups[2])$replicate)
        )
        assertthat::assert_that(length(keep_replicates) > 1)
        keep_samples <-
          subset(design_filter_table, replicate %in% keep_replicates)$sample
        assertthat::assert_that(length(keep_replicates) > 1)
        design_deseq2 <- design_deseq2[keep_samples, ]
        count_deseq2 <- count_deseq2[, keep_samples]
      }

      if (is.null(private$opts$fit_formula)) {
        design_formula <- as.formula(
          paste("~", paste(groupvars, collapse = " + "), sep = " ")
        )
      } else {
        design_formula <- as.formula(private$opts$fit_formula)
      }

      # 1. Remove spaces
      formula_nospace <- gsub(
        " ", "",
        paste(as.character(design_formula), collapse = "")
      )

      # 2. Check starts with "~group"
      if (!startsWith(formula_nospace, "~group")) {
        emsg <- "The formula does not start with '~group'."
        logging::logerror(emsg)
        stop()
      }

      # Extract the right side of the formula
      right_side <- gsub("^~group", "", formula_nospace)
      if (right_side != "") {
        elements <- unlist(strsplit(right_side, "\\+"))
        elements <- elements[2:length(elements)]
        for (element in elements) {
          # 3. Separate value with and without ":"
          if (grepl(":", element)) {
            # 5. Split each value with ":" by ":" and check
            parts <- unlist(strsplit(element, ":"))
            if (! "group" %in% parts) {
              emsg <-
                paste("Invalid interaction term (group not present):", element)
              logging::logerror(emsg)
              stop()
            }
            var_name <- parts[parts != "group"]
            if (length(var_name) != 1) {
              emsg <- paste(
                "Invalid interaction term (more than 1 interaction):",
                element
              )
              logging::logerror(emsg)
              stop()
            }
            if (!(var_name %in% add_vars)) {
              emsg <- paste(
                "Invalid interaction term (undeclared var):",
                element
              )
              logging::logerror(emsg)
              stop()
            }
          } else {
            # 4. Check that values without ":" are in add_vars
            if (!(element %in% add_vars)) {
              emsg <- paste("Variable", element, "not found in add_vars")
              logging::logerror(emsg)
              stop()
            }
          }
        }
      }

      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = count_deseq2,
        colData = design_deseq2,
        design = design_formula
      )
      if (private$opts$force_deseq2_norm) {
        norm_fact <- results$len
        # DESeq2 authors recommend providing a matrix with row-wise geometric
        # means of 1, so that the mean of normalized counts for a gene is close
        # to the mean of the unnormalized counts.
        # TODO: Check if there are only non-zero value ? otherwise the geometric
        # means can be problematic
        if (only_paired) {
          norm_fact <- norm_fact[, keep_samples]
        }
        norm_fact <- norm_fact / exp(rowMeans(log(norm_fact)))
        dds <- DESeq2::estimateSizeFactors(dds, normMatrix = norm_fact)
      } else {
        norm_fact <- results$intra_norm_fact %*% diag(results$inter_norm_fact)
        colnames(norm_fact) <- colnames(results$intra_norm_fact)
        if (only_paired) {
          norm_fact <- norm_fact[, keep_samples]
        }
        norm_fact <- norm_fact / exp(rowMeans(log(norm_fact)))
        DESeq2::normalizationFactors(dds) <- 1 / norm_fact
      }

      design_n_ctrl <- nrow(subset(design_deseq2, group == results$ctrl_group))
      design_n_test <- nrow(design_deseq2) - design_n_ctrl

      dds$group <- relevel(dds$group, ref = results$ctrl_group)
      # fit set to local. If not, in case of no convergence in other algorithms,
      # there is anyway a fallback to "local".
      dds <- DESeq2::DESeq(dds, fitType = "local", quiet = TRUE)
      # Get the result table before to make some option possible with LFC
      # shrinkage like independentFiltering
      res <- DESeq2::results(
        dds,
        name = DESeq2::resultsNames(dds)[2],
        independentFiltering = auto_filtering,
      )
      if (is.null(shrinkage)) {
        res_shrink <- res
      } else {
        res_shrink <- DESeq2::lfcShrink(
          dds,
          coef = 2,
          res = res,
          type = shrinkage,
          quiet = TRUE
        )
      }
      data_res <- dplyr::arrange(
        as.data.frame(res_shrink),
        .data$pvalue
      )[, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

      data_res$status <- "analyzed"

      if (! is.null(private$opts$min_base_mean)) {
        data_res$padj <- NA
        select_rows <- (data_res$baseMean > private$opts$min_base_mean) &
          (! is.na(data_res$pvalue))
        data_res$padj[select_rows] <- p.adjust(data_res$pvalue[select_rows],
          method = "BH")
      }

      data_res[is.na(data_res$padj), "status"] <- "filtered"
      data_res[is.na(data_res$pvalue), "status"] <- "outlier"
      data_res[data_res$baseMean == 0, "status"] <- "undetected"
      t_levels <- c("undetected", "filtered", "outlier", "analyzed")
      data_res$status <- factor(data_res$status, levels = t_levels)

      data_res$tag_id <- row.names(data_res)
      data_res$batch <- in_batch
      data_res$group <- in_group

      #return nested
      data_res %>%
        dplyr::group_by(.data$batch, .data$group) %>%
        tidyr::nest() %>%
        mutate(
          n_test = design_n_test,
          n_ctrl = design_n_ctrl
        )

    }

  )
)
