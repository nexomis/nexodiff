#' @include utils.r
#' @include class_PairwiseDesignWithAnnotation.r
#' @include class_ExprData.r
#' @include class_PairwiseComp.r

NULL

#' @title Dearseq Results from Pairwise Design
#' @description A class representing DESeq2 results from a pairwise design.
#'
#' @details
#' This class generates results from expression data associated with a pairwise
#' design using the dearseq package.
#'
#' This class inherits from the "PairwiseComp" class.
#' @export

PairwiseDearseq <- R6::R6Class("PairwiseDearseq", # nolint
  inherit = PairwiseComp,
  public <- list(

    #' @description
    #' Initialize `PairwiseDearseq` object.
    #'
    #' @param expr_data ExprData : expression data
    #' @param add_vars Character vactor with additional variable name in deisgn.
    #' Those given variable will be add to the analysis formula. However the
    #' fold-change in results will be the coefficient from the variable `group`
    #' @param auto_filtering whether or not to perform independent filtering.
    #' @param ncpus number of cpus to use for computation (default = 1)
    #' @param seed seed for random computation in Dearseq
    #' @param force_deseq2_norm renorm with deseq2 in place of internorm
    #'
    #' @return A new `PairWiseCompDESeq2` object.
    initialize = function(expr_data, add_vars = c(), auto_filtering = TRUE,
      ncpus = 1, seed = 1234567, force_deseq2_norm = FALSE
    ) {
      private$expr_data <- expr_data
      private$opts <- list()
      private$opts$add_vars <- add_vars
      private$opts$auto_filtering <- auto_filtering

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

      results <- private$expr_data$extract_pairwise_data_with_design(
        in_batch = in_batch, in_group = in_group)

      groupvars <- c("group", add_vars)
      design_deseq2 <- results$design_table

      for (var in groupvars) {
          if (! var %in% colnames(design_deseq2)) {
            logging::logerror("`add_vars` not recognized in design")
            stop()
          }
          design_deseq2[, var] <- as.factor(design_deseq2[, var])
      }

      design_deseq2 <- dplyr::select(
        design_deseq2, dplyr::all_of(groupvars))

      counts <- round(results$raw)
      mode(counts) <- "integer"

      design_formula <- as.formula(
          paste("~", paste(groupvars, collapse = " + "), sep = " ")
      )

      dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = counts,
          colData = design_deseq2,
          design = design_formula
      )
      if (private$opts$force_deseq2_norm) {
        norm_fact <- results$len
        # DESeq2 authors recommend providing a matrix with row-wise geometric
        # means of 1, so that the mean of normalized counts for a gene is close
        # to the mean of the unnormalized counts.
        norm_fact <- norm_fact / exp(rowMeans(log(norm_fact)))
        dds <- DESeq2::estimateSizeFactors(dds, normMatrix = norm_fact)
      } else {
        norm_fact <- results$intra_norm_fact %*% diag(results$inter_norm_fact)
        norm_fact <- norm_fact / exp(rowMeans(log(norm_fact)))
        DESeq2::normalizationFactors(dds) <- 1 / norm_fact
      }
      dds$group <- relevel(dds$group, ref = results$ctrl_group)
      # fit set to local. If not, in case of no convergence in other algorithms,
      # there is anyway a fallback to "local".
      dds <- dearseq::dear_seq(
        object = dds,
        covariates = add_vars,
        variables2test = "group")
      # Get the result table before to make some option possible with LFC
      # shrinkage like independentFiltering
      res <- DESeq2::results(
        dds,
        name = DESeq2::resultsNames(dds)[2],
        independentFiltering = auto_filtering,
        # sample_group = ?,
        preprocessed = TRUE
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
        tidyr::nest()

    }

  )
)
