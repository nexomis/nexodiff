#' @include utils.r
#' @include class_PairwiseDesignWithAnnotation.r
#' @include class_ExprData.r

NULL

#' Class representing raw Log Fold change results from a pairwaise design
#'
#' @description
#'
#' This class generates raw Log Fold change without statistical inference or
#' tests
#'
#' @export

PairwiseLFC <- R6::R6Class("PairwiseLFC", # nolint
  public = list(
    #' @description
    #' Initialize `PairwiseLFC` object.
    #'
    #' @param expr_data ExprData : expression data
    #' @param ncpus number of cpus to use for computation (default = 1)
    #' @param method_ctrl method to compute the reference expression
    #' - `ctrl_each` : samples are compared to each ctrl
    #' - `ctrl_med` : samples are compared to the median of ctrl
    #' - `ctrl_mean` : samples are compared to the mean of ctrl
    #' - `ctrl_geo` : samples are compared to the geo mean of ctrl
    #' - `ctrl_geo_mod` : samples are compared to the modified geo mean of ctrl
    #' @param compute_chunk a vector of index to select the groups to compute
    #' this option is present to limit compution time.
    #' @return A new `PairwiseLFC` object.
    initialize = function(expr_data, ncpus = 1,
      method_ctrl = "ctrl_med", compute_chunk = NULL
    ) {

      private$expr_data <- expr_data
      private$opts <- list()

      if (! method_ctrl %in%
        c("ctrl_mean", "ctrl_med", "ctrl_each", "ctrl_geo", "ctrl_geo_mod")) {
        logging::logerror("argument `method_ctrl` not recognized")
        stop()
      }

      private$opts$method_ctrl <- method_ctrl

      data_design <- expr_data$get_design()$get_pairwise_design()

      data_design <- data_design[, c("batch", "group")] %>%
          dplyr::distinct()

      if (method_ctrl == "ctrl_each") {
        make_one_pairwise_comp <- function(x, y) {
          private$make_one_pairwise_comp_ctrl_each(x, y)
        }
      } else if (method_ctrl == "ctrl_mean") {
        make_one_pairwise_comp <- function(x, y) {
          private$make_one_pairwise_comp_ctrl_mean(x, y, mean, method_ctrl)
        }
      } else if (method_ctrl == "ctrl_med") {
        make_one_pairwise_comp <- function(x, y) {
          private$make_one_pairwise_comp_ctrl_mean(x, y, median, method_ctrl)
        }
      } else if (method_ctrl == "ctrl_geo") {
        make_one_pairwise_comp <- function(x, y) {
          private$make_one_pairwise_comp_ctrl_mean(x, y,
          function(z) (exp(mean(log(z)))),
          method_ctrl)
        }
      } else if (method_ctrl == "ctrl_geo_mod") {
        make_one_pairwise_comp <- function(x, y) {
          private$make_one_pairwise_comp_ctrl_mean(x, y,
          function(z) (geom_mean_modified(z, 1e-05)),
          method_ctrl)
        }
      }

      if (is.null(compute_chunk)) {
        i_chunks <- seq_len(nrow(data_design))
      } else {
        i_chunks <- intersect(compute_chunk, seq_len(nrow(data_design)))
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

      private$lfc_results <- furrr::future_map2_dfr(
        data_design$batch[i_chunks],
        data_design$group[i_chunks],
        make_one_pairwise_comp
      )

      future::plan(future::sequential)


    },

    #' @description
    #' Get LFC results
    #'
    #' @return A nested dataframe.
    get_results = function() {
      private$lfc_results
    }

  ),
  private <- list(

    expr_data = NULL,
    lfc_results = NULL,
    opts = NULL,

    # function for computing the ratios
    make_one_pairwise_comp_ctrl_each = function(in_batch, in_group) {

      results <- private$expr_data$extract_pairwise_data_with_design(
        in_batch = in_batch, in_group = in_group)
      design_table <- results$design_table
      # Too bad that's recomputed for each group for ctrl
      # However

      norm <- results$raw *
        (results$intra_norm_fact %*% diag(results$inter_norm_fact))

      combinations <- purrr::cross_df(list(
            test_sample = results$test_samples,
            ctrl_sample = results$ctrl_samples
          ),
          .filter = function(x, y) x == y
        ) %>%
          dplyr::mutate(
            key = paste0(
              pmin(.data$test_sample, .data$ctrl_sample),
              pmax(.data$test_sample, .data$ctrl_sample),
            sep = "")
        ) %>%
          dplyr::distinct(key, .keep_all = TRUE)

      lfc_func <- function(x, y) {
        df <- data.frame(value = log2(norm[, x] / norm[, y]))
        row.names(df) <- row.names(norm)
        names(df) <- c(paste(x, y, sep = "_vs_"))
        return(df)
      }

      data_res <- purrr::map2_dfc(
        combinations$test_sample,
        combinations$ctrl_sample,
        lfc_func
      )

      data_res$batch <- in_batch
      data_res$group <- in_group

      #return nested
      data_res %>%
        dplyr::group_by(.data$batch, .data$group) %>%
        tidyr::nest()
    },

    make_one_pairwise_comp_ctrl_mean = function(in_batch, in_group,
      mean_funct, method_ctrl) {

      results <- private$expr_data$extract_pairwise_data_with_design(
        in_batch = in_batch, in_group = in_group)
      design_table <- results$design_table
      # Too bad that's recomputed for each group for ctrl
      norm <- results$raw *
        (results$intra_norm_fact %*% diag(results$inter_norm_fact))

      ctrl_data <- unlist(lapply(seq(1,nrow(norm)), function(i){
        mean_funct(norm[i, results$ctrl_samples])
      }))
      lfc_func <- function(x) {
        df <- data.frame(value = log2(norm[, x] / ctrl_data))
        names(df) <- c(paste(x, method_ctrl, sep = "_vs_"))
        row.names(df) <- row.names(norm)
        return(df)
      }

      data_res <- purrr::map_dfc(
        results$test_sample,
        lfc_func
      )

      data_res$batch <- in_batch
      data_res$group <- in_group

      #return nested
      data_res %>%
        dplyr::group_by(.data$batch, .data$group) %>%
        tidyr::nest()

    }
  )
)