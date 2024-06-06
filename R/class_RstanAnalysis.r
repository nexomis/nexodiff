#' @include utils.r
#' @include class_PairwiseLFC.r
#' @include class_ExprData.r

NULL
#' Class runing the bayesian model for the distribution parameters inference
#'
#' @description
#'
#' This class generates results using rStan in order to inference diferent
#' distribution parameters from a PairwiseLFC object
#'
#' @export
RStanAnalysis<-R6::R6Class("RStanAnalysis", #nolint
    public <- list(
      #' @description
      #' Initialize `RStanAnalysis` object.
      #'
      #' @param lfc_data PairwiseLFC : log fold change daa
      #' @param option_data is a dataframe with all the parametrization
      #' of the models that must be run. Each row correpond to a model.
      #' * "distribution": 0 (gaussian), 1 (cauchy) or 2 (logistic);
      #'   default (1,1)
      #' * "base_alpha_equals_0": 0 or 1; default (1,1)
      #' * "alpha_shift_equals_0" 0 or 1; default (1,1)
      #' * "base_mu_equals_0": 0 or 1; default (1,0)
      #' * "mu_shift_equals_0": 0 or 1; default (0,0)
      #' * "sigma_ratio_equals_1": 0 or 1; default (0,0)
      #' @param ncpus number of cpus to use for computation (default = 1)
      #' @param seed seed for random computation in rStan (default=1234567)
      #' @param iter number of iteration in the rStan sampling (default=3000)
      #' @param refresh can be used to control how often the progress of the
      #' sampling is reported (i.e. show the progress every refresh
      #' iterations). By default, refresh = max(iter/10, 1). The progress
      #' indicator is turned off if refresh <= 0 (default=0).
      #' @param estimation_approach approach to find one estimates.
      #' @param chains number of chains to run the rStan sampling (default=2)
      #' @param ... other parameters pass to rStan sampling (see the
      #' documentation)
      #'
      #' @return A new `RStanAnalysis` object.
      initialize = function(lfc_data,
        option_data = data.frame("distribution" = rep(1, 2),
          "base_alpha_equals_0" = c(1, 1),
          "alpha_shift_equals_0" = c(1, 1),
          "base_mu_equals_0" = c(1, 0),
          "mu_shift_equals_0" = c(0, 0),
          "sigma_ratio_equals_1" = c(0, 0)),
        ncpus = 1, seed = 1234567,
        iter = 3000, refresh = 0,
        chains = 2, estimation_approach = c("sampling", "optimizing"), ...) {
          private$lfc_data <- lfc_data$get_results()
          private$option_data <- option_data
          private$general_model <- rstan::stan_model(system.file("extdata",
            "rstan_models/general_model.stan", package = "nexodiff"))
          private$estimation_approach <- estimation_approach
          if (ncpus == 1) {
              future::plan(
              strategy = future::sequential()
              )
          } else if (parallelly::supportsMulticore()) {
              future::plan(
              strategy = future::multicore,
              workers = min(ncpus, parallelly::availableCores(omit = 1))
              )
          } else {
              future::plan(
              strategy = future::multisession,
              workers = min(ncpus, parallelly::availableCores(omit = 1))
              )
          }

          if (estimation_approach == "sampling") {
            private$rstan_results <-
              furrr::future_map(unique(private$lfc_data$batch),
                                function(in_batch) {
                    private$main_sampling(setup = private$set_setup(in_batch,
                                                          private$lfc_data),
                                iter = iter,
                                refresh = refresh,
                                cores = min(chains, ncpus),
                                chains = chains, seed = seed, ...)
              })
          }else if (estimation_approach == "optimizing") {
            private$rstan_results <-
              furrr::future_map(unique(private$lfc_data$batch),
                                function(in_batch) {
                    private$main_optimizing(setup = private$set_setup(in_batch,
                                                          private$lfc_data),
                                seed = seed, ...)
              })
          }
          names(private$rstan_results) <- unique(private$lfc_data$batch)
      },

      #' @description
      #' Get the raw results (full fit object) from the bayesian
      #' inference using rStan.
      #'
      #' @return A nested list of lists: one list per batch that contains
      #' one list per model. Each one of the model lists corresponds to the
      #' fit object of rStan.
      get_results = function() {
        private$rstan_results
      },

      #' @description
      #' Get the options used to parametrize the models.
      #'
      #' @return A dataframe with all the parametrization
      #' of the models used to crete the object. Each row correpond to a model
      get_options = function() {
        private$option_data
      },

      #' @description
      #' Get the formated data for the prediction
      #'
      #' @param batchs vector with the batchs names that should
      #' be formatted. If NULL all the batchs are included in the
      #' final output (default=NULL).
      #' @param probs vector of the quantiles vector of
      #' probabilities for which the quantiles should be return.
      #'
      #' @return list of dataframes (one per model in the options
      #' data dataframe) with all the selected batchs as rows and
      #' the the predictors (descriptors of the hyperparameters
      #' of the models) as columns. The default descriptors are:
      #' * mean
      #' * sd
      #' * quantiles associeted whit the selected probabilities (`probs`).
      get_formatted_results = function(batchs = NULL, probs = c(0.25, 0.75)) {
        code_distribution <- list("0" = "gaussian",
          "1" = "cauchy",
          "2" = "logistic")
        code_variables <- list("base_alpha_equals_0" = "alpha_base[1]",
          "alpha_shift_equals_0" = c("alpha_shift[1]", "alpha_2[1]"),
          "base_mu_equals_0" = "mu_base[1]",
          "mu_shift_equals_0" = c("mu_shift[1]", "mu_2[1]"),
          "sigma_ratio_equals_1" = c("sigma_ratio[1]", "sigma_2[1]"))
        list_df <- list()
        if (is.null(batchs)) {
          batchs <- names(private$rstan_results)
        }
        if (private$estimation_approach == "sampling") {
          extract_info <- function(rstan_output, variables, probs) {
            results <- rstan::summary(rstan_output, probs = probs)$summary
            results <- results[variables,
            -which(colnames(results) %in% c("se_mean", "n_eff", "Rhat"))]
            rownames(results) <- stringr::str_remove(rownames(results),
            "\\[1\\]")
            rownames(results) <- stringr::str_replace(rownames(results),
            "2", "treat")
            df_line <- setNames(do.call(c, as.data.frame(results)),
              paste(rownames(results),
              rep(colnames(results), each = nrow(results)), sep = "_"))
            return(df_line)
          }
        }else if (private$estimation_approach == "optimizing") {
          extract_info <- function(rstan_output, variables, ...) {
            results <- rstan_output$par[variables]
            names(results) <- stringr::str_remove(names(results),
            "\\[1\\]")
            names(results) <- stringr::str_replace(names(results),
            "2", "treat")
            df_line <- setNames(do.call(c, as.data.frame(results)),
              paste0(names(results), "_PE"))
            return(df_line)
          }
        }
        for (i in seq_len(nrow(private$option_data))) {
            name <- paste("df",
              code_distribution[[
                as.character(private$option_data$distribution[i])]],
              sum(private$option_data$distribution[1:i] ==
                private$option_data$distribution[i]),
              sep = "_")
            variables <- unlist(
              map(
                colnames(private$option_data)[
                  which(private$option_data[i, ] == 0)],
                ~code_variables[[.x]]
              )
            )
            variables <- c(variables, "sigma_base")
            df <- map_dfr(
              batchs,
              ~extract_info(rstan_output = private$rstan_results[[.x]][[i]],
                variables = variables,
                probs = probs)
            )
            df <- as.data.frame(df)
            rownames(df) <- batchs
            list_df[[name]] <- df
          }
          return(list_df)
      }
    ),
    private <- list(
      lfc_data = NULL,
      option_data = NULL,
      rstan_results = NULL,
      general_model = NULL,
      estimation_approach = NULL,
      # Function that runs the rstan sampling using the compiled general_model
      rstan_sampling = function(setup, distribution, base_alpha_equals_0,
        alpha_shift_equals_0, base_mu_equals_0, mu_shift_equals_0,
        sigma_ratio_equals_1, ...) {
          rstan::sampling(private$general_model,
            data = list(N = setup[[1]], values = setup[[2]],
              group_design = setup[[3]], distribution = distribution,
              base_alpha_equals_0 = base_alpha_equals_0,
              alpha_shift_equals_0 = alpha_shift_equals_0,
              base_mu_equals_0 = base_mu_equals_0,
              mu_shift_equals_0 = mu_shift_equals_0,
              sigma_ratio_equals_1 = sigma_ratio_equals_1),
            ...
          )
      },
      # Function that runs the rstan optimizing (point estimate)
      # using the compiled general_model
      rstan_optimizing = function(setup, distribution, base_alpha_equals_0,
        alpha_shift_equals_0, base_mu_equals_0, mu_shift_equals_0,
        sigma_ratio_equals_1, ...) {
          rstan::optimizing(private$general_model,
            data = list(N = setup[[1]], values = setup[[2]],
              group_design = setup[[3]], distribution = distribution,
              base_alpha_equals_0 = base_alpha_equals_0,
              alpha_shift_equals_0 = alpha_shift_equals_0,
              base_mu_equals_0 = base_mu_equals_0,
              mu_shift_equals_0 = mu_shift_equals_0,
              sigma_ratio_equals_1 = sigma_ratio_equals_1),
            ...
          )
      },
      # Function that calls the rstan_sampling or rstan_optimizing for each one
      # of the models defined in the option_data dataframe
      main_sampling = function(setup, ...) {
        res <- purrr::pmap(
          private$option_data,
          function(distribution, base_alpha_equals_0, alpha_shift_equals_0,
            base_mu_equals_0, mu_shift_equals_0, sigma_ratio_equals_1) {
              private$rstan_sampling(setup, distribution = distribution,
                base_alpha_equals_0 = base_alpha_equals_0,
                alpha_shift_equals_0 = alpha_shift_equals_0,
                base_mu_equals_0 = base_mu_equals_0,
                mu_shift_equals_0 = mu_shift_equals_0,
                sigma_ratio_equals_1 = sigma_ratio_equals_1,
                ...
              )
          }
        )
        return(res)
      },
      main_optimizing = function(setup, ...) {
        res <- purrr::pmap(
          private$option_data,
          function(distribution, base_alpha_equals_0, alpha_shift_equals_0,
            base_mu_equals_0, mu_shift_equals_0, sigma_ratio_equals_1) {
              private$rstan_optimizing(setup, distribution = distribution,
                base_alpha_equals_0 = base_alpha_equals_0,
                alpha_shift_equals_0 = alpha_shift_equals_0,
                base_mu_equals_0 = base_mu_equals_0,
                mu_shift_equals_0 = mu_shift_equals_0,
                sigma_ratio_equals_1 = sigma_ratio_equals_1,
                ...
              )
          }
        )
        return(res)
      },
      # Setup of the input for the bayesian inference
      set_setup = function(in_batch, lfc_data) {
        lfc_nt <-
          lfc_data[lfc_data$batch == in_batch & lfc_data$group == "NT", ]
        lfc_nt <- unlist(lfc_nt$data[[1]], use.names = FALSE)
        lfc_nt <- na.omit(lfc_nt)
        lfc_nt <- lfc_nt[!is.infinite(lfc_nt)]
        lfc_t <- lfc_data[lfc_data$batch == in_batch & lfc_data$group == "T", ]
        lfc_t <- unlist(lfc_t$data[[1]], use.names = FALSE)
        lfc_t <- na.omit(lfc_t)
        lfc_t <- lfc_t[!is.infinite(lfc_t)]
        setup <- list(
          N = (length(lfc_nt) + length(lfc_t)),
          values = c(lfc_nt, lfc_t),
          group_design = c(rep(1, length(lfc_nt)), rep(2, length(lfc_t))))
        return(setup)
      }
    )
)
