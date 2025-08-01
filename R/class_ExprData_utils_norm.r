#' @include nexodiff-package.R
#' @include utils.r
NULL

#' Helper function to set mean function based on method name
#'
#' @param method_name Method name for mean calculation
#' @return Function for computing mean
set_mean_function <- function(method_name) {
  rfun <- NULL
  if (method_name == "median") {
    rfun <- stats::median
  } else if (method_name == "geometric") {
    rfun <- function(x) exp(mean(log(x)))
  } else if (method_name == "nz.geometric") {
    rfun <- function(x) exp(mean(log(x[x != 0])))
  } else if (method_name == "mod.geometric") {
    rfun <- function(x) exp(mean(log(x + 1e-05)))
  } else if (method_name == "arithmetic") {
    rfun <- mean
  } else {
    logging::logerror("Unrecognized mean function method")
    stop()
  }
  rfun
}

#' Helper function for median normalization factor calculation
#'
#' @param norm normalized matrix
#' @param sample_ref reference sample values
#' @return normalization factors
calc_norm_fact_median <- function(norm, sample_ref) {
  # filter out when sample_ref = 0
  keep_row <- sample_ref != 0
  norm <- norm[keep_row, ] / sample_ref[keep_row]
  norm_fact <- 1 / matrixStats::colMedians(norm)
  names(norm_fact) <- colnames(norm)
  norm_fact
}

#' Helper function to get reference for normalization
#'
#' @param norm normalized matrix
#' @param ctrl_samples control sample names
#' @param norm_ref reference specification
#' @param norm_ref_mean_fun mean function for reference calculation
#' @return reference values
get_reference <- function(norm, ctrl_samples, norm_ref, norm_ref_mean_fun) {
  if (all(norm_ref == "all")) {
    norm_ref_samples <- colnames(norm)
  } else if (all(norm_ref == "ctrl")) {
    norm_ref_samples <- ctrl_samples
  } else if (all(norm_ref %in% colnames(norm)) |
    all(norm_ref %in% seq_len(ncol(norm)))
  ) {
    norm_ref_samples <- norm_ref
  } else {
    logging::logerror("`norm_ref` argument not recognized")
    stop()
  }
  norm_ref_values <- unlist(lapply(seq(1, nrow(norm)), function(i) {
    norm_ref_mean_fun(norm[i, norm_ref_samples])
  }))
  norm_ref_values
}

#' Helper function to get trimmed matrix
#'
#' @param raw raw count matrix
#' @param norm_fact normalization factors
#' @param a_trim_norm whether to use normalized values for trimming
#' @param a_trim_value threshold for trimming
#' @param a_trim_mean_fun mean function for trimming
#' @param replace_zero_by value to replace zeros
#' @return trimmed matrix
get_trimmed_mat <- function(
  raw, norm_fact, a_trim_norm, a_trim_value,
  a_trim_mean_fun, replace_zero_by
) {

  trimmed_mat <- raw * norm_fact

  if (a_trim_norm) {
    base_mat <- trimmed_mat
  } else {
    base_mat <- raw
  }
  base_mat_mean <- unlist(lapply(seq(1, nrow(base_mat)), function(i) {
    a_trim_mean_fun(base_mat[i, ])
  }))
  trimmed_mat[trimmed_mat == 0] <- replace_zero_by
  trimmed_mat[base_mat_mean > a_trim_value, ]
}

#' Helper function for ExprData$compute_and_set_intra_norm_fact
#'
#' @param method see ExprData$compute_and_set_intra_norm_fact
#' @param selected_ids private$selected_ids (see ExprData)
#' @param raw_matrix raw count matrix
#' @param len_matrix length matrix
#' @return intra normalization factor matrix
compute_intra_norm_factor <- function(
  method, selected_ids, raw_matrix, len_matrix
) {
  intra_norm_fact <- matrix(rep(1, length(raw_matrix)),
    nrow = nrow(raw_matrix),
    ncol = ncol(raw_matrix),
    dimnames = list(
      row.names(raw_matrix),
      colnames(raw_matrix)
    )
  )

  if (method == "tpm") {
    intra_norm_fact[selected_ids, ] <-
      sweep(1e+06 / len_matrix,
        2,
        colSums(raw_matrix / len_matrix),
        FUN = "/"
      )
  } else if (method == "fpk") {
    intra_norm_fact[selected_ids, ] <-
      1e+03 / len_matrix
  } else if (method == "fpm") {
    intra_norm_fact[selected_ids, ] <-
      sweep(
        matrix(1e+06,
          nrow = length(selected_ids),
          ncol = ncol(raw_matrix),
          dimnames = list(
            selected_ids,
            colnames(raw_matrix)
          )
        ),
        2,
        colSums(raw_matrix),
        FUN = "/"
      )
  } else if (method == "fpkm") {
    intra_norm_fact[selected_ids, ] <-
      sweep(1e+09 / len_matrix,
        2,
        colSums(raw_matrix),
        FUN = "/"
      )
  } else if (method == "none") {
    intra_norm_fact[selected_ids, ][] <- 1
  } else {
    logging::logdebug("Unrecognized normalization method")
    stop()
  }

  # There will be NaN value where colSums(raw_matrix) or
  # colSums(raw_matrix/len_matrix) is 0
  # We are going to replace those NaN value by 1
  intra_norm_fact[is.nan(intra_norm_fact)] <- 1

  intra_norm_fact
}

#' Helper function to check args
#' @param method see ExprData$compute_and_set_intra_norm_fact
#' @param norm_scale see ExprData$compute_and_set_intra_norm_fact
#' @param norm_by see ExprData$compute_and_set_intra_norm_fact
#' @param norm_ref see ExprData$compute_and_set_intra_norm_fact
#' @param norm_ref_mean see ExprData$compute_and_set_intra_norm_fact
#' @param m_trim_prop see ExprData$compute_and_set_intra_norm_fact
#' @param m_trim_mean see ExprData$compute_and_set_intra_norm_fact
#' @param a_trim_mean see ExprData$compute_and_set_intra_norm_fact
#' @param a_trim_norm see ExprData$compute_and_set_intra_norm_fact
#' @param ncpus see ExprData$compute_and_set_intra_norm_fact
#' @param a_trim_value see ExprData$compute_and_set_intra_norm_fact
#' @param replace_zero_by see ExprData$compute_and_set_intra_norm_fact
inter_norm_check <- function(
  method, norm_scale, norm_by, norm_ref, norm_ref_mean, m_trim_prop,
  m_trim_mean, a_trim_value, a_trim_norm, a_trim_mean, replace_zero_by, ncpus
){

  if (method == "tmm" && m_trim_mean == "median") {
    logging::logwarning(paste(
      "tmm method with median as mean function...",
      "Is that not simply the median method ?"
    ))
  }

  if (norm_scale == "group") {
    if (norm_ref != "ctrl") {
      logging::loginfo(
        "For group-scale normalization, `norm_ref` is forced to 'ctrl'."
      )
      norm_ref <- "ctrl"
    }
  }

  if (method == "tmm" && m_trim_prop == 0) {
    logging::logerror("`m_trim_prop` cannot be 0 with `method` tmm")
    stop()
  }

  if (m_trim_prop <= 0 || m_trim_prop >= 0.5) {
    logging::logerror("`m_trim_prop must` be between 0 and 0.5 excluded")
    stop()
  }

  if (norm_by == "group") {
    logging::logerror("`norm_by` group is not yet implemented")
    stop()
  }

  if (norm_ref == "ctrl" && norm_scale == "design") {
    logging::logerror(
      "`norm_ref` ctrl is not compatible with `norm_scale` design"
    )
    stop()
  }

  if (ncpus > parallelly::availableCores() - 1) {
    logging::logwarn("The number of core asked was more than available")
  }

}

#' Helper function for ExprData$compute_norm_fact
#'
#' @param selected_ids private$selected_ids (see ExprData)
#' @param design private$design (see ExprData)
#' @param intra_norm_fact private$intra_norm_fact (see ExprData)
#' @param inter_norm_fact private$inter_norm_fact (see ExprData)
#' @param inter_norm_fact_opts private$inter_norm_fact_opts (see ExprData)
#' @param in_batch see ExprData$compute_norm_fact
#' @param in_group see ExprData$compute_norm_fact
#' @param inter_norm see ExprData$compute_norm_fact
#' @param intra_norm see ExprData$compute_norm_fact
#' @return normalization factor matrix
compute_norm_fact_helper <- function(
  selected_ids, design, intra_norm_fact, inter_norm_fact,
  inter_norm_fact_opts, in_batch = NULL, in_group = NULL,
  inter_norm = FALSE, intra_norm = TRUE
) {
  row_ids <- selected_ids
  col_ids <- design$extract_sample_names(in_batch, in_group)

  if (intra_norm) {
    norm_fact <- intra_norm_fact[row_ids, col_ids]
  } else {
    norm_fact <- NULL
  }

  norm_fact_inter <- NULL

  if (inter_norm) {
    if (!is.null(inter_norm_fact)) {
      if (inter_norm_fact_opts$norm_scale == "design") {
        norm_fact_inter <- inter_norm_fact[col_ids]
      } else {
        if (is.null(in_batch)) {
          logging::logerror(
            "`in_batch` cannot be null when norm scale is not design"
          )
          stop()
        }
        if (inter_norm_fact_opts$norm_scale == "batch") {
          norm_fact_inter <- inter_norm_fact[[in_batch]][col_ids]
        } else {
          if (is.null(in_group)) {
            logging::logerror(
              "`in_group` cannot be null when norm scale is group"
            )
            stop()
          }
          norm_fact_inter <-
            inter_norm_fact[[in_batch]][[in_group]][col_ids]
        }
      }
    } else {
      norm_fact_inter <- rep(1, length(col_ids))
      names(norm_fact_inter) <- col_ids
    }
  }
  if (!is.null(norm_fact_inter)) {
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
}
