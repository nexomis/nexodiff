#' @include nexodiff-package.R
#' @include utils.r
NULL

#' Helper function for median normalization factor calculation
#'
#' @param tgt_vector target expression vector
#' @param ref_vector reference expression vector
#' @param mean_fun mean function
#' @return normalization factors
#' @keywords internal
calc_norm_fac <- function(tgt_vector, ref_vector, norm_mean_fun,
                          trim_extreme) {
  keep_row <- !(tgt_vector == 0 & ref_vector == 0)
  ratios <- tgt_vector[keep_row] / ref_vector[keep_row]

  if (trim_extreme && (any(is.infinite(ratios)) || any(ratios == 0))) {
    n_inf <- sum(is.infinite(ratios))
    n_zero <- sum(ratios == 0)
    n_trim <- max(n_inf, n_zero)

    if (n_trim > 0 && length(ratios) > 2 * n_trim) {
      sorted_ratios <- sort(ratios)
      ratios <- sorted_ratios[(n_trim + 1):(length(sorted_ratios) - n_trim)]
    } else if (n_trim > 0) {
      logging::logerror(
        "Cannot trim extreme values: insufficient data points"
      )
      stop(1)
    }
  }

  1 / norm_mean_fun(ratios)
}

#' Main helper function for normalization factor calculation
#' @param mat expression matrix
#' @param ref_mat reference expression matrix
#' @param ref_mean_fun mean function for reference
#' @param norm_mean_fun mean function for normalization
#' @param tgt_mean_fun mean function for target
#' @param a_mean_fun mean function for A vector
#' @param a_trim_value trim value for A vector
#' @param m_trim_prop trim proportion for M vector
#' @param grouped_samples grouped samples
#' @keywords internal
main_calc_norm_fac <- function(
  mat, ref_mat, ref_mean_fun, norm_mean_fun, tgt_mean_fun,
  a_mean_fun, a_trim_value, m_trim_prop, grouped_samples, trim_extreme
) {
  # get reference
  assert_that(nrow(mat) == nrow(ref_mat))
  ref_vector <- apply(ref_mat, 1, ref_mean_fun)
  # get ref trimmed ids
  unlist(purrr::map(
    seq_along(grouped_samples),
    function(i) {
      tgt_mat <- mat[, grouped_samples[[i]], drop = FALSE]
      if (length(grouped_samples[[i]]) == 1) {
        tgt_vec <- tgt_mat[, 1]
      } else {
        tgt_vec <- apply(tgt_mat, 1, tgt_mean_fun)
      }
      trimmed_ref <- ref_vector
      if (a_trim_value != 0 || m_trim_prop != 0) {
        trimmed_ids <- get_trimmed_ids(
          ref_vector, tgt_vec, a_mean_fun, a_trim_value, m_trim_prop
        )
        tgt_vec <- tgt_vec[trimmed_ids]
        trimmed_ref <- trimmed_ref[trimmed_ids]
      }
      setNames(rep(calc_norm_fac(
        tgt_vec, trimmed_ref, norm_mean_fun, trim_extreme
      ), length(grouped_samples[[i]])), grouped_samples[[i]])
    }
  ))
}

#' Helper function to get reference ids
#' @param ids possible sample ids
#' @param ctrl_samples control sample ids
#' @param ref_type type of reference
#' - "all": use all samples as reference
#' - "ctrl": use control samples as reference
#' - "specified": use specified samples as reference
#' @param ref_samples specified reference samples
#' @keywords internal
get_reference_ids <- function(
  ids,
  ctrl_samples,
  ref_type = c("all", "ctrl", "specified"),
  ref_samples = NULL
) {
  ref_type <- match.arg(ref_type)
  switch(
    ref_type,
    all = ids,
    ctrl = ctrl_samples,
    specified = {
      assert_that(!is.null(ref_samples))
      if (is.numeric(ref_samples)) {
        assert_that(all(ref_samples %in% seq_along(ids)))
        ref_samples <- ids[ref_samples]
      } else {
        assert_that(all(ref_samples %in% ids))
        ref_samples
      }
    }
  )
}

#' Helper function to get trimmed ids
#' @param trim_ref reference expression vector
#' @param trim_grp expression vector
#' @param mean_fun mean function to compute average
#' @param a_trim_value trim value
#' @param m_trim_prop trim proportion top and bottom
#' @return trimmed ids
#' @keywords internal
get_trimmed_ids <- function(
  trim_ref, trim_grp, mean_fun, a_trim_value, m_trim_prop
) {
  assert_that(! is.null(names(trim_ref)))
  if (a_trim_value > 0) {
    a_vector <- apply(
      matrix(c(trim_ref, trim_grp), byrow = FALSE, ncol = 2,
             dimnames = list(names(trim_ref), c("ref", "grp"))), 1, mean_fun
    )
    ids <- names(a_vector[a_vector > a_trim_value])
  } else {
    ids <- names(trim_ref)
  }
  if (m_trim_prop > 0) {
    assert_that(m_trim_prop < 0.5)
    ratios <- trim_grp[ids] / trim_ref[ids]
    n_trim <- floor(length(ids) * m_trim_prop)
    ordered_ids <-
      names(ratios)[order(ratios)][(n_trim + 1):(length(ratios) - n_trim)]
    ids <- ids[ids %in% ordered_ids]
  }
  ids
}

#' Helper function to do inter normalization
#' @param raw raw matrix
#' @param norm_fact intra norm fact matrix
#' @param ref_type reference type
#' @param ref_samples reference samples
#' @param design private$design
#' @param norm_scale scale of normalization
#' @param norm_by normalization by group or sample
#' @param ref_mean mean function for reference
#' @param norm_mean mean function for normalization
#' @param tgt_mean mean function for target
#' @param a_mean mean function for A
#' @param a_trim_value trim value for A
#' @param m_trim_prop trim proportion for M
#' @return inter norm fact matrix
#' @keywords internal
compute_inter_norm <- function(
  raw, norm_fact, ref_type, ref_samples, design, norm_scale, norm_by, ref_mean,
  norm_mean, tgt_mean, a_mean, a_trim_value, m_trim_prop, trim_extreme = FALSE
) {
  ref_mean_fun <- set_mean_function(ref_mean)
  norm_mean_fun <- set_mean_function(norm_mean)
  tgt_mean_fun <- set_mean_function(tgt_mean)
  a_mean_fun <- set_mean_function(a_mean)
  assert_that(norm_scale %in% c("design", "batch", "group"))
  assert_that(norm_by %in% c("sample", "group"))
  batches <- design$list_batches()
  batch2ctrl <- design$find_control_group_per_batches()
  mat <- norm_fact * raw
  switch(
    norm_scale,
    design = {
      ctrl_ids <- unique(unlist(purrr::map(
        batches,
        ~ design$extract_sample_names(.x, batch2ctrl[[.x]])
      )))
      ref_ids <- get_reference_ids(colnames(mat), ctrl_ids, ref_type,
                                   ref_samples)
      data_design <- unique(design$get_pairwise_design()[c("batch", "group")])
      grouped_samples <-  switch(
        norm_by,
        group = purrr::map2(
          data_design$batch,
          data_design$group,
          ~ design$extract_sample_names(.x, .y)
        ),
        sample = as.list(colnames(mat))
      )
      main_calc_norm_fac(
        mat, mat[, ref_ids, drop = FALSE], ref_mean_fun, norm_mean_fun,
        tgt_mean_fun, a_mean_fun, a_trim_value, m_trim_prop, grouped_samples,
        trim_extreme
      )
    },
    batch = {
      batch2groups <- design$list_groups_per_batches(include_ctrl = TRUE)
      setNames(purrr::map(
        batches,
        function(batch) {
          ctrl_grp <- batch2ctrl[[batch]]
          ctrl_ids <- design$extract_sample_names(
            in_batch = batch, in_group = ctrl_grp
          )
          tgt_ids <- design$extract_sample_names(
            in_batch = batch
          )
          ref_ids <- get_reference_ids(tgt_ids, ctrl_ids, ref_type, ref_samples)
          grouped_samples <- switch(
            norm_by,
            group = purrr::map(
              batch2groups[[batch]],
              ~ design$extract_sample_names(batch, .x)
            ),
            sample = as.list(tgt_ids)
          )
          main_calc_norm_fac(
            mat[, tgt_ids, drop = FALSE], mat[, ref_ids, drop = FALSE],
            ref_mean_fun, norm_mean_fun, tgt_mean_fun, a_mean_fun,
            a_trim_value, m_trim_prop, grouped_samples, trim_extreme
          )
        }
      ), batches)
    },
    group =  {
      batch2groups <- design$list_groups_per_batches(include_ctrl = TRUE)
      setNames(purrr::map(
        batches,
        function(batch) {
          ctrl_grp <- batch2ctrl[[batch]]
          ctrl_ids <- design$extract_sample_names(
            in_batch = batch, in_group = ctrl_grp
          )
          setNames(purrr::map(
            batch2groups[[batch]],
            function(group) {
              grp_ids <- design$extract_sample_names(
                in_batch = batch, in_group = group
              )
              tgt_ids <- unique(c(grp_ids, ctrl_ids))
              ref_ids <- get_reference_ids(tgt_ids, ctrl_ids, ref_type,
                                           ref_samples)
              grouped_samples <- switch(
                norm_by,
                group = list(
                  grp_ids, ctrl_ids
                ),
                sample = as.list(tgt_ids)
              )
              main_calc_norm_fac(
                mat[, tgt_ids, drop = FALSE], mat[, ref_ids, drop = FALSE],
                ref_mean_fun, norm_mean_fun, tgt_mean_fun, a_mean_fun,
                a_trim_value, m_trim_prop, grouped_samples, trim_extreme
              )
            }
          ), batch2groups[[batch]])
        }
      ), batches)
    }
  )
}



#' Helper function for ExprData$compute_and_set_intra_norm_fact
#'
#' @param method see ExprData$compute_and_set_intra_norm_fact
#' @param selected_ids private$selected_ids (see ExprData)
#' @param raw_matrix raw count matrix
#' @param len_matrix length matrix
#' @return intra normalization factor matrix
#' @keywords internal
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
    logging::logerror("Unrecognized normalization method")
    stop(1)
  }

  # There will be NaN value where colSums(raw_matrix) or
  # colSums(raw_matrix/len_matrix) is 0
  # We are going to replace those NaN value by 1
  intra_norm_fact[is.nan(intra_norm_fact)] <- 1

  intra_norm_fact
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
#' @param include_ctrl if TRUE and if `in_batch` and `in_group` are
#' specified, control samples from the specified batch will be included.
#' @param rescale_inter_norm Rescale inter-sample normalization factors so
#' that their geometric mean is 1. Default is TRUE.
#' @return normalization factor matrix
#' @keywords internal
compute_norm_fact_helper <- function(
  selected_ids, design, intra_norm_fact, inter_norm_fact,
  inter_norm_fact_opts, in_batch = NULL, in_group = NULL,
  inter_norm = FALSE, intra_norm = TRUE, include_ctrl = FALSE,
  rescale_inter_norm = TRUE
) {

  assert_that((!inter_norm) || (inter_norm && intra_norm),
    msg = "`intra_norm` cannot be FALSE if `inter_norm` is TRUE"
  )

  row_ids <- selected_ids
  col_ids <- design$extract_sample_names(in_batch, in_group)

  if (include_ctrl && !is.null(in_batch) && !is.null(in_group)) {
    batch2ctrl <- design$find_control_group_per_batches()
    ctrl_group <- batch2ctrl[[in_batch]]
    if (!is.null(ctrl_group)) {
      ctrl_ids <- design$extract_sample_names(in_batch, ctrl_group)
      col_ids <- unique(c(col_ids, ctrl_ids))
    }
  }

  if (intra_norm) {
    norm_fact <- intra_norm_fact[row_ids, col_ids, drop = FALSE]
  } else {
    norm_fact <- matrix(
      1, dimnames = list(row_ids, col_ids), nrow = length(row_ids),
      ncol = length(col_ids)
    )
  }

  norm_fact_inter <- NULL

  if (inter_norm) {
    assert_that(! is.null(inter_norm_fact))
    norm_scale <- inter_norm_fact_opts$norm_scale
    norm_fact_inter <- switch(
      norm_scale,
      design = {
        inter_norm_fact[col_ids]
      },
      batch = {
        assert_that(!is.null(in_batch),
          msg = "`in_batch` cannot be null when norm scale is 'batch'"
        )
        inter_norm_fact[[in_batch]][col_ids]
      },
      group = {
        assert_that(!is.null(in_batch),
          msg = "`in_batch` cannot be null when norm scale is 'group'"
        )
        assert_that(!is.null(in_group),
          msg = "`in_group` cannot be null when norm scale is 'group'"
        )
        inter_norm_fact[[in_batch]][[in_group]][col_ids]
      },
      {
        logging::logerror("Invalid `norm_scale`: %s",
                          inter_norm_fact_opts$norm_scale)
        stop(1)
      }
    )
    if (rescale_inter_norm && !is.null(norm_fact_inter)) {
      norm_fact_inter <- norm_fact_inter /
        exp(mean(log(norm_fact_inter)))
    }
  }

  if (!is.null(norm_fact_inter)) {
    norm_fact <- sweep(norm_fact, 2,
                       norm_fact_inter[colnames(norm_fact)], FUN = "*")
  }
  norm_fact
}
