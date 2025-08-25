#' @include utils.r
#' @include class_PairwiseDesign.r
#' @include class_Annotation.r
#' @include class_ExprData.r

NULL

#' Pairwise Comparison Class
#'
#' @description A class representing a nested representation of pairwise
#' comparisons given a pairwise design.


#' @details
#' This class serves as a superclass for all methods of pairwise comparisons
#' considered in this package. It defines the shared public and private methods
#' common to all subclass.
#'
#' Pairwise comparisons are used in experimental designs with test vs control
#' samples, comparing their expression tags to determine deregulation between
#' conditions.
#'
#' @section Attributes:
#' - `expr_data`: An ExprData object.
#' - `results`: A nested data.frame with columns "batch", "group", and a nested
#' "data" column containing:
#'   - baseMean: Numeric vector with the base mean of expression value
#'   - log2FoldChange: Numeric vector with the log2FC
#'   - lfcSE: Numeric vector with the log2FC standard error
#'   - pvalue: Numeric vector with the p-values
#'   - padj: Numeric vector with the adjusted p-values
#'   - status: Factor vector with four levels: undetected, filtered, outlier,
#' analyzed
#'   - tag_id: Character vector with the tag id
#' - `opts`: Important options defined during initialization.
#' @param in_batch Vector of batch codes; select or report only samples in these
#' batches
#' @param in_group Vector of group codes; select or report only samples in these
#' groups
#' @param id Character vector; see `to` argument of method
#' `generate_translate_dict` from class \link{Annotation}
#' @param use_padj Whether to use the adjusted p-value or not
#' @param type Type of deregulation: "deregulated", "upregulated", or
#' "downregulated"
#' @param log2_expr Whether to use log2(x+2) transform for expression data
#' @param lfc_abs_lim Threshold for log fold change
#' @param min_signif Threshold for significance
#' @param tag_id_select Tag ID used for selection. Tag id can be "gid", "tgid",
#' "txid", "symbol" or "uniprot"
#' @param tag_id_show Tag ID used for plots. Tag id can be "gid", "tgid",
#' "txid", "symbol" or "uniprot"
#' @param select_ids Vector of IDs (`tag_id_select`) to select the genes
#' @param facet_scales scales = "fixed" by default
#' @param facet_space space = "fixed" by default
#' @param safe_translate If TRUE, the function will keep original ids
#' if the translation fails
#' @export
PairwiseComp <- R6::R6Class("PairwiseComp", # nolint

  public <- list(

    #' @description
    #' Get results
    #'
    #' @return A nested dataframe.
    get_results = function() {
      private$results
    },


    #' @description Retrieve the results table after filtering for specific
    #' groups/batches and with the option to add new columns to extend the tag
    #' ID.
    #'
    #' @details
    #' This function filters the results table based on specified groups and
    #' batches. Note that the term "filter" in this context does not refer to a
    #' statistical filter.
    #'
    #' @param verbose If TRUE, the column names are replaced with more
    #' meaningful names
    #' @param add_ids Character vector; see `to` argument of method
    #' `generate_translate_dict` from class \link{Annotation}
    #' Available options are: "gid", "symbol", "uniprot", "protein_names",
    #' "type", and "tax_name"
    #' @return A data.frame containing filtered results with optional additional
    #' ID columns
    filter_and_get_results = function(in_batch, in_group, verbose = FALSE,
      add_ids = c(
        "gid", "symbol", "uniprot", "protein_names", "type", "tax_name"
      ),
      safe_translate = FALSE) {
      # TODO error handling if batch/group does not exists
      table <- private$results %>%
        dplyr::filter(
          (.data$batch == in_batch) &
          (.data$group == in_group)
        ) %>%
        tidyr::unnest(cols = c("data")) %>%
        dplyr::ungroup()

      annot <- private$expr_data$get_annotation()

      table <- as.data.frame(table)
      base_id <- private$expr_data$get_main_etag()
      table[, base_id] <- table$tag_id
      table$tag_id <- NULL
      for (id in add_ids){
        if (id != base_id) {
          if (safe_translate) {
            table[, id] <- safe_translate_ids(
              annot$generate_translate_dict(
                base_id, id
              ), table[, base_id]
          )
        } else {
            table[, id] <-
              private$expr_data$get_annotation()$generate_translate_dict(
                base_id, id
              )[table[, base_id]]
        }
      }
      }
      if (verbose) {
        verbose_translate <- c(
          "Expression Mean",
          "Log2(Fold-Change)",
          "Standard Deviation",
          "P value",
          "P value adjusted",
          "Gene status",
          "Gene ID (ENTREZ)",
          "Gene Symbol",
          "Protein ID (Uniprot)",
          "Protein names",
          "RNA type",
          "Taxonomic name",
          "Taxonomic ID (NCBI)"
        )
        names(verbose_translate) <- c(
          "baseMean",
          "log2FoldChange",
          "lfcSE",
          "pvalue",
          "padj",
          "status",
          "gid",
          "symbol",
          "uniprot",
          "protein_names",
          "type",
          "tax_name",
          "tax_id"
        )
        for (i in seq_len(ncol(table))) {
          if (names(table)[i] %in% names(verbose_translate)) {
            names(table)[i] <- verbose_translate[names(table)[i]]
          }
        }
      }
      table
    },

    #' @description
    #' Generate a list of deregulated genes for a specific comparison based
    #' on given criteria.
    #' @param ranking A character string specifying the value to use for
    #' ranking. Can be either "log2FoldChange" (default) or "z" for a z-score
    #' calculated from the p-value and fold change direction.
    #' @param top_x An integer. If not NULL (the default), the list is
    #' restricted to the top `top_x` genes after ranking.
    #' @return A character vector
    generate_a_list = function(in_batch, in_group, id = "symbol",
      use_padj = TRUE, type = "deregulated", lfc_abs_lim = 1,
      min_signif = 0.05, ranking = "log2FoldChange", top_x = NULL
    ) {
      # Initial data retrieval and filtering for analyzed genes
      filtered <- self$filter_and_get_results(
        in_batch, in_group,
        add_ids = c(id)
      ) %>%
        dplyr::filter(.data$status == "analyzed")

      # Significance filtering
      if (use_padj) {
        filtered <- filtered %>%
          dplyr::filter(.data$padj < min_signif)
      } else {
        filtered <- filtered %>%
          dplyr::filter(.data$pvalue < min_signif)
      }

      # LFC and deregulation type filtering
      if (type == "deregulated") {
        filtered <- filtered %>%
          dplyr::filter(abs(.data$log2FoldChange) > lfc_abs_lim)
      } else if (type == "upregulated") {
        filtered <- filtered %>%
          dplyr::filter(.data$log2FoldChange > lfc_abs_lim)
      } else if (type == "downregulated") {
        filtered <- filtered %>%
          dplyr::filter(.data$log2FoldChange < -lfc_abs_lim)
      } else {
        logging::logerror(
          "Unexpected type of deregulation: must be one of 'deregulated',",
          "'upregulated', 'downregulated'"
        )
        stop()
      }

      # Rank and select top genes if top_x is specified
      if (!is.null(top_x)) {
        if (ranking == "z") {
          filtered <- filtered %>%
            dplyr::mutate(
              z = sign(.data$log2FoldChange) * qnorm(1 - .data$pvalue)
            )
        } else if (ranking != "log2FoldChange") {
          logging::logerror(
            "Unexpected ranking method: must be 'log2FoldChange' or 'z'"
          )
          stop()
        }

        if (type == "deregulated") {
          filtered <- filtered %>%
            dplyr::arrange(dplyr::desc(abs(.data[[ranking]])))
        } else if (type == "upregulated") {
          filtered <- filtered %>%
            dplyr::arrange(dplyr::desc(.data[[ranking]]))
        } else { # downregulated
          filtered <- filtered %>%
            dplyr::arrange(.data[[ranking]])
        }

        filtered <- filtered %>%
          dplyr::slice_head(n = as.integer(top_x))
      }

      as.vector(filtered[[id]])

    },

    #' @title Generate Crossed Lists from Nested Results
    #' @description Generate crossed lists of genes from all comparisons in
    #' nested results using multiple filtering criteria.
    #'
    #' @param cross_id Character vector; see `id` argument in
    #' \code{\link{generate_a_list}} method
    #' @param cross_type Vector of deregulation types; see `type` argument in
    #' \code{\link{generate_a_list}} method
    #' Available options are: "deregulated", "upregulated", "downregulated"
    #' (default: c("deregulated", "upregulated", "downregulated"))
    #' @param cross_lfc_abs_lim Vector of log fold change thresholds; see
    #' `lfc_abs_lim` argument in \code{\link{generate_a_list}} method
    #' @param cross_min_signif Vector of significance thresholds; see
    #' `min_signif` argument in \code{\link{generate_a_list}} method
    #' @param top_x An integer. If not NULL (the default), the list is
    #' restricted to the top `top_x` genes after ranking.
    #' @return A nested data.frame containing crossed lists of genes for each
    #' combination of filtering criteria
    cross_args_and_generate_lists = function(
      cross_id = c("symbol", "uniprot"),
      use_padj = TRUE,
      cross_type = c("deregulated", "upregulated", "downregulated"),
      cross_lfc_abs_lim = c(log2(1.5), 1),
      cross_min_signif = c(0.01, 0.05),
      top_x = NULL) {
      design <- private$expr_data$get_design()$get_pairwise_design() %>%
        dplyr::filter(! .data$ctrl) %>%
        dplyr::select(tidyselect::all_of(c("batch", "group"))) %>%
        dplyr::distinct()
      init_l <- list(
        design_id = seq_len(nrow(design)),
        type = cross_type,
        lfc_abs_lim = cross_lfc_abs_lim,
        min_signif = cross_min_signif,
        top_x = top_x
      )
      crossed_df <- purrr::cross_df(init_l)
      crossed_df$group <- design$group[crossed_df$design_id]
      crossed_df$batch <- design$batch[crossed_df$design_id]
      crossed_df$design_id <- NULL
      crossed_df <- crossed_df %>%
        dplyr::relocate(
          tidyselect::all_of(c(
            "batch",
            "group",
            "lfc_abs_lim",
            "min_signif",
            "type"
          ))
        )
      crossed_df$data <- purrr::map(
        seq_len(nrow(crossed_df)),
        function(i) {
          purrr::map_dfc(
            cross_id,
            function(x) {
              data <- tibble::tibble(
                var <- self$generate_a_list(
                  in_batch = crossed_df[[i, "batch"]],
                  in_group = crossed_df[[i, "group"]],
                  id = x,
                  type = crossed_df[[i, "type"]],
                  lfc_abs_lim = crossed_df[[i, "lfc_abs_lim"]],
                  min_signif = crossed_df[[i, "min_signif"]],
                  use_padj = use_padj
              )
)
              names(data) <- x
              data
           }
          )
        }
      )
      crossed_df
    },

    #' @description
    #' Write comparison results to xlsx
    #' @param file_suffix suffix for output file
    #' @param output_folder output folder
    #' @return NULL
    write_to_xlsx = function(file_suffix = "_pairwiseComp.xlsx",
      output_folder = ".", in_batch = NULL) {
      gpb <- private$expr_data$get_design(
        )$list_groups_per_batches()
      if (! dir.exists(output_folder)) {
        dir.create(output_folder)
      }
      for (batch in names(gpb)) {
        if (! is.null(in_batch)) {
          if (! batch %in% in_batch) {
            next
          }
        }
        v_tables <- list()
        for (group in gpb[[batch]]) {
          v_tables[[group]] <- self$filter_and_get_results(
            batch,
            group,
            verbose = TRUE
          )
          class(v_tables[[group]][["P value"]]) <- "scientific"
          class(v_tables[[group]][["P value adjusted"]]) <- "scientific"
          v_tables[[group]][["Gene ID (ENTREZ)"]] <-
            as.integer(v_tables[[group]][["Gene ID (ENTREZ)"]])
        }
        options(openxlsx.numFmt = "#,#0.00")
        openxlsx::write.xlsx(v_tables,
          paste(output_folder, "/", batch, file_suffix, sep = ""),
          asTable = TRUE, tableStyle = "TableStyleLight1"
        )
      }
    },

    #' @description
    #' Get the summary in number of deregulated genes from the comparisons
    #' @param ... passed to `cross_args_and_generate_lists` method
    #' @return tibble with summary stats
    generate_summary = function(...) {

      comp_lists <- self$cross_args_and_generate_lists(...)

      comp_lists$n_deg <- purrr::map_int(
        comp_lists$data,
        function(x) (nrow(x))
      )

      comp_lists[, c("batch", "group", "lfc_abs_lim",
        "min_signif", "type", "n_deg")]
    },

    #' @description
    #' Plot the summary in number of deregulated genes from the comparisons
    #' @param ... passed to `generate_summary` method
    #' @return ggplot2 object with summary plot
    plot_summary = function(...) {

      comp_summary <- self$generate_summary(...)

      ggplot2::ggplot(comp_summary,
        ggplot2::aes(.data$n_deg, .data$type, label = .data$n_deg)) +
        ggplot2::geom_bar(stat = "identity") +
        ggh4x::facet_nested(
          formula("2^lfc_abs_lim + min_signif + batch ~ group"),
          scales = "free") +
        ggplot2::geom_label() +
        nexodiff::THEME_NEXOMIS
    },

    #' @description
    #' Return data for plots
    #' @param max_tags maximum number of tags to show
    #' @param in_batches vector of batches code to keep (default all)
    #' @param select_batches vector of batch to mark for selection
    #' @param hard_select whether to remove unselected data
    #' @return long data frame for plot
    extract_data_for_plot = function(use_padj = TRUE, lfc_abs_lim = 1,
      min_signif = 0.05, tag_id_select = "symbol", tag_id_show = "symbol",
      select_ids = NULL, select_batches = NULL, max_tags = 15,
      in_batches = NULL, hard_select = FALSE, safe_translate = TRUE
    ) {
      if (length(select_ids) > max_tags) {
        logging::logerror(
          "number of tags in `select_ids` must be inferior to `max_tags`"
        )
        stop()
      }

      if (max_tags > 15) {
        logging::logwarn("`max_tags` > 15 is not ideal for visibility")
      }
      base_id <- private$expr_data$get_main_etag()

      known_tag_ids <- private$expr_data$get_annotation()$get_from_ids()
      annot <- private$expr_data$get_annotation()
      assert_that(
        tag_id_show %in% annot$get_to_ids(base_id) &
          tag_id_select %in% annot$get_to_ids(base_id)
      )
      keep_vars <- c("batch", "group", "baseMean", "log2FoldChange", "status",
        "lfcSE", "tag_id"
      )
      if (use_padj) {
        keep_vars <- c(keep_vars, "padj")
        pval_name <- "padj"
      } else {
        keep_vars <- c(keep_vars, "pvalue")
        pval_name <- "pvalue"
      }

      data <- private$results %>%
        tidyr::unnest(cols = c("data")) %>%
        dplyr::ungroup() %>%
        dplyr::select(tidyselect::all_of(keep_vars)) %>%
        dplyr::filter(.data$status != "outlier" & .data$status != "undetected"
        ) %>%
        dplyr::rename_with(function(x) ("pval"), tidyselect::all_of(pval_name)
        ) %>%
        dplyr::mutate(pval = dplyr::if_else(
          is.na(.data$pval),
          1,
          .data$pval
        )) %>%
        dplyr::mutate(status = as.character(.data$status)) %>%
        dplyr::mutate(status = dplyr::if_else(
          .data$pval < min_signif &
            abs(.data$log2FoldChange) > lfc_abs_lim &
            .data$status == "analyzed",
          "deregulated",
          .data$status
        ))
      data <- data %>%
        dplyr::mutate(z = sign(.data$log2FoldChange) * qnorm(1 - .data$pval))

      if (!is.null(in_batches)) {
        data <- data %>%
          dplyr::filter(.data$batch %in% in_batches)
      }

      base_id <- private$expr_data$get_main_etag()
      data[, base_id] <- data$tag_id
      data$tag_id <- NULL

      for (id in c(tag_id_select, tag_id_show)){
        if (id %in% annot$get_to_ids(base_id)){
          if (safe_translate) {
            data[, id] <- safe_translate_ids(
              annot$generate_translate_dict(
                base_id, id
              ), as.vector(data[, base_id])[[1]]
            )
          } else {
            data[, id] <- annot$generate_translate_dict(
              base_id, id
            )[as.vector(data[, base_id])[[1]]]
          }
        } else {
          logging::logwarn(
            paste("Cannot find id", id, "in translate dict")
          )
        }
      }

      data <- data %>%
        dplyr::rename_with(function(x) ("tag_id_select"),
          tidyselect::all_of(tag_id_select)
        )
      if (tag_id_select != tag_id_show) {
        data <- data %>%
          dplyr::rename_with(function(x) ("tag_id_show"),
            tidyselect::all_of(tag_id_show)
          )
      } else {
        data$tag_id_show <- data$tag_id_select
      }

      if (! is.null(select_ids)) {
        data$selected <- as.vector(data[, "tag_id_select"])[[1]] %in% select_ids
        data <- data %>%
          dplyr::mutate(status = dplyr::if_else(
            .data$selected,
            paste(.data$status, "selected", sep = "_"),
            .data$status
          )) %>%
          dplyr::arrange(.data$selected)
        if (hard_select & (length(select_ids) > 0)) {
          data <- dplyr::filter(
            data,
            .data[["selected"]]
          )
        }
      }
      if (! is.null(select_batches)) {
        data$selected_batch <- as.vector(data[, "batch"])[[1]] %in%
          select_batches
        data <- data %>%
          dplyr::mutate(status_batch = dplyr::if_else(
            .data$selected_batch,
            "selected",
            "not_selected"
          ))
      }
      data
    },

    #' @description
    #' Return MA plot(s)
    #' @param ... options passed to plot_de method
    #' @return The plot
    plot_ma = function(...) {
      self$plot_de(plot_type = "ma", ...)
    },

    #' @description
    #' Return VULCANO plot(s)
    #' @param ... options passed to plot_de method
    #' @return The plot
    plot_vulcano = function(...) {
      self$plot_de(plot_type = "vulcano", ...)
    },

    #' @description
    #' Return LFC per group plot(s)
    #' @param ... options passed to plot_de method
    #' @return The plot
    plot_lfc_per_group = function(...) {
      self$plot_de(plot_type = "lfc_per_group", ...)
    },

    #' @description
    #' Return LFC per group with tags as facet plot(s)
    #' @param ... options passed to plot_de method
    #' @param geoms vector with geoms see method `plot_de`
    #' @return The plot
    plot_lfc_per_group_facet_tags = function(geoms = c("bar", "errorbar"), ...
    ) {
      self$plot_de(plot_type = "lfc_per_group_facet_tags", geoms = geoms, ...)
    },

    #' @description
    #' Return plot(s)
    #' @param plot_type plot to draw either:
    #' * "ma" : MA-plot
    #' * "vulcano" : vulcano-plot
    #' * "lfc_per_group" : lfc per group
    #' * "lfc_per_group_facet_tags" : lfc per group with tag facet
    #' @param lfc_limits values to draw hline and or vline for lfc
    #' @param max_nrow maximum of facet row for the nested wrap facetting
    #' @param geoms vector with geoms to add:
    #' * point : compatible with all
    #' * bar : only for "lfc_per_group_facet_tags"
    #' * line : only for "lfc_per_group*"
    #' * errorbar : only for "lfc_per_group_facet_tags"
    #' @param show_selected_ids whether or not to show the selected ids
    #' (not compatible with `lfc_per_group*`)
    #' @param tag_ids_size size when tag ids are plotted in plot area
    #' @param tag_ids_alpha alpha when tag ids are plotted in plot area
    #' @param text_in_box whether or not to use box to plot text
    #' (show_selected_ids)
    #' @param batch_layout Define the type of layout for batch
    #' * "facet_nested_wrap"
    #' * "facet_grid_x"
    #' * "facet_grid_y"
    #' * "color_selected" represent all batch on the same plot with a focus on
    #' a selected batch;only with plot_type "lfc_per_group_facet_tags"
    #' with `color_selected` we must use `select_batches` from
    #' `extract_data_for_plot`
    #' @param ... options passed to extract_data_for_plot method
    #' @return The plot
    plot_de = function(plot_type, lfc_limits = NULL, geoms = "point",
      show_selected_ids = TRUE, text_in_box = FALSE, max_nrow = 5,
      tag_ids_size = 4, tag_ids_alpha = 0.5, log2_expr = TRUE,
      batch_layout = "facet_grid_x", facet_scales = "fixed",
      facet_space = "fixed", ...) {

      if (batch_layout == "color_selected" &
        plot_type != "lfc_per_group_facet_tags") {
        logging::error(paste(
          "`batch_layout=\"color_selected\"` is only compatible with",
          "`plot_type=\"lfc_per_group_facet_tags\"`"
        ))
        stop()
      }

      data <- dplyr::filter(self$extract_data_for_plot(...), .data$baseMean > 0)

      if (
        (show_selected_ids & (! "selected" %in% names(data))) |
        plot_type %in% c("lfc_per_group_facet_tags", "lfc_per_group")
        ) {
        show_selected_ids <- FALSE
      }

      batch2label <- private$expr_data$get_design()$get_b_labels()
      group2label <- private$expr_data$get_design()$get_g_labels()
      facet_nrow <- length(unique(names(batch2label)))
      facet_labeller <- ggplot2::labeller(
        batch = ggplot2::as_labeller(batch2label),
        group = ggplot2::as_labeller(group2label)
      )

      # order batch and groups like labels
      data$batch <- factor(data$batch, levels = names(batch2label))
      data$group <- factor(data$group, levels = names(group2label))

      if (plot_type %in% c("ma", "vulcano")) {
        if (batch_layout == "facet_nested_wrap") {
          facet_formula <- formula("~ batch + group")
        } else if (batch_layout == "facet_grid_y") {
          facet_formula <- formula("batch ~ group")
        } else if (batch_layout == "facet_grid_x") {
          facet_formula <- formula("group ~ batch")
        } else if (batch_layout == "color_selected") {
          facet_formula <- formula("~ group")
        }
      }

      if (plot_type == "ma") {
        if (log2_expr) {
          data$x <- log2(data$baseMean)
          x_lab <- "M - Log2(BaseMean)"
        } else {
          data$x <- data$baseMean
          x_lab <- "Expression BaseMean"
        }
        ltitle <- "Analysis outcome"
        data$aes <- data$status
        data$aes_group <- data$tag_id_select
        data$y <- data$log2FoldChange
        y_lab <- "A - Log2(FoldChange)"
      } else if (plot_type == "vulcano") {
        data$aes <- data$status
        ltitle <- "Analysis outcome"
        data$aes_group <- data$tag_id_select
        data$x <- data$log2FoldChange
        data$y <- - log10(data$pval)
        x_lab <- "Log2(FoldChange)"
        y_lab <- "- Log10(P-value)"
      } else if (
        plot_type %in% c("lfc_per_group_facet_tags", "lfc_per_group")
      ) {
        y_lab <- "Log2(FoldChange)"
        x_lab <- "Groups"
        data$y <- data$log2FoldChange
        data$x <- factor(group2label[data$group], levels = group2label)
        if (plot_type == "lfc_per_group") {
          facet_formula <- formula("~ batch")
        } else if (plot_type == "lfc_per_group_facet_tags") {
          if (! "selected" %in% names(data)) {
            logging::logerror(paste(
              "no tags have been selected; tag selection is required for",
              "plot_type lfc_per_group_facet_tags see the documentation of",
              "PairwiseComp$extract_data_for_plot"
            ))
            stop()
          }
          data <- data[data$selected, ]
          if ("errorbar" %in% geoms) {
            data$ymin <- data$y - data$lfcSE
            data$ymax <- data$y + data$lfcSE
          }

          if (batch_layout == "facet_nested_wrap") {
            facet_formula <- formula("~ batch + tag_id_show")
          } else if (batch_layout == "facet_grid_y") {
            facet_formula <- formula("batch ~ tag_id_show")
          } else if (batch_layout == "facet_grid_x") {
            facet_formula <- formula("tag_id_show ~ batch")
          } else {
            facet_formula <- formula("~ tag_id_show")
          }

          facet_nrow <- length(unique(data$tag_id_select))
        }
      } else {
        logging::logerror("`plot_type` argument is not recognized")
        stop()
      }
      if (plot_type %in% c("ma", "vulcano") & (! identical("point", geoms))) {
        logging::logwarn("For MA and VULCANO plots `geoms` must be point only")
        geoms <- "point"
      }

      # bar and error_bar only compatible with lfc_per_group_facet
      if (
        ("bar" %in% geoms | "error_bar" %in% geoms) &
        plot_type != "lfc_per_group_facet_tags"
      ) {
          logging::logwarn(paste(
            "bar and error_bar `geoms` are only compatible with",
            "lfc_per_group_facet_tags"
          ))
          geoms <- setdiff(geoms, c("bar", "error"))
      }

      if (batch_layout == "color_selected") {
        data <- dplyr::arrange(data, .data$status_batch)
        data$aes <- data$status_batch
        data$aes_group <- data$batch
        ltitle <- "Analysis batch"
      } else {
        data$aes <- data$status
        data$aes_group <- data$tag_id_select
        ltitle <- "Analysis outcome"
      }
      g <-
        ggplot2::ggplot(
          data,
          ggplot2::aes(
            .data$x, .data$y,
            group = .data$aes_group,
            fill = .data$aes,
            color = .data$aes,
            shape = .data$aes,
            alpha = .data$aes,
            size = .data$aes
          )
        )

      if ("bar" %in% geoms) {
        g <- g + ggplot2::geom_bar(
          stat = "identity",
          linewidth = 1,
        )
      }
      if ("line" %in% geoms) {
        g <- g + ggplot2::geom_line(linewidth = 1)
      }
      if ("errorbar" %in% geoms) {
        g <- g + ggplot2::geom_errorbar(
          ggplot2::aes(
            ymin = .data$ymin,
            ymax = .data$ymax
          ),
          color = "darkred",
          linewidth = 1,
          alpha = 1,
          width = 0.5
        )
      }
      if ("point" %in% geoms) {
        g <- g + ggplot2::geom_point()
      }
      if (plot_type %in% c("ma", "vulcano")) {
        g <- g +
          ggplot2::geom_rug(
            mapping = ggplot2::aes(size = NULL),
            alpha = 0.01, linewidth = 0.2, color = "black"
          )
      }
      if (batch_layout == "facet_nested_wrap") {
        g <- g +
          ggh4x::facet_nested_wrap(
            facet_formula,
            nrow = min(max_nrow, facet_nrow),
            labeller = facet_labeller,
            scales = facet_scales
          )
      } else {
        g <- g +
          ggh4x::facet_grid2(
            facet_formula,
            labeller = facet_labeller,
            scales = facet_scales,
            space = facet_space
          )
      }

      if (batch_layout == "color_selected") {
        llabels <- c(
          selected = "selected",
          not_selected = "others"
        )
        g <- g +
          ggplot2::scale_fill_manual(
            values = c(
              not_selected = "skyblue",
              selected = "darkblue"
            ),
            labels = llabels) +
          ggplot2::scale_color_manual(
            values = c(
              not_selected = "skyblue",
              selected = "black"
            ),
            labels = llabels) +
          ggplot2::scale_shape_manual(
            values = c(
              not_selected = 4,
              selected = 4
            ),
            labels = llabels) +
          ggplot2::scale_alpha_manual(
            values = c(
              not_selected = 0.2,
              selected = 0.8
            ),
            labels = llabels) +
          ggplot2::scale_size_manual(
            values = c(
              not_selected = 2,
              selected = 2.5
            ),
            labels = llabels)
      } else {
        llabels <- c(
          filtered = "Unknown (filtered)",
          filtered_selected = "Unknown (filtered) [goi]",
          analyzed = "Not significant",
          deregulated = "Modulated",
          deregulated_selected = "Modulated [goi]",
          analyzed_selected = "Not significant [goi]"
        )
        g <- g +
          ggplot2::scale_fill_manual(
            values = c(
              filtered = "skyblue",
              filtered_selected = "darkblue",
              analyzed = "grey",
              deregulated = "orange",
              deregulated_selected = "red",
              analyzed_selected = "beige"),
            labels = llabels) +
          ggplot2::scale_color_manual(
            values = c(
              filtered = "skyblue",
              filtered_selected = "black",
              analyzed = "grey",
              deregulated = "orange",
              deregulated_selected = "black",
              analyzed_selected = "black"),
            labels = llabels) +
          ggplot2::scale_shape_manual(
            values = c(
              filtered = 4,
              filtered_selected = 4,
              analyzed = 3,
              deregulated = 21,
              deregulated_selected = 21,
              analyzed_selected = 23),
            labels = llabels) +
          ggplot2::scale_alpha_manual(
            values = c(
              filtered = 0.2,
              filtered_selected = 0.8,
              analyzed = 0.2,
              deregulated = 0.2,
              deregulated_selected = 0.8,
              analyzed_selected = 0.8),
            labels = llabels) +
          ggplot2::scale_size_manual(
            values = c(
              filtered = 2,
              filtered_selected = 2.5,
              analyzed = 2,
              deregulated = 2,
              deregulated_selected = 2.5,
              analyzed_selected = 2.5),
            labels = llabels)
      }
      g <- g +
        THEME_NEXOMIS +
        ggplot2::labs(
          x = x_lab,
          y = y_lab,
          fill = ltitle,
          color = ltitle,
          alpha = ltitle,
          size = ltitle,
          shape = ltitle
        )
      if (plot_type %in% c("lfc_per_group_facet_tags", "lfc_per_group")) {
        g <- g +
          ggplot2::scale_y_continuous(
            labels = scales::label_scientific(digits = 2))
      } else {
        g <- g +
          ggplot2::scale_x_continuous(
            labels = scales::label_scientific(digits = 2)) +
          ggplot2::scale_y_continuous(
            labels = scales::label_scientific(digits = 2))
      }
      if (plot_type == "vulcano" & (! is.null(lfc_limits))) {
        g <- g +
          ggplot2::geom_vline(xintercept = lfc_limits, color = "darkred")
      } else if (! is.null(lfc_limits)) {
        g <- g +
          ggplot2::geom_hline(yintercept = lfc_limits, color = "darkred")
      }
      if (show_selected_ids) {
        if (text_in_box) {
          g <- g + ggrepel::geom_label_repel(
            data = data[data$selected, ],
            ggplot2::aes(label = .data$tag_id_show),
            color = "black",
            size = tag_ids_size,
            fill = "antiquewhite",
            alpha = tag_ids_alpha
          )
        } else {
          g <- g + ggrepel::geom_text_repel(
            data = data[data$selected, ],
            ggplot2::aes(label = .data$tag_id_show),
            color = "black",
            size = tag_ids_size,
            alpha = tag_ids_alpha
          )
        }
      }
      g
    },


    #' @description
    #' Return tag hclust
    #' "NA" values are replaced with 0
    #' @param data Optional data similar to extract_data_for_plot output
    #' @param meth_dist method to compute distance, see `dist()` function
    #' @param meth_value value to compute distance
    #' * "z"
    #' * "log2FoldChange"
    #' @param meth_clust method to compute clustering, see `hclust()` function
    #' @param max_abs_value maximum absolute value to consider
    #' @param ... options passed to extract_data_for_plot method
    #' @return The plot
    build_tag_hclust = function(data = NULL, meth_dist = "minkowski",
      meth_value = "z", meth_clust = "centroid", max_abs_value = 100,
      ...
    ) {

      if (is.null(data)) {
        data <- self$extract_data_for_plot(...)
      }

      data <- data %>%
        dplyr::mutate(tag_value = .data[[meth_value]]) %>%
        dplyr::mutate(tag_value = dplyr::if_else(
          .data$tag_value > max_abs_value, max_abs_value, .data$tag_value)) %>%
        dplyr::mutate(tag_value = dplyr::if_else(
          .data$tag_value < - max_abs_value, - max_abs_value, .data$tag_value))

      tag_x_group_wide <- data %>%
        dplyr::select(tidyselect::all_of(
          c("tag_id_show", "batch", "group", "tag_value")
        )) %>%
        tidyr::pivot_wider(
          names_from = tidyselect::all_of(c("batch", "group")),
          values_from = tidyselect::all_of("tag_value"),
        )


      tag_x_group_wide <- as.data.frame(tag_x_group_wide)
      row.names(tag_x_group_wide) <- tag_x_group_wide$tag_id_show
      tag_x_group_wide$tag_id_show <- NULL

      tag_x_group_wide <- as.matrix(tag_x_group_wide)
      tag_x_group_wide[is.na(tag_x_group_wide)] <- 0


      tag_dist <- dist(tag_x_group_wide, method = meth_dist)


      tag_clust <- hclust(tag_dist, method = meth_clust)

      return(tag_clust)
    },

    #' @description
    #' Return HEATmap for a batch
    #' @param ... options passed to extract_data_for_plot method
    #' @param meth_tag_dist `meth_dist` for `build_tag_hclust()`
    #' @param meth_tag_value `meth_value` for `build_tag_hclust()`
    #' @param meth_tag_clust `meth_clust` for `build_tag_hclust()`
    #' @param max_abs_value_hclust `max_abs_value` for `build_tag_hclust()`
    #' @param meth_tag_dist2 `meth_dist` for `build_tag_hclust()` when finding
    #' gene cluster if cut_tree_k is not null
    #' @param meth_tag_value2 `meth_value` for `build_tag_hclust()` when finding
    #' gene cluster if cut_tree_k is not null
    #' @param meth_tag_clust2 `meth_clust` for `build_tag_hclust()` when finding
    #' gene cluster if cut_tree_k is not null
    #' @param max_abs_value_hclust2 `max_abs_value` for `build_tag_hclust()`
    #' when finding gene cluster if cut_tree_k is not null
    #' @param max_abs_value_plot maximum value for values to plot
    #' @param show_selected_ids whether to show or not the selected ids.
    #' @param plot_value value to plot
    #' @param cut_tree_k if provided with cut the heatmap into pieces based
    #' on the hclust
    #' * "log2FoldChange"
    #' * "z" : z-score
    #' @return The plot
    plot_heatmap = function(meth_tag_dist = "minkowski",
      meth_tag_clust = "centroid", meth_tag_value = "z",
      max_abs_value_hclust = 100, plot_value = "z",
      max_abs_value_plot = 5, cut_tree_k = NULL,
      show_selected_ids = FALSE, meth_tag_dist2 = "minkowski",
      meth_tag_value2 = "z", meth_tag_clust2 = "ward.D2",
      max_abs_value_hclust2 = 100, ...
    ) {

      data <- self$extract_data_for_plot(...)

      if ("selected" %in% names(data)) {
        data <- data[data$selected, ]
      }

      tag_clust <- self$build_tag_hclust(
        data = data,
        meth_dist = meth_tag_dist,
        meth_value = meth_tag_value,
        meth_clust = meth_tag_clust,
        max_abs_value = max_abs_value_hclust
      )

      if (! is.null(cut_tree_k)) {
        tag_clust2 <- self$build_tag_hclust(
          data = data,
          meth_dist = meth_tag_dist2,
          meth_value = meth_tag_value2,
          meth_clust = meth_tag_clust2,
          max_abs_value = max_abs_value_hclust2
        )
      }

      tag_levels <- tag_clust$labels[tag_clust$order]

      data$tag_id_show <- factor(data$tag_id_show, levels = tag_levels)

      data <- data %>%
        dplyr::mutate(tag_value = .data[[plot_value]]) %>%
        dplyr::mutate(tag_value = dplyr::if_else(
          .data$tag_value > max_abs_value_plot,
          max_abs_value_plot, .data$tag_value)) %>%
        dplyr::mutate(tag_value = dplyr::if_else(
          .data$tag_value < - max_abs_value_plot,
          - max_abs_value_plot, .data$tag_value))

      batch2label <- private$expr_data$get_design()$get_b_labels()
      group2label <- private$expr_data$get_design()$get_g_labels()

      facet_labeller <- ggplot2::labeller(
        batch = ggplot2::as_labeller(batch2label)
      )

      # order batch and groups like labels
      data$batch <- factor(data$batch, levels = names(batch2label))
      data$group <- factor(group2label[data$group], levels = group2label)

      make_heatmap <- function(data) {
        g <-
          ggplot2::ggplot(
            data,
            ggplot2::aes(.data$group, .data$tag_id_show, fill = .data$tag_value)
          ) +
          ggplot2::scale_fill_gradient2(
            low = DOWN_COLOR,
            high = UP_COLOR,
            na.value = NA_COLOR,
            mid = MID_COLOR,
            midpoint = 0
          ) +
          ggh4x::facet_grid2(
            formula(" . ~ batch"),
            labeller = facet_labeller,
            space = "free",
            scales = "free_x"
          ) +
          ggplot2::geom_tile() +
          ggplot2::labs(fill = plot_value, y = "Expression tag", x = "Group") +
          THEME_NEXOMIS
          if (!("selected" %in% names(data)) & (show_selected_ids)) {
            g <- g + ggplot2::theme(axis.text.y = ggplot2::element_blank())
          }
        g
      }

      results <- list(
        main = make_heatmap(data),
        clust = tag_clust
      )

      if (! is.null(cut_tree_k)) {
        cut_heatmap <- list()
        clusters <- cutree(tag_clust2, k = cut_tree_k)
        for (i in 1:cut_tree_k) {
          cluster_tags <- names(clusters)[clusters == i]
          cut_heatmap[[i]] <- make_heatmap(dplyr::filter(data,
            .data$tag_id_show %in% cluster_tags
          ))
        }
        results[["clust2"]] <- tag_clust2
        results[["cut_heatmap"]] <- cut_heatmap
      }
      results
    }
  ),
  private = list(
    expr_data = NULL,
    results = NULL,
    opts = NULL
  )
)