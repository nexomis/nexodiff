#' @include nexodiff-package.R
#' @include utils.r
NULL

#' Helper function to prepare data for complex plots
#'
#' @param selected_ids Vector of selected expression tag IDs
#' @param design Design object from ExprData
#' @param in_data Expression data matrix
#' @param in_batch Batch filter (optional)
#' @param in_group Group filter (optional)
#' @param df_design_filter Data frame to filter samples
#' @param plot_scale Scale for plotting ("design", "batch", or "group")
#' @param include_ctrl_at_group_scale Whether to include controls at group scale
#' @param tr_fn Data transformation function
#' @return List with prepared data and metadata
prepare_complex_plot_data <- function(
  selected_ids, design, in_data, in_batch = NULL, in_group = NULL,
  df_design_filter = NULL, plot_scale = "group",
  include_ctrl_at_group_scale = FALSE, tr_fn = NULL
) {
  
  # Apply data transformation if provided
  if (!is.null(tr_fn)) {
    tr_fn_df <- function(df) {
      rn <- row.names(df)
      cn <- colnames(df)
      df <- tr_fn(as.matrix(df))
      row.names(df) <- rn
      colnames(df) <- cn
      df
    }
    in_data <- tr_fn_df(in_data)
  }
  
  # Prepare design data
  data_design <- design$get_pairwise_design()[, c("batch", "group")] %>%
    dplyr::distinct()
  
  if (!is.null(in_batch)) {
    data_design <- dplyr::filter(data_design, .data$batch %in% in_batch)
  }
  if (!is.null(in_group)) {
    data_design <- dplyr::filter(data_design, .data$group %in% in_group)
  }
  
  # Apply design filter if provided
  if (!is.null(df_design_filter)) {
    if (plot_scale == "batch" && (!"batch" %in% names(df_design_filter))) {
      logging::logerror(
        "at batch scale, df_design_filters must precise batch for filters")
      stop()
    } else if (plot_scale == "group" &&
               (!all(c("batch", "group") %in% names(df_design_filter)))) {
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
    df_design_filter <- design$get_pairwise_design()
  }
  
  list(
    data = in_data,
    data_design = data_design,
    df_design_filter = df_design_filter,
    selected_ids = selected_ids
  )
}

#' Helper function to prepare tag selection
#'
#' @param main_etag Main expression tag ID
#' @param annotation Annotation object
#' @param selected_ids Currently selected IDs
#' @param tags Vector of tag ids to use for the plot
#' @param tag_type Type of tag to use
#' @return Vector of selected tag IDs
prepare_tag_selection <- function(
  main_etag, annotation, selected_ids, tags = NULL, tag_type = NULL
) {
  # Prepare tag selection
  if (is.null(tag_type)) {
    dict_ids <- selected_ids
    names(dict_ids) <- dict_ids
  } else {
    dict_ids <- annotation$generate_translate_dict(main_etag, tag_type)
  }
  
  if (is.null(tags)) {
    selected_tag_ids <- names(dict_ids)
  } else {
    selected_tag_ids <- names(dict_ids)[dict_ids %in% tags]
  }

  results <- intersect(selected_tag_ids, selected_ids)
  assert_that(length(results) > 0)
  results

}

#' Helper function to create PCA plots
#'
#' @param in_data Expression data matrix
#' @param in_title Plot title
#' @param data_design Design data frame
#' @param prcomp_args Arguments for prcomp
#' @param prcomp_autoplot_args Arguments for autoplot
#' @param ggplot_mod Additional ggplot modifications
#' @param color_palette Color palette name
#' @return ggplot object
create_prcomp_plot <- function(
  in_data, in_title, data_design, prcomp_args = list(),
  prcomp_autoplot_args = list(), ggplot_mod = NULL,
  color_palette = "BrBg"
) {
  prcomp_args$x <- t(in_data)
  pca_res <- do.call(prcomp, prcomp_args)
  prcomp_autoplot_args$object <- pca_res
  prcomp_autoplot_args$data <- data_design
  g <- do.call(ggplot2::autoplot, prcomp_autoplot_args) +
    ggplot2::ggtitle(in_title) +
    ggplot2::scale_fill_brewer(palette = color_palette) +
    ggplot2::scale_color_brewer(palette = color_palette) +
    THEME_NEXOMIS +
    ggplot2::theme(legend.position = "bottom")
  if (!is.null(ggplot_mod)) {
    g <- g + ggplot_mod
  }
  return(g)
}

#' Helper function to create correlation plots
#'
#' @param in_data Expression data matrix
#' @param in_title Plot title
#' @param ggplot_mod Additional ggplot modifications
#' @return ggplot object (ggmatrix)
create_corr_plot <- function(
  in_data, in_title, ggplot_mod = NULL
) {
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
  if (!is.null(ggplot_mod)) {
    g <- g + ggplot_mod
  }
  return(GGally::ggmatrix_gtable(g))
}

#' Helper function to create hierarchical clustering plots
#'
#' @param in_data Expression data matrix
#' @param in_title Plot title
#' @param data_design Design data frame
#' @param dist_method Distance method
#' @param hclust_method Hierarchical clustering method
#' @param dim_reduce Dimensionality reduction parameter
#' @param clust_bar_var Variables for color bars
#' @param height_main Main plot height
#' @param width_main Main plot width
#' @param color_palette Color palette name
#' @param prcomp_args Arguments for prcomp
#' @param ggplot_mod Additional ggplot modifications
#' @return ggplot object or arranged plots
create_hclust_plot <- function(
  in_data, in_title, data_design, dist_method = "euclidean",
  hclust_method = "ward.D2", dim_reduce = NULL, clust_bar_var = c(),
  height_main = 10, width_main = 4, color_palette = "BrBg",
  prcomp_args = list(), ggplot_mod = NULL
) {
  if (is.null(dim_reduce)) {
    d <- philentropy::distance(
      t(in_data), method = dist_method, as.dist.obj = TRUE,
      use.row.names = TRUE
    )
  } else {
    prcomp_args$x <- t(in_data)
    prcomp_args$rank. <- as.integer(dim_reduce)
    pca_res <- do.call(prcomp, prcomp_args)
    d <- philentropy::distance(
      pca_res$x, method = dist_method, as.dist.obj = TRUE,
      use.row.names = TRUE
    )
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
    if (!var %in% names(data_design)) {
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
      ggplot2::ylab(var) +
      ggplot2::geom_tile(color = "black") + ggplot2::theme_minimal() +
      ggplot2::coord_cartesian(
        xlim = c(-1, nrow(data_design) + 1),
        expand = FALSE) +
      ggplot2::scale_fill_brewer(palette = color_palette)
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

#' Helper function to arrange complex plots
#'
#' @param graphs List of ggplot objects
#' @return Arranged plot grid
arrange_complex_plots <- function(graphs) {
  p <- cowplot::plot_grid(
    plotlist = graphs,
    nrow = as.integer(sqrt(length(graphs)))
  )
  return(p)
}

#' Helper function to check plot complex arguments
#'
#' @param plot_scale Scale for plotting ("design", "batch", or "group")
#' @param inter_norm Whether inter normalization is applied
#' @param inter_norm_fact_opts Inter normalization options
#' @param in_batch Batch filter
#' @return NULL (stops execution if checks fail)
plot_complex_check <- function(
  plot_scale, inter_norm, inter_norm_fact_opts, in_batch
) {
  if (plot_scale == "batch" && inter_norm && 
      inter_norm_fact_opts$norm_scale == "group") {
    logging::logerror(paste(
      "plot_scale=batch is not compatible with group-scale inter",
      "normalization")
    )
    stop()
  }
  
  if (plot_scale == "design" && inter_norm && 
      inter_norm_fact_opts$norm_scale != "design") {
    logging::logerror(paste(
      "plot_scale=design is not compatible with group-scale or",
      "batch-scale in normalization")
    )
    stop()
  }
  
  if (plot_scale == "design" && is.null(in_batch)) {
    logging::logerror("`in_batch` cannot be null when plot_scale is design")
    stop()
  }
}
