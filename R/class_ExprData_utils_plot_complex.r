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
#' @keywords internal
prepare_complex_plot_data <- function(
  selected_ids, design, in_data, in_batch = NULL, in_group = NULL,
  df_design_filter = NULL, plot_scale = "group",
  include_ctrl_at_group_scale = FALSE, tr_fn = NULL
) {

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

  assert_that(all(c("group", "batch") %in% df_design_filter) || is.null)

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
#' @keywords internal
prepare_tag_selection <- function(
  main_etag, annotation, selected_ids, tags = NULL, tag_type = NULL
) {
  # Prepare tag selection
  if (is.null(tags)) {
    return(selected_ids)
  }
  if (! is.null(tag_type)) {
    tags <- annotation$generate_translate_dict(main_etag, tag_type)[tags]
  }
  intersect(tags, selected_ids)
}

#' Helper function to create PCA plots
#'
#' @param in_data Expression data matrix
#' @param in_title Plot title
#' @param data_design Design data frame
#' @param pca_plot_dims Integer vector specifying which principal components to plot
#' @param mshape Character string specifying variable for point shapes
#' @param mcolor Character string specifying variable for point colors
#' @param prcomp_args Arguments for prcomp
#' @param ggplot_mod Additional ggplot modifications
#' @param color_palette Color palette name
#' @param point_size Numeric value for point size (default: 8)
#' @return ggplot object
#' @keywords internal
create_prcomp_plot <- function(
  in_data, in_title, data_design, pca_plot_dims = c(1, 2),
  mshape = "batch", mcolor = "group", prcomp_args = list(),
  ggplot_mod = NULL, color_palette = "BrBG", point_size = 8
) {
  # Validate pca_plot_dims
  assert_that(
    is.numeric(pca_plot_dims) &
      (length(pca_plot_dims) >= 2) &
      all(pca_plot_dims > 0) &
      all((pca_plot_dims - as.integer(pca_plot_dims)) == 0)
  )
  pca_plot_dims <- as.integer(pca_plot_dims)

  # Check if mapping variables exist and are unique
  if (!mshape %in% names(data_design) ||
        length(unique(data_design[[mshape]])) < 2) {
    mshape <- NULL
  }
  if (!mcolor %in% names(data_design) ||
        length(unique(data_design[[mcolor]])) < 2) {
    mcolor <- NULL
  }

  # Run PCA
  prcomp_args$x <- t(in_data)
  pca_res <- do.call(prcomp, prcomp_args)
  pca_data <- as.data.frame(pca_res$x)
  colnames(pca_data) <- paste0("PC", seq_len(ncol(pca_data)))
  pca_data$sample <- row.names(pca_data)

  plot_data <- dplyr::inner_join(
    pca_data,
    data_design,
    by = dplyr::join_by("sample")
  )

  assert_that(
    length(intersect(plot_data$sample, pca_data$sample)) ==
      length(pca_data$sample)
  )

  # Calculate variance explained
  variance_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

  # Generate all pairwise combinations
  dim_combinations <- t(combn(pca_plot_dims, 2))
  colnames(dim_combinations) <- c("x_dim", "y_dim")

  # Create long format data for all combinations
  long_data <- purrr::map_dfr(
    seq_len(nrow(dim_combinations)), 
    function(i) {
      x_dim <- dim_combinations[i, "x_dim"]
      y_dim <- dim_combinations[i, "y_dim"]
      plot_data$pair <- sprintf(
        "X: PC%d (%.1f%%) Y: PC%d (%.1f%%)",
        x_dim, variance_explained[x_dim],
        y_dim, variance_explained[y_dim]
      )
      plot_data$x <- plot_data[[paste0("PC", x_dim)]]
      plot_data$y <- plot_data[[paste0("PC", y_dim)]]
      plot_data$x_lab <- sprintf("PC%d (%.1f%%)", x_dim,
                                 variance_explained[x_dim])
      plot_data$y_lab <- sprintf("PC%d (%.1f%%)", y_dim,
                                 variance_explained[y_dim])

      plot_data
    }
  )
  # Create base plot
  p <- ggplot2::ggplot(
    data = long_data,
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_point(
      ggplot2::aes_string(shape = mshape, color = mcolor),
      size = point_size
    ) +
    ggplot2::ggtitle(in_title) +
    THEME_NEXOMIS +
    ggplot2::theme(legend.position = "bottom")

  # Add faceting for multiple pairs, or axis labels for a single pair
  if (nrow(dim_combinations) > 1) {
    p <- p + ggplot2::facet_wrap(
      ~pair, scales = "free", ncol = ceiling(sqrt(nrow(dim_combinations)))
    ) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL)
  } else {
    p <- p +
      ggplot2::xlab(long_data$x_lab[1]) +
      ggplot2::ylab(long_data$y_lab[1])
  }

  # Add shape mapping if applicable
  if (!is.null(mshape)) {
    p <- p + ggplot2::scale_shape_manual(values = c(15, 16, 17, 3, 4, 8, 0, 1,
                                                    2))
  }
  if (!is.null(mcolor)) {
    p <- p + ggplot2::scale_color_brewer(palette = color_palette)
  }

  if (!is.null(ggplot_mod)) {
    p <- p + ggplot_mod
  }
  p
}

#' Helper function to create correlation plots
#'
#' @param in_data Expression data matrix
#' @param in_title Plot title
#' @param ggplot_mod Additional ggplot modifications
#' @return ggplot object (ggmatrix)
#' @keywords internal
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
  invisible(GGally::ggmatrix_gtable(g))
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
#' @keywords internal
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
#' @keywords internal
arrange_complex_plots <- function(graphs) {
  p <- cowplot::plot_grid(
    plotlist = graphs,
    nrow = as.integer(sqrt(length(graphs)))
  )
  p
}

