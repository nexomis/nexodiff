#' Helper function for ExprData$plot_sum_per_type_per_sample
#'
#' @param data data frame with columns: batch, group, sample, type, total
#' @param batch2label batch labels mapping
#' @param group2label group labels mapping
#' @param same_scale see ExprData$plot_sum_per_type_per_sample
#' @param horizontal_bar see ExprData$plot_sum_per_type_per_sample
#' @param log2_expr see ExprData$plot_sum_per_type_per_sample
#' @param exclude_type see ExprData$plot_sum_per_type_per_sample
#' @return ggplot2 graph
plot_sum_per_type_helper <- function(
  data, batch2label, group2label, same_scale = TRUE,
  horizontal_bar = TRUE, log2_expr = FALSE, exclude_type = c()
) {

  data <- dplyr::filter(data, ! .data$type %in% exclude_type)

  y_label <- "Sum of expression"
  if (log2_expr) {
    y_label <- paste(y_label, "(log2-transformed)")
  }
  if (! same_scale) {
    f_scales <- "free"
  } else if (horizontal_bar) {
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
      labels = scales::label_scientific(digits = 2)
    ) +
    ggplot2::scale_fill_brewer(palette = "BrBG") +
    THEME_NEXOMIS +
    ggplot2::labs(
      fill = "RNA type",
      x = "Samples",
      y = y_label
    )
  if (horizontal_bar) {
    g <- g + ggplot2::coord_flip()
  }
  g
}

#' Helper function for ExprData$plot_dist_per_sample
#'
#' @param data data frame with columns: batch, group, sample, value
#' @param batch2label batch labels mapping
#' @param group2label group labels mapping
#' @param geoms see ExprData$plot_dist_per_sample
#' @param mean_fun see ExprData$plot_dist_per_sample
#' @param log2_expr see ExprData$plot_dist_per_sample
#' @return ggplot2 graph
plot_dist_per_sample_helper <- function(
  data, batch2label, group2label,
  geoms = c("boxplot"), mean_fun = NULL, log2_expr = TRUE
) {

  if (length(setdiff(geoms, c("boxplot", "violin", "histo"))) > 0) {
    logging::logerror("unrecognized value for geoms")
    stop()
  }

  y_label <- "Expression"
  if (log2_expr) {
    y_label <- paste(y_label, "(log2-transformed)")
  }

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
      labels = scales::label_scientific(digits = 2)
    ) +
    THEME_NEXOMIS +
    ggplot2::labs(
      x = "Samples",
      y = y_label
    )
  g
}

#' Helper function for ExprData$plot_complex (make_plot_complex)
#'
#' @param in_data expression data matrix
#' @param in_title plot title
#' @param design_data design data frame
#' @param plot_type see ExprData$plot_complex
#' @param ggplot_mod see ExprData$plot_complex
#' @param prcomp_args see ExprData$plot_complex
#' @param prcomp_autoplot_args see ExprData$plot_complex
#' @param dist_method see ExprData$plot_complex
#' @param hclust_method see ExprData$plot_complex
#' @param dim_reduce see ExprData$plot_complex
#' @param clust_bar_var see ExprData$plot_complex
#' @param height_main see ExprData$plot_complex
#' @param width_main see ExprData$plot_complex
#' @param color_palette see ExprData$plot_complex
#' @return ggplot2 graph or plot grid
make_plot_complex <- function(
  in_data, in_title, design_data, plot_type, ggplot_mod,
  prcomp_args, prcomp_autoplot_args, dist_method, hclust_method, dim_reduce,
  clust_bar_var, height_main, width_main, color_palette
) {

  var_base <- names(design_data)
  data_grouped <- design_data %>%
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
      df
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
    GGally::ggmatrix_gtable(g)
  } else if (plot_type == "prcomp") {
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
    if (! is.null(ggplot_mod)) {
      g <- g + ggplot_mod
    }
    g
  } else if (plot_type == "hclust") {
    if (is.null(dim_reduce)) {
      d <- philentropy::distance(
        t(in_data), method = dist_method, as.dist.obj = TRUE,
        use.row.names = TRUE
      )
    } else {
      prcomp_args$x <- t(in_data)
      prcomp_args$rank. <- as.integer(dim_reduce)
      pca_res <- do.call(prcomp, prcomp_args)
      str(pca_res)
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
      if (! var %in% names(data_design)) {
        logging::logerror("var in clust_bar_var arg is not recognized")
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
          expand = FALSE
        ) +
        ggplot2::scale_fill_brewer(palette = color_palette)
      list_legs[[nth]] <-
        ggpubr::as_ggplot(ggpubr::get_legend(
          var_plot +
            ggplot2::theme(
              plot.margin = ggplot2::margin(0, 0, 0, 0, "cm")
            ) +
            ggplot2::guides(
              fill = ggplot2::guide_legend(title = var, ncol = 1)
            )
        ))
      ws <- c(ws, 1 + length(unique(data_design[[var]])))
      list_plots[[length(list_plots) + 1]] <- var_plot +
        ggplot2::theme(
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(
            angle = 0, vjust = 0.5, hjust = 1
          ),
          axis.ticks = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          legend.position = "none",
          line = ggplot2::element_blank()
        )
    }
    if (length(list_legs) == 0) {
      egg::ggarrange(
        plots = list_plots,
        ncol = 1,
        heights = hs, draw = FALSE
      )
    } else {
      egg::ggarrange(
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
      )
    }
  }
}
