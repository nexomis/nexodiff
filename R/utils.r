#' @include nexodiff-package.R

NULL

#' ggplot2 theme for nexomis
#'
#' Simple theme
#' @export
THEME_NEXOMIS <- ggplot2::theme_classic() + ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "#FAFAFA"),
  axis.text.x = ggplot2::element_text(
    angle = 45, size = 10, vjust = 1, hjust = 1
  ),
  axis.text.y = ggplot2::element_text(
    angle = 0, size = 10, vjust = 1, hjust = 1
  ),
  strip.text = ggplot2::element_text(size = 12, color = "black"),
  legend.title = ggplot2::element_text(size = 12, color = "black"),
  legend.text = ggplot2::element_text(size = 10, color = "black"),
  plot.title = ggtext::element_textbox(hjust = 0.5)
)

#' Color palette for nexomis: Down blue
#' @export
DOWN_COLOR <- "#313695" # nolint
#' Color palette for nexomis: Up red
#' @export
UP_COLOR <- "#A50026" # nolint
#' Color palette for nexomis: Mid yellow-white
#' @export
MID_COLOR <- "#FFFFBF" # nolint
#' Color palette for nexomis: Down blue
#' @export
NA_COLOR <- "#474747" # nolint
#' Color palette gradient for nexomis
#' from UP red to DOWN blue through MID yellow-white
#' @export
UP_TO_DOWN_GRADIENT_COLORS <- c( # nolint
   UP_COLOR, "#D73027", "#F46D43",
  "#FDAE61", "#FEE090", MID_COLOR,
  "#E0F3F8", "#ABD9E9", "#74ADD1",
  "#4575B4", DOWN_COLOR)

#' Generate a simple chord plot
#' @param data a tibble data frame with nested gene list
#' @param superset_colname name of the column with the superset code
#' @param set_colname name of the column with the set code
#' @param data_colname name of the column with the nested lists
#' @param data_list_colname  name of the subcolumn within the nested lists that
#' hold the ids of interest.
#' @param scale_groups whether to scale the sizes of groups on the graph
#' @param isolate_supersets  whether to isolate the sets within the same
#' superset
#' @export
simple_chord <- function(data, superset_colname = "batch",
  set_colname = "group", data_colname = "data", data_list_colname = "gid",
  scale_groups = TRUE, isolate_supersets = TRUE) {

  data$unique_name <-
    paste(data[[superset_colname]], data[[set_colname]], sep = ":")

  data$size <- purrr::map_int(
    data[[data_colname]],
    function(x) {
      length(x[[data_list_colname]])
    }
  )

  mat <- as.matrix(purrr::map_dfc(
    seq_len(nrow(data)),
    function(j) {
      tibble::tibble(
        purrr::map_int(
          seq_len(nrow(data)),
          function(i) {
            if (isolate_supersets &
              data[[superset_colname]][[i]] == data[[superset_colname]][[j]]
            ) {
              return(0)
            }else {
              return(length(intersect(
                data[[data_colname]][[i]][[data_list_colname]],
                data[[data_colname]][[j]][[data_list_colname]]
              )))
            }
          }
        ),
        .name_repair = ~ c(data$unique_name[j])
      )
    }
  ))

  row.names(mat) <- colnames(mat)

  supersets <- data[[superset_colname]]
  names(supersets) <- data$unique_name

  super_sets <- unique(data[[superset_colname]])
  n_col <- length(super_sets)
  color_per_supersets <- RColorBrewer::brewer.pal(
    max(3, n_col), name = "Accent")[seq(1, n_col)]

  names(color_per_supersets) <- super_sets

  color_per_sets <- color_per_supersets[data[[superset_colname]]]
  names(color_per_sets) <- data$unique_name

  unique_name2set <- data[[set_colname]]
  names(unique_name2set) <- data$unique_name

  circlize::circos.clear()
  circlize::chordDiagram(
    mat,
    grid.col = color_per_sets,
    annotationTrack = c("grid"),
    symmetric = TRUE,
    group = supersets,
    scale = scale_groups)
  circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim <- circlize::get.cell.meta.data("xlim")
    ylim <- circlize::get.cell.meta.data("ylim")
    sector_name <- circlize::get.cell.meta.data("sector.index")
    circlize::circos.text(mean(xlim), ylim[1], unique_name2set[sector_name],
      facing = "inside", niceFacing = TRUE, adj = c(1, - 1))

  }, bg.border = NA)

}

#' Modified geometric mean function
#' @description
#' Modified geometric mean function to avoid zero values.
#' In brief, the function adds a delta value to the dataset to avoid zeros.
#' It optimize the delta value to be added to the dataset until the difference
#' between (i) the geometric mean of the non zero values in the dataset
#' and the geometric mean of the same values plus delta is less than epsilon.
#' see https://arxiv.org/abs/1806.06403
#' Thanks to Roberto de la Cruz for the implementation (06/05/2020)
#' @param dataset vector of values
#' @param epsilon epsilon value to optimize the delta value added to the dataset
#' @export
geom_mean_modified <- function(dataset, epsilon) {
  if (sum(dataset) == 0) {
    return(0)
  }

  dataset_nozeros <- dataset[dataset > 0]
  geomean_nozeros <- exp(mean(log(dataset_nozeros)))
  epsilon <- epsilon * geomean_nozeros

  #Simple bisection  method to calculate delta: ( (Eq. I) is increasing as
  #consequence of the Superaddivity of the Geometric Mean)
  deltamin <- 0
  deltamax <- (geomean_nozeros + epsilon)
  while (exp(mean(log(dataset_nozeros + deltamax))) - deltamax < epsilon) {
    #Just for data set with very small standard deviation
    deltamin <- deltamax
    deltamax <- deltamax * 2
  }
  delta <- (deltamin + deltamax) / 2
  aux_exp <- exp(mean(log(dataset_nozeros + delta))) - delta
  #Define aux_exp to not repeat operations
  while ((aux_exp - geomean_nozeros) > epsilon) {
    if ((aux_exp < geomean_nozeros)) {
      deltamin <- delta
    } else {
      deltamax <- delta
    }
    delta <- (deltamin + deltamax) / 2
    aux_exp <- exp(mean(log(dataset_nozeros + delta))) - delta
  }
  exp(mean(log(dataset + delta))) - delta
}

#' Helper function to set mean function based on method name
#'
#' @param method_name Method name for mean calculation
#' @return Function for computing mean
#' @keywords internal
set_mean_function <- function(
  method_name = c(
    "median", "geometric", "nz.geometric", "mod.geometric", "arithmetic"
  )
) {
  method_name <- match.arg(method_name)
  switch(
    method_name,
    median = stats::median,
    geometric = function(x) exp(mean(log(x))),
    nz.geometric = function(x) exp(mean(log(x[x != 0]))),
    mod.geometric = function(x) geom_mean_modified(x, 1e-5),
    arithmetic = mean
  )
}

# Utils function to generate file paths for counts files and check
# their existence
get_counts_file_path <- function(src_dir, sample_name) {
  # Construct the potential file paths
  path_abundance <- file.path(src_dir, sample_name, "abundance.h5")
  path_simple <- file.path(src_dir, paste0(sample_name, ".h5"))

  # Check file existence and decide which path to return
  if (file.exists(path_abundance)) {
    return(path_abundance)
  } else if (file.exists(path_simple)) {
    return(path_simple)
  } else {
    stop("No valid file found for sample: ", sample_name)
  }
}

#' Remove version suffixes from a vector of IDs
#'
#' @description
#' This function takes a character vector of IDs and removes a specified suffix
#' pattern. It also checks if removing the suffix results in new duplicate IDs
#' that were not duplicates before the suffix removal.
#'
#' @param ids A character vector of identifiers.
#' @param suffix_pattern A string representing the regular expression pattern
#'   for the suffix to remove. Defaults to `"\\.\\d+$"` (a hyphen followed by one
#'   or more digits at the end of the string). If `NULL` or `FALSE` or an empty
#'   string, the original IDs are returned unmodified.
#'
#' @return A character vector of IDs with the specified suffixes removed.
#'   Issues a warning if new duplicates are created from previously unique
#'   suffixed IDs.
#' @export
#' @examples
#' ids_to_clean <- c("ID-1", "ID-2.1", "ID-2.2", "ID-3", "ID-4.01", "ID-4.001")
#' remove_id_version_suffix(ids_to_clean, suffix_pattern = "\\.\\d+$")
#' # c("ID-1", "ID-2", "ID-2", "ID-3", "ID-4", "ID-4")
#' # Warning: Removing ID version suffixes created new duplicate IDs...
#'
#' remove_id_version_suffix(ids_to_clean, suffix_pattern = "-\\d+$")
#' # c("ID", "ID-2.1", "ID-2.2", "ID", "ID-4.01", "ID-4.001")
#' # Warning: Removing ID version suffixes created new duplicate IDs...
#'
#' remove_id_version_suffix(c("TX.1", "TX.2", "TX.3"), suffix_pattern = NULL)
#' # c("TX.1", "TX.2", "TX.3")
remove_id_version_suffix <- function(ids, suffix_pattern = "\\.\\d+$") {
  if (is.null(suffix_pattern) || suffix_pattern == FALSE || suffix_pattern == "") {
    return(ids)
  }

  if (!is.character(ids)) {
    logging::logerror("Input 'ids' must be a character vector.")
    stop("Input 'ids' must be a character vector.")
  }
  if (length(ids) == 0) {
    return(character(0))
  }

  original_ids <- ids
  modified_ids <- stringr::str_remove(original_ids, suffix_pattern)

  if (sum(duplicated(original_ids)) != sum(duplicated(modified_ids))) {
    stop("Removing ID version suffixes created new duplicate IDs.")
  }

  return(modified_ids)
}

#' Safely translate IDs using a named vector
#'
#' Translates a vector of IDs using a named vector (dictionary). If an ID
#' is not found in the named vector (resulting in NA), the original ID is
#' kept.
#'
#' @param named_vector A named vector where names are the original IDs and
#'   values are the translated IDs.
#' @param vector_to_translate A character vector of IDs to translate.
#' @return A character vector with translated IDs, preserving original IDs
#'   where translation was not possible.
#' @export
#' @examples
#' translation_map <- c("ID1" = "NewID1", "ID2" = "NewID2")
#' ids_to_translate <- c("ID1", "ID3", "ID2")
#' safe_translate_ids(translation_map, ids_to_translate)
#' # Should return: c("NewID1", "ID3", "NewID2")
safe_translate_ids <- function(named_vector, vector_to_translate) {
  translated_ids <- named_vector[as.character(vector_to_translate)]
  na_indices <- is.na(translated_ids)
  translated_ids[na_indices] <- vector_to_translate[na_indices]
  return(translated_ids)
}
