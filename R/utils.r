#' @include nexodiff-package.R
#' @importFrom magrittr %>%
#' @importFrom data.table .SD
#' @importFrom data.table :=
#' @import dplyr
#' @import ggfortify

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

DOWN_COLOR <- "#313695" # nolint
UP_COLOR <- "#A50026" # nolint
MID_COLOR <- "#FFFFBF" # nolint
NA_COLOR <- "#474747" # nolint

UP_TO_DOWN_GRADIENT_COLORS <- c( # nolint
  UP_COLOR, "#D73027", "#F46D43",
  "#FDAE61", "#FEE090", MID_COLOR,
  "#E0F3F8", "#ABD9E9", "#74ADD1",
  "#4575B4", DOWN_COLOR)

#' Simulate data from config
#'
#' Data are generated in a temporary location and paths are returned
#'
#' @param set_id set_id for data simulation
#' @return create all files and return paths
#' @export
make_test_data_from_xlsx <- function(set_id) {

  config_dir_path <- system.file("extdata",
    paste("sim_inputs", set_id, sep = "/"),
    package = "nexodiff"
  )
  
  csv_config_dir <- file.path(config_dir_path, "config")

  tmp_dir <- tempdir(TRUE)

  design_file_path <- paste(tmp_dir, "design.csv", sep = "/")
  annotation_file_path <- paste(tmp_dir, "annotation.txt", sep = "/")
  id_mapping_file_path <- paste(tmp_dir, "id_mapping.tab", sep = "/")
  results <- list(
    config = csv_config_dir,
    design = design_file_path,
    annotation = annotation_file_path,
    id_mapping = id_mapping_file_path,
    src_dir = tmp_dir
  )
  
  # Read annotation from CSV file
  annotation_csv_path <- file.path(csv_config_dir, "annotation.csv")
  annotation <- readr::read_csv(annotation_csv_path, col_names = FALSE, 
                               show_col_types = FALSE)
  names(annotation) <- c("gene", "tx", "type", "tax_id", "tax_name")
  annotation$txid <- paste(
    annotation$gene,
    annotation$tx,
    sep = "_"
  )
  annotation$gid <- as.integer(as.factor(annotation$gene))
  annotation <- annotation[, c(
    "txid",
    "gid",
    "type",
    "gene",
    "tax_id",
    "tax_name"
  )]
  annotation$tax_name <- stringr::str_replace(annotation$tax_name, " ", "_")
  readr::write_delim(
    annotation,
    quote = "none",
    file = annotation_file_path,
    delim = " ",
    col_names = FALSE
  )
  id_mapping <- dplyr::distinct(annotation[, c("gene", "gid")])
  id_mapping$status <- "reviewed"
  id_mapping$prot_name <- paste("long protein name for", id_mapping$gene)
  id_mapping$entry <- paste("UN", id_mapping$gid, sep = "")
  id_mapping$gid <- paste(id_mapping$gid, ";", sep = "")
  id_mapping$annot <- "5 out of 5"
  id_mapping <- id_mapping[, c(

    "entry",
    "status",
    "prot_name",
    "gene",
    "gid",
    "annot"
  )]

  names(id_mapping) <- c(
    "Entry",
    "Status",
    "Protein names",
    "Gene names",
    "Cross-reference (GeneID)",
    "Annotation"
  )

  readr::write_delim(
    id_mapping,
    quote = "none",
    file = id_mapping_file_path,
    delim = "\t",
    col_names = TRUE
  )

  # Read design from CSV file
  design_csv_path <- file.path(csv_config_dir, "design.csv")
  design <- readr::read_csv(design_csv_path, show_col_types = FALSE)

  readr::write_delim(
    design,
    quote = "none",
    file = design_file_path,
    delim = ";",
    col_names = TRUE
  )

  # Read raw data from CSV file
  raw_csv_path <- file.path(csv_config_dir, "tx_raw.csv")
  raw <- readr::read_csv(raw_csv_path, show_col_types = FALSE)

  # Read length data from CSV file
  len_csv_path <- file.path(csv_config_dir, "tx_len.csv")
  len <- readr::read_csv(len_csv_path, show_col_types = FALSE)

  for (sample in names(raw)){
    h5_file_name <- paste(tmp_dir, "/", sample, ".h5", sep = "")
    invisible(suppressWarnings(file.remove(h5_file_name)))
    rhdf5::h5createFile(h5_file_name)
    rhdf5::h5createGroup(h5_file_name, "aux")
    rhdf5::h5write(raw[, sample], h5_file_name, "est_counts")
    rhdf5::h5write(len[, sample], h5_file_name, "aux/eff_lengths")
    rhdf5::h5write(annotation$txid, h5_file_name, "aux/ids")
  }
  results
}

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

# geometric mean from https://arxiv.org/abs/1806.06403
geom_mean_modified <- function(dataset, epsilon) {
  #Roberto de la Cruz, 06/05/2020
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

# set mean function
set_mean_function <- function(mean_arg) {
  if (mean_arg == "geometric") {
    return(function(x) (exp(mean(log(x)))))
  } else if (mean_arg == "mod.geometric") {
    return(function(x) (geom_mean_modified(x, 1e-05)))
  }else if (mean_arg == "arithmetic") {
    return(function(x) (mean(x)))
  } else if (mean_arg == "median") {
    return(function(x) (median(x)))
  } else if (mean_arg == "nz.geometric") {
    return(function(x) {
      if (sum(x) == 0) {
        res <- 0
      } else {
        res <- exp(mean(log(x)))
      }
      res
    })
  } else {
    logging::error("mean function not recognized")
  }
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
