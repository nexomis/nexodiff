#' Simulate data from config
#'
#' Data are generated in a temporary location and paths are returned
#'
#' @param set_id set_id for data simulation
#' @param tmp_dir where to create the dataset
#' @return create all files and return paths
#' @export
make_test_data <- function(set_id, tmp_dir) {

  config_dir_path <- system.file("extdata",
    paste("sim_inputs", set_id, sep = "/"),
    package = "nexodiff"
  )

  csv_config_dir <- file.path(config_dir_path, "config")

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
  annotation$gene <- paste("SYM", annotation$gene, sep = "")
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
    rhdf5::h5write(raw[, sample][[1]], h5_file_name, "est_counts")
    rhdf5::h5write(len[, sample][[1]], h5_file_name, "aux/eff_lengths")
    rhdf5::h5write(annotation[, "txid"][[1]], h5_file_name, "aux/ids")
  }
  results
}
