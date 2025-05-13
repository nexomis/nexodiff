#' @include nexodiff.r
#' 
NULL

#' Parse GFF file to annotation format
#'
#' This function parses a GFF file and extracts gene annotation information.
#' It processes the attributes column (9th column) to extract transcript ID,
#' gene ID, gene symbol, and taxonomic ID.
#' @param gff_file Path to the GFF file or URL
#' @param type_filter Vector of regex patterns to filter the type column (3rd column)
#'
#' @return A data frame with columns: txid, gid, type, symbol, tax_id
#'
#' @export
parse_gff_to_annotation <- function(
  gff_file,
  type_filter = c(".*RNA", ".*transcript")
) {
  
  if (grepl("^(http|ftp|https)://", gff_file)) {
    temp_file <- tempfile(fileext = ".gz")
    download.file(gff_file, temp_file, mode = "wb")
    gff_file <- temp_file
  }
  
  if (tools::file_ext(gff_file) == "gz") {
    tool <- "zcat"
  } else {
    tool <- "cat"
  }
  
  logging::logdebug("parse_gff_to_annotation: Reading GFF data using fread")
  gff_data <- data.table::fread(
    cmd = paste(tool, gff_file, "| grep -v '^#' | cut -f 1,3,9"),
    sep = "\t",
    header = FALSE,
    col.names = c("source", "type", "attributes"),
    fill = TRUE
  )
  
  gff_data_region <- data.table::copy(gff_data)
  gff_data_region <- gff_data_region[type == "region"]
  gff_data_region[, tax_id := stringr::str_replace(stringr::str_extract(
    gff_data_region$attributes, "taxon:(\\d+)"), "taxon:", "")]
  gff_data_region <- gff_data_region[!is.na(tax_id)]
  gff_data_region <- unique(gff_data_region[, .(source, tax_id)])
  source2taxon <- gff_data_region$tax_id
  names(source2taxon) <- gff_data_region$source
  
  if (!is.null(type_filter) && length(type_filter) > 0) {
    unique_types <- unique(gff_data$type)
    type_pattern <- paste(type_filter, collapse = "|")
    filtered_types <- unique_types[grepl(type_pattern, unique_types, perl = TRUE)]
    gff_data <- gff_data[type %in% filtered_types]
  }
  
  parsed_data <- parse_gff_attributes(
    gff_data$attributes, gff_data$type
  )

  data.table::setDT(parsed_data)

  data.table::setnames(parsed_data, 
                      c("txid", "gid", "type", "symbol"))

  logging::logdebug("parse_gff_to_annotation: Fetching scientific names from EBI API")
  query = paste0(paste(unique(source2taxon), collapse = "%20OR%20tax_id%3D"))
  url = paste0(
    "https://www.ebi.ac.uk/ena/portal/api/search?result=taxon&query=tax_id%3D",
    query,
    "&fields=scientific_name,tax_id&format=tsv",
    sep = ""
  )
  logging::logdebug("parse_gff_to_annotation: Fetching Done")

  tax2name_dt <- data.table::fread(
    cmd = paste0("curl -Ss '", url, "'", sep = ""),
    header = TRUE,
    sep = "\t"
    )

  tax2name <- tax2name_dt$scientific_name
  names(tax2name) <- as.character(tax2name_dt$tax_id)

  parsed_data[, tax_id := source2taxon[gff_data$source]]
  parsed_data[, tax_name := tax2name[parsed_data$tax_id]]
  result <- unique(parsed_data)
  return(result)
}
