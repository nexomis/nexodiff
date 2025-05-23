#' Extract id_mapping table from UniProt API
#'
#' @description
#' This function extracts the id_mapping table from the UniProt API based on the specified taxonomic ID.
#'
#' @param tax_id Taxonomic ID to filter the UniProt database.
#'
#' @return A data frame containing the id_mapping table.
#'
#' @examples
#' \dontrun{
#' id_mapping <- fetch_id_mapping(tax_id = 9606)
#' }
#'
#' @export
fetch_id_mapping <- function(tax_id) {
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accessio",
    "n%2Creviewed%2Cprotein_name%2Cgene_names%2Cxref_geneid%2Cannotation_score",
    "&format=tsv&query=(*)+AND+(model_organism%3A", tax_id, ")",  sep = ""
  )
  id_mapping <- data.table::fread(
    url,
    sep = "\t",
    header = TRUE,
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  return(id_mapping)
}
