#' Extract id_mapping table from UniProt API
#'
#' @description
#' This function extracts the id_mapping table from the UniProt API based on the specified taxonomic ID.
#' It fetches protein information and then maps UniProt IDs to the target gene ID type using
#' the UniProt ID mapping API.
#'
#' @param tax_id Taxonomic ID to filter the UniProt database.
#' @param source_db The source database of gene ID to fetch. Either "ncbi" (for GeneID) or "ensembl".
#'
#' @return A data frame containing the id_mapping table with columns:
#'   uniprot, status, names, gene_names, gid, score
#'
#' @examples
#' \dontrun{
#' id_mapping <- fetch_id_mapping(tax_id = 9606)
#' id_mapping_ensembl <- fetch_id_mapping(tax_id = 9606, source_db = "ensembl")
#' }
#'
#' @export
fetch_id_mapping <- function(tax_id, source_db = c("ncbi", "ensembl"), txid2gid = NULL) {
  source_db <- match.arg(source_db)
  
  # Map style to UniProt database name
  target_db <- switch(source_db,
    "ncbi" = "xref_geneid",
    "ensembl" = "xref_ensembl",
    stop("Unknown source_db: ", source_db)
  )
  
  # Step 1: Fetch protein info without xref (keeps protein names and other info)
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accessio",
    "n%2Creviewed%2Cprotein_name%2Cgene_names%2C",target_db,"%2Cannotation_score",
    "&format=tsv&query=%28organism_id%3A", tax_id, "%29"
  )

  old_timeout <- getOption('timeout')
  options(timeout = 3600.0)
  on.exit(options(timeout = old_timeout))

  protein_info <- data.table::fread(
    url,
    sep = "\t",
    header = TRUE,
    fill = TRUE,
    stringsAsFactors = FALSE
  )

  if (source_db == "ncbi") {
    data.table::setnames(protein_info, 
                         c("Entry", "Reviewed", "Protein names", "Gene Names", "GeneID", "Annotation"),
                         c("uniprot", "status", "names", "gene_names", "gid", "score"),
                         skip_absent = TRUE)
  
  } else if (source_db == "ensembl") {

    # using apply family we need to transform each Ensembl entry with the following.
    # - drop what is between []
    # - split into a character vector of txids separated by ;
    # - remove the versioning
    # - map gid to txid (ensure removing names)
    # - keep unique and reformat with appending ";" after each entry and using a space separator to get one string.

    protein_info[, Ensembl := stringi::stri_replace_all(Ensembl, "", regex = "\\[[^;]*?\\]")]
    
    protein_info[, Ensembl := vapply(
      Ensembl,
      function(ensembl_str) {
        if (is.na(ensembl_str) || ensembl_str == "") {
          return(NA_character_)
        }
                
        # Step 2: Split by semicolon to get transcript IDs
        tx_ids <- unlist(strsplit(ensembl_str, ";", fixed = TRUE))
        tx_ids <- trimws(tx_ids)
        tx_ids <- tx_ids[tx_ids != ""]
        
        if (length(tx_ids) == 0) {
          return(NA_character_)
        }
        
        # Step 3: Remove versioning (e.g., ENST00000637218.2 -> ENST00000637218)
        tx_ids_clean <- sub("\\.[0-9]+$", "", tx_ids)
        
        # Step 4: Map transcript IDs to gene IDs using txid2gid dictionary
        if (is.null(txid2gid)) {
          # If no mapping provided, return NA
          return(NA_character_)
        }
        
        # Look up gene IDs and remove names (unname)
        gene_ids <- unname(txid2gid[tx_ids_clean])
        
        # Remove NA values
        gene_ids <- gene_ids[!is.na(gene_ids)]
        
        if (length(gene_ids) == 0) {
          return(NA_character_)
        }
        
        # Step 5: Keep unique gene IDs
        gene_ids <- unique(gene_ids)
        
        # Step 6: Reformat with ";" after each entry, space separated
        # Format: "gene1; gene2; gene3;"
        paste0(paste(gene_ids, collapse = "; "), ";")
      },
      character(1)
    )]

    # Rename columns to standard format
    data.table::setnames(protein_info, 
                         c("Entry", "Reviewed", "Protein names", "Gene Names", "Ensembl", "Annotation"),
                         c("uniprot", "status", "names", "gene_names", "gid", "score"),
                         skip_absent = TRUE)    
    
  } else {
    stop("Unknown source_db: ", source_db)
  }
  
  return(protein_info)
}
