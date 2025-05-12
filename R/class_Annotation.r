#' @include utils.r

NULL

#' Class representing annotation for expression data.
#'
#' @description
#' An R6 class to represent transcriptome annotation with protein mapping.
#' Specifically, it provides the following features:
#' - Transcriptome annotation:
#'   - Association of transcript IDs with gene IDs, RNA types, and species
#'   - Association of gene IDs with gene symbols, protein IDs (UniProt), etc.
#' @details
#' For more information on how to specify the arguments, please refer to the
#' vignette defining file formats.
#' @export
Annotation <- R6::R6Class( # nolint
  "Annotation",
  public <- list(
    #' @description
    #' Initialize a new `Annotation` object.
    #'
    #' This method initializes a new `Annotation` object.
    #' @param annotation Path to the annotation file in GFF format / or as
    #'  described in the vignette. It can be a vector of files.
    #' @param format_gff If `TRUE`, the annotation file is in GFF format.
    #' @param idmapping Path to the idmapping file, which allows the
    #' mapping of uniprot id to gene id. If not give, the idmapping file will be
    #' downloaded from uniprot using taxon ID from annotation file.
    #' @return A new `Annotation` object.
    initialize = function(
      annotation = NULL,
      format_gff = TRUE,
      idmapping = NULL,
      log_level = "WARN"
    ) {
      logging::basicConfig(log_level)

      if (format_gff) {
        logging::logdebug("Annotation$initialize: Parsing GFF annotation")
        # If the annotation file is in GFF format, parse it using the parse_gff
        # to_annotation function
        annotation <- data.table::rbindlist(
          lapply(annotation, parse_gff_to_annotation)
        )
      } else {
        # If the annotation file is not in GFF format, read it using
        # data.table::fread
        logging::logdebug("Annotation$initialize: Reading non-GFF annotation")
        annotation <- data.table::rbindlist(lapply(annotation, function(x) {
          data.table::fread(
            x, sep = " ",
            col.names = c(
              "txid", "gid", "type", "symbol", "tax_id", "tax_name"
            ),
            colTypes = c("c", "c", "c", "c", "c", "c")
          )
        }))
      }

      logging::logdebug("Annotation$initialize: Calling private$initialize_annotation")
      private$initialize_annotation(annotation, idmapping)
    },

    #' @description
    #' generate a translate dictionary (named vector) from and id to another
    #' based on the annotation
    #' @param from the id type corresponding to the keys
    #' @param to the id type corresponding to the values
    #' @details
    #' Below is precised what is available and id type for value depending
    #' on id type for keys
    #' * Entrez Gene ID `gid` (Or equivalent on Ensembl,UCSC...)
    #'   * Gene symbol `symbol`
    #'   * Protein name `protein_names`
    #'   * Protein ID `uniprot`
    #'   * Taxon ID `tax_id`
    #'   * Taxon name `tax_name`
    #'   * Typed gene ID `tgid`
    #' * RefSeq transcript ID `txid` (Or equivalent on Ensembl,UCSC...)
    #'   * Entrez Gene ID `gid`
    #'   * rna type `type`
    #'   * Taxon ID `tax_id`
    #'   * Taxon name `tax_name`
    #'   * Typed gene ID `tgid`
    #' * Typed gene ID `tgid` i.e. gene id with type concatenated (sep is _)
    #' (one gene can have several type, some coding, some not coding)
    #'   * Gene symbol `symbol`
    #'   * Protein name `protein_names` (warning even if not coding type)
    #'   * Protein ID `uniprot` (warning even if not coding type)
    #'   * rna type `type`
    #'   * Taxon ID `tax_id`
    #'   * Taxon name `tax_name`
    #' * Uniprot protein ID `uniprot`
    #'   * Entrez Gene ID `gid`
    #'   * Protein name `protein_names`
    #'   * Gene symbol `symbol`
    #' @return a named vector that can translate from an id type to another
    generate_translate_dict = function(from, to) {
      if (! from %in% names(private$annotations)) {
        logging::logerror(
          "The `from` argument is not recognized in annotations"
        )
        logging::logerror(str(private$annotations))
        stop()
      }
      if (! to %in% names(private$annotations[[from]])) {
        logging::logerror(
          "The `to` argument is not recognized in annotations"
        )
        logging::logerror(str(private$annotations[[to]]))
        stop()
      }
      private$annotations[[from]][[to]]
    },

    #' Get the list of id type available for translate dicttionary as keys
    #'
    #' @return vecotr with potentioal key "from"
    get_from_ids = function() {
      names(private$annotations)
    },

    #' Get the list of id type available for translate dicttionary as values
    #' 
    #' @param from the key of the translate dictionary
    #' @return vecotr with potentioal key "to"
    get_to_ids = function(from) {
      if (! from %in% names(private$annotations)) {
        logging::logerror(
          "The `from` argument is not recognized in annotations"
        )
      }
      names(private$annotations[[from]])
    }

  ),
  private <- list(
    # Annotation dict
    annotations = NULL,

    # initialize annotation
    initialize_annotation = function(annotation, idmapping) {
      # Parse annotation
      annotations <- list(
        "gid" = list(),
        "txid" = list(),
        "tgid" = list(),
        "uniprot" = list()
      )

      # Remove the last digits in transcript ids that are version dependent.
      annotation[, txid := stringr::str_remove(txid, "-\\d+$")]

      if (is.null(idmapping)) {
        idmapping <- data.table::rbindlist(lapply(
          unique(annotation$tax_id),
          fetch_id_mapping
        ))
      }

      dest_vector <- c("type", "tax_id", "tax_name", "gid")

      for (dest in dest_vector) {
          annotations[["txid"]][[dest]] <- annotation[[dest]]
          names(annotations[["txid"]][[dest]]) <- annotation$txid
      }

      annotation_gene <- unique(annotation[, !c("txid")])
      dest_vector <- c("symbol", "tax_id", "tax_name")
      for (dest in dest_vector) {
        annotations[["gid"]][[dest]] <- annotation_gene[[dest]]
        names(annotations[["gid"]][[dest]]) <- annotation_gene$gid
      }

      logging::logdebug("Annotation$private$initialize_annotation: Run create_id_mappings function")
      # Get unique gene IDs and symbols from the annotation data
      annotation_dt <- unique(annotation[, .(gid, symbol)])
      # Create the id mappings
      id_mappings <- create_uniprot_gene_id_mappings(idmapping, annotation_dt)

      annotations[["uniprot"]][["gid"]] <- id_mappings[["uniprot2gid"]]
      annotations[["uniprot"]][["protein_names"]] <- id_mappings[["uniprot2name"]]
      annotations[["gid"]][["uniprot"]] <- id_mappings[["gid2uniprot"]]

      annotations[["gid"]][["protein_names"]] <-   annotations[["uniprot"]][["protein_names"]][annotations[["gid"]][["uniprot"]]]
      names(annotations[["gid"]][["protein_names"]]) <- names(annotations[["gid"]][["uniprot"]])

      annotations[["uniprot"]][["symbol"]] <- annotations[["gid"]][["symbol"]][annotations[["uniprot"]][["gid"]]]
      names(annotations[["uniprot"]][["symbol"]]) <- names(annotations[["uniprot"]][["gid"]])

      df <- data.frame(
        txid = names(annotations[["txid"]][["gid"]]),
        gid = annotations[["txid"]][["gid"]],
        type = annotations[["txid"]][["type"]]
      )

      df$tgid <- stringr::str_c(df$gid, df$type, sep = "_")

      annotations[["txid"]][["tgid"]] <- df$tgid
      names(annotations[["txid"]][["tgid"]]) <- df$txid

      df_gid_tgid <- unique(df[, c("gid", "tgid")])
      df_type_tgid <- unique(df[, c("type", "tgid")])

      annotations[["tgid"]] <- list()

      annotations[["tgid"]][["txid"]] <- df[["txid"]]
      names(annotations[["tgid"]][["txid"]]) <- df[["tgid"]]

      annotations[["tgid"]][["gid"]] <- df_gid_tgid[["gid"]]
      names(annotations[["tgid"]][["gid"]]) <- df_gid_tgid[["tgid"]]

      annotations[["tgid"]][["type"]] <- df_type_tgid[["type"]]
      names(annotations[["tgid"]][["type"]]) <- df_type_tgid[["tgid"]]

      dest_vector <-
        c("protein_names", "uniprot", "symbol", "tax_id", "tax_name")
      for (dest in dest_vector) {
        annotations[["tgid"]][[dest]] <-
          annotations[["gid"]][[dest]][df_gid_tgid[["gid"]]]
        names(annotations[["tgid"]][[dest]]) <- df_gid_tgid[["tgid"]]
      }

      private$annotations <- annotations

    }
  )
)

#' utils function for this class only
#' scores min given for unique first and then for dups before merging
#' switch to data.table
#' @param idmapping data.table with id mappings from uniprot KB
#' @param annotation data.table with annotation from gff files
#' @param reviewed_min_scores vector of two integers with min scores for
#' reviewed entries, first for unique and then for dups
#' @param unreviewed_min_scores vector of two integers with min scores for 
#' unreviewed entries, first for unique and then for dups
create_uniprot_gene_id_mappings <- function(
  idmapping,
  annotation,
  unique_reviewed_min_score = 1L,
  dups_reviewed_min_score = 4L,
  unique_unreviewed_min_score = 2L,
  dups_unreviewed_min_score = 5L
  ) {

  if (ncol(idmapping) == 6) {
    cnames <- c("uniprot", "status", "names", "symbol", "gid", "score")
  }else if (ncol(idmapping) == 5) {
    cnames <- c("uniprot", "status", "names", "gid", "score")
  }else {
    # If the idmapping data.table has a different number of columns, stop the execution
    stop("Uniprot mapping file not recognized")
  }

  data.table::setnames(idmapping, cnames)
  status_levels <- c("unreviewed", "reviewed")
  idmapping[, status := factor(status, levels = status_levels)]

  if (!is.numeric(idmapping$score)) {
    # Convert the score column to integer for old uniprotKB format
    idmapping[, score := as.integer(gsub(" out of 5", "", score))]
  }

  # NB: We do not consider protein spanning different genes that may be
  # fusion protein or readthrough transcripts
  idmapping <- idmapping[stringr::str_count(idmapping$gid, ";") == 1, ]
  idmapping[, gid := stringr::str_replace(gid, stringr::fixed(";"), "")]

  uniprot2gid <- idmapping$gid
  names(uniprot2gid) <- idmapping$uniprot
  uniprot2name <- idmapping$names
  names(uniprot2name) <- idmapping$uniprot

  is_dups <- duplicated(idmapping$gid)
  unique_idmapping <- idmapping[
    (!is_dups) & (
     ((status == "reviewed") & (score >= unique_reviewed_min_score)) |
      ((status == "unreviewed") & (score >= unique_unreviewed_min_score))
    )
    , ]
  dups_idmapping <- idmapping[
    (is_dups) & (
     ((status == "reviewed") & (score >= dups_reviewed_min_score)) |
      ((status == "unreviewed") & (score >= dups_unreviewed_min_score))
    )
    , ]

  is_still_dups <- duplicated(dups_idmapping$gid)
  unique_idmapping <- data.table::rbindlist(
    list(unique_idmapping, dups_idmapping[ !is_still_dups,])
  )
  dups_idmapping <- dups_idmapping[is_still_dups, ]
  data.table::setorder(dups_idmapping, -status, -score)

  dups_idmapping <- dups_idmapping[, .SD[1], by = gid]
  data.table::setcolorder(dups_idmapping, names(unique_idmapping))
  final_idmapping <- data.table::rbindlist(
    list(unique_idmapping, dups_idmapping)
  )

  final_idmapping[, score:= NULL] # nolint object_usage_linter
  final_idmapping[, status:= NULL] # nolint object_usage_linter

  gid2uniprot <- final_idmapping$uniprot
  names(gid2uniprot) <- final_idmapping$gid

  return(list(
    uniprot2gid = uniprot2gid,
    gid2uniprot = gid2uniprot,
    uniprot2name = uniprot2name
  ))
}
