#' @include class_ExprData.r

NULL

#' Class representing an expression dataset at gene-level
#'
#' @description
#' Extension of expression dataset at gene-level.
#'
#' @export

ExprDataGene <- R6::R6Class("ExprDataGene", # nolint
  inherit = ExprData,
  public = list(

    #' @description
    #' Initialize `ExprDataGene` object.
    #'
    #' The `ExprDataGene` object is initiated from a `ExprDatatranscript` object
    #'
    #' Raw count are summarized per genes.
    #'
    #' Lengths are computed using a weighted mean. Weights are raw count divided
    #' by length
    #'
    #' @param expr_data_transcript Expression data at transcript level.
    #' (`ExprDataTranscript` object)
    #' @return A new `ExprDataGene` object.
    initialize = function(expr_data_transcript) {

      private$design <- expr_data_transcript$get_design()
      private$annotation <- expr_data_transcript$get_annotation()
      private$with_fixed_length <- expr_data_transcript$is_with_fixed_length()
      private$intra_norm_fact_method <- "none"
      private$inter_norm_fact_opts <- NULL
      private$inter_norm_fact <- NULL
      private$main_etag <- "tgid"
      private$at_gene_level <- TRUE
      # Summarize raw counts

      txid2tgid <- private$annotation$generate_translate_dict(
        expr_data_transcript$get_main_etag(), private$main_etag
      )

      private$raw <- as.data.frame(expr_data_transcript$get_raw())
      private$raw$gene <-
        txid2tgid[row.names(private$raw)]
      private$raw <- as.data.frame(
        private$raw %>%
          dplyr::group_by(gene) %>%
          dplyr::summarise(dplyr::across(tidyselect::everything(), sum))
      )
      row.names(private$raw) <- private$raw$gene
      private$raw$gene <- NULL
      private$raw <- as.matrix(private$raw)

      # Initialize selected etags and norm_fact
      private$selected_ids <- row.names(private$raw)

      private$intra_norm_fact <- matrix(rep(1, length(private$raw)),
        nrow = nrow(private$raw),
        ncol = ncol(private$raw),
        dimnames = list(
          row.names(private$raw),
          colnames(private$raw)
        )
      )

      # Initialize lengths
      if (expr_data_transcript$is_with_fixed_length()) {
        private$len <- private$intra_norm_fact # all lengths are 1
      } else {
        # Ponderation of lengths with normalized values (raw/length)
        norm <- expr_data_transcript$get_raw() / expr_data_transcript$get_len()
        private$len <- as.data.frame(expr_data_transcript$get_len() * norm)

        # Summarize normalized values for later depondaration
        norm <- as.data.frame(norm)
        norm$gene <- txid2tgid[row.names(norm)]
        g_norm <- as.data.frame(
          norm %>%
            dplyr::group_by(gene) %>%
            dplyr::summarise(dplyr::across(tidyselect::everything(), sum))
        )
        row.names(g_norm) <- g_norm$gene
        g_norm$gene <- NULL
        g_norm <- as.matrix(g_norm)
        # WARNING G_NORM CAN BE ZERO !!! DIVISION BY ZERO

        # Summarize ponderated lengths
        private$len$gene <-
          txid2tgid[row.names(private$len)]
        private$len <- as.data.frame(
          private$len %>%
            dplyr::group_by(gene) %>%
            dplyr::summarise(dplyr::across(tidyselect::everything(), sum))
        )
        row.names(private$len) <- private$len$gene
        private$len$gene <- NULL
        # deponderate
        private$len <- as.matrix(private$len) / g_norm
        # replace length of 0 by 1 to avoid dividing by 0 later
        private$len[private$len == 0] <- 1
        # remove nan values
        private$len[is.nan(private$len)] <- 1
      }
    }
  )
)
