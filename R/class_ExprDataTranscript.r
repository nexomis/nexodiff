#' @include class_ExprData.r

NULL

#' Class representing an expression dataset at transcript-level
#'
#' @description
#' Extension of expression dataset at transcript-level
#' @export

ExprDataTranscript <- R6::R6Class("ExprDataTranscript", # nolint
  inherit = ExprData,
  public = list(
    #' @description
    #' Initialize `ExprDataTranscript` object
    #' @param design A pairwise design object \link{PairwiseDesign}
    #' @param annotation An annotation object \link{Annotation}
    #' @param with_fixed_length A boolean indicating whether the NGS library is
    #' of a fixed length like 3 prime or not
    #' @param log_level Level of logging (see logging package). Default = WARN
    #' @param format Format of files. Default=kallisto
    #' @param suffix_pattern See clean_txid_versions method
    #' Default = txid
    #' @return A new `ExprDataTranscript` object.
    initialize = function(
      design, annotation, with_fixed_length = FALSE, log_level = "WARN",
      format = "kallisto", suffix_pattern = "\\.\\d+$"
    ) {
      logging::basicConfig(log_level)

      private$design <- design$clone()
      private$annotation <- annotation$clone()
      private$at_gene_level <- FALSE
      private$intra_norm_fact_method <- "none"
      private$inter_norm_fact_opts <- NULL
      private$inter_norm_fact <- NULL
      private$with_fixed_length <- with_fixed_length
      private$main_etag <- "txid"
      # todo try to detect if fixed length ...

      # Import data

      logging::loginfo("Import transcript abundance files (.h5)")

      files <- private$design$build_file_paths()

      if (!all(file.exists(files))) {
        logging::logerror("The following files are missing")
        logging::logerror(files[!file.exists(files)])
        stop()
      }

      if (format == "kallisto") {

        private$raw <- NULL
        private$selected_ids <- NULL

        for (sample_name in names(files)) {

          fpath <- files[[sample_name]]
          ids <- as.character(rhdf5::h5read(fpath, "aux/ids"))

          # Check for duplicates in ids
          assert_that(!any(duplicated(ids)))

          if (is.null(private$selected_ids)) {
            private$selected_ids <- ids
            private$raw <- matrix(1, nrow = length(ids),
                                  ncol = length(files))
            rownames(private$raw) <- ids
            colnames(private$raw) <- names(files)
            private$len <- private$raw
            slug_ids <- digest::digest(sort(ids))
          } else {
            assert_that(digest::digest(sort(ids)) == slug_ids)
          }

          # Read counts and verify dimensions before assignment
          counts <- rhdf5::h5read(fpath, "est_counts")
          assert_that(length(counts) == length(ids))

          private$raw[ids, sample_name] <- counts

          if (! with_fixed_length) {
            lengths <- rhdf5::h5read(fpath, "aux/eff_lengths")
            assert_that(length(lengths) == length(ids))
            private$len[ids, sample_name] <- lengths
          }
        }

      } else {
        logging::logerror("wrong value for `format` argument")
        stop()
      }

      private$intra_norm_fact <- matrix(rep(1, length(private$raw)),
        nrow = nrow(private$raw),
        ncol = ncol(private$raw),
        dimnames = list(
          row.names(private$raw),
          colnames(private$raw)
        )
      )

      private$len[is.nan(private$len)] <- 1

      self$clean_id_versions(suffix_pattern = suffix_pattern)
    }
  )
)
