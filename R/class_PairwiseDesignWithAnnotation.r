#' @include utils.r

NULL

#' Class representing metadata for expression data: Design and Gene annotation
#'
#' An R6 class to represent metadata for expression data. It combines a pairwise
#' design of experiments with a 2-level structure and transcriptome annotation.
#' Specifically, it provides the following features:
#' - Pairwise design of experiment with a 2-level structure:
#'   - Samples are regrouped into *groups*
#'   - Groups are regrouped into *batches*, where one group is designed as the
#' control group used as a reference for pairwise comparisons
#' - Transcriptome annotation:
#'   - Association of transcript IDs with gene IDs, RNA types, and species
#'   - Association of gene IDs with gene symbols, protein IDs (UniProt), etc.
#'
#' @param in_batch A character vector of batch codes. Select or report only
#' samples in those batches.
#' @param in_group A character vector of group codes. Select or report only
#' samples in those groups.
#' @param in_run A character vector of run IDs. Select or report only samples
#' in those runs.
#'
#' For more information on how to specify the arguments, please refer to the
#' vignette defining file formats.
#'
#' @export
PairwiseDesignWithAnnotation <- R6::R6Class( # nolint
  "PairwiseDesignWithAnnotation",
  public <- list(
    #' Initialize a new `PairwiseDesignWithAnnotation` object.
    #'
    #' This method initializes a new `PairwiseDesignWithAnnotation` object.
    #' The object represents metadata for expression data, including pairwise
    #' design of experiments and transcriptome annotation.
    #'
    #' @param pairwise_design_file Path to the pairwise design file, which must
    #' be in .yml or .csv format.
    #' @param annotation_file Path to the annotation file (mapping txid to gid).
    #' @param idmapping_file Path to the idmapping file, which allows the
    #' mapping of uniprot id to gene id.
    #' @param src_dir Path(s) specified only in case of csv file with the path
    #' for each run to the folder with expression files. In case of several runs
    #' defined in the design file, a named vector is required. Can be specified
    #' for YML also but only 1 value for all runs.
    #'
    #' @return A new `PairwiseDesignWithAnnotation` object.
    #' @export
    initialize = function(pairwise_design_file, annotation_file,
      idmapping_file, src_dir = NULL) {

      logging::logdebug("Parse pairwise design")

      file_ext <- tolower(tools::file_ext(pairwise_design_file))
      if (file_ext %in% c("yml", "yaml")) {
        private$initialize_pairwise_design_yml(pairwise_design_file,
          in_src_dir = src_dir)
        if (! is.null(src_dir)) {
          logging::logwarn("`src_dir` not used for yml file")
        }
      }else if (file_ext %in% "csv") {
        if (is.null(src_dir)) {
          logging::logerror(
            "`src_dir` argument is required for csv design file")
          stop()
        }
        for (i in seq(1, length(src_dir))) {
          src_dir[i] <- normalizePath(src_dir[i])
        }
        private$src_dir <- src_dir
        private$initialize_pairwise_design_csv(pairwise_design_file)
      } else {
        logging::logerror("Unrecognized pairwise design file format")
        stop()
      }

      private$initialize_annotation(annotation_file, idmapping_file)
      # RAW NAMES MUST NOT BE USED BECAUSE SAMPLE MIGHT NOT BE UNIQUE

      private$selected_samples <- private$pairwise_design$sample
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

    #' @description
    #' get batch labels
    #' @return named vector with batch label (batch code as names)
    get_b_labels = function() {
      private$b_labels
    },

    #' @description
    #' get group labels
    #' @return named vector with group label (group code as names)
    get_g_labels = function() {
      private$g_labels
    },

    #' @description
    #' get the paths of files with raw expression data for all samples
    #' @return named vector with paths and sample name as keys
    build_file_paths = function() {
      if (length(private$src_dir) > 1) {
        files <- c()
        for (iter_run_id in names(private$src_dir)) {
          iter_sample_info <- unique(subset(self$get_pairwise_design(),
            run_id == iter_run_id)$sample_base)
          iter_files <- mapply(
            get_counts_file_path,
            rep(private$src_dir[iter_run_id], length(iter_sample_info)),
            iter_sample_info,
            SIMPLIFY = TRUE
          )
          names(iter_files) <- unique(subset(self$get_pairwise_design(),
            run_id == iter_run_id)$sample)
          files <- c(files, iter_files)
        }
      }else {
        iter_sample_info <- unique(self$get_pairwise_design()$sample_base)
        files <- mapply(
          get_counts_file_path,
          rep(private$src_dir, length(iter_sample_info)),
          iter_sample_info,
          SIMPLIFY = TRUE
        )
        names(files) <- unique(self$get_pairwise_design()$sample)
      }
      files
    },

    #' @description
    #' get sample names with simple design-based filtering
    #' @param basename whether to return the basename in place of the complete
    #' name (with run id)
    #' @return sample names in character vector
    extract_sample_names = function(in_batch = NULL, in_group = NULL,
      in_run = NULL, basename = FALSE
    ) {
      if (basename) {
        var <- "sample_base"
      } else {
        var <- "sample"
      }
      unique(as.character(dplyr::filter(
        private$pairwise_design,
          (.data$sample %in% private$selected_samples) &
          (.data$batch %in% in_batch | is.null(in_batch)) &
          (.data$group %in% in_group | is.null(in_group)) &
          (.data$run_id %in% in_run | is.null(in_run))
      )[, var]))
    },

    #' @description
    #' Get a the groups foreach batch
    #' @param include_ctrl whether to include control group
    #' @return list of vectors with batch as key and groups as value
    list_groups_per_batches = function(include_ctrl = FALSE) {
      results <- list()
      data_design <- self$get_pairwise_design()

      if (!include_ctrl) {
        data_design <- data_design[!data_design$ctrl, ]
      }

      for (iter_batch in unique(self$get_pairwise_design()$batch)) {
        results[[iter_batch]] <- unique(
          data_design[
            data_design$batch == iter_batch,
           "group"]
        )
      }
      results
    },


    #' @description
    #' Get a all batchs
    #' @return vector with batch code
    list_batches = function() {
      unique(self$get_pairwise_design()$batch)
    },

    #' @description
    #' get the control group per batch
    #' @return get named vector with batch as key and ctrl group as value
    find_control_group_per_batches = function() {
      data_design <- self$get_pairwise_design()
      data_design <- data_design %>%
        dplyr::filter(.data$ctrl) %>%
        dplyr::select(tidyselect::all_of(c("batch", "group"))) %>%
        dplyr::distinct()
      results <- as.vector(data_design$group)
      names(results) <- as.vector(data_design$batch)
      results
    },

    #' @description
    #' get the pairwise_design table with design-based filtering
    #' @return pairwise_design data.frame
    get_pairwise_design = function(in_batch = NULL,
      in_group = NULL, in_run = NULL) {
      as.data.frame(
        dplyr::filter(
          private$pairwise_design,
            (.data$sample %in% private$selected_samples) &
            (.data$batch %in% in_batch | is.null(in_batch)) &
            (.data$group %in% in_group | is.null(in_group)) &
            (.data$run_id %in% in_run | is.null(in_run))
        )
      )
    },

    #' @description
    #' get the id to be used for paired analysis per sample
    #' only return for samples in batch with a paired design
    #' @return named character vector sample_name -> paired_id
    list_paired_id_per_sample = function() {
      index <- self$get_pairwise_design()$sample %in%
        self$extract_sample_names(in_batch = self$is_paired())
      results <- self$get_pairwise_design()[index, "replicate"]
      names(results) <- self$get_pairwise_design()[index, "sample"]
      results
    },


    #' @description
    #' Select samples based on design. The results will be an intersect with
    #' the previous selection. You can reset the object if it's not desired.
    #' @return vector of the selected samples
    filter_and_set_selected_samples = function(in_batch = NULL, in_group = NULL,
      in_run = NULL) {
      private$selected_samples <-
        intersect(private$selected_samples,
          self$extract_sample_names(in_batch, in_group, in_run)
        )
      private$selected_samples
    },

    #' @description
    #' Set selected samples. The results will be an intersect with
    #' the previous selection. You can reset the object if it's not desired.
    #' @param sample_list vector of sample codes to select
    #' @return vector of the selected samples
    set_selected_samples = function(sample_list) {
      private$selected_samples <-
        intersect(private$selected_samples,
          sample_list
        )
      private$selected_samples
    },

    #' @description
    #' get the list of batch with a paired design (DEPRECATED)
    #' @return character vector
    is_paired = function() {
      private$paired
    },

    #' @description
    #' Reset object, this affects:
    #' * samples selection (as initially)
    #' @return vector of the selected samples
    reset = function() {
      private$selected_samples <- private$pairwise_design$sample
      private$selected_samples
    }

  ),
  private = list(
    # Annotation dict
    annotations = NULL,
    # pairwise_design
    pairwise_design = NULL,
    # batch label dict
    b_labels = NULL,
    # group label dict
    g_labels = NULL,
    # list of batch with a paired design
    paired = NULL,
    # absolute path for paths of run
    src_dir = NULL,
    # sample names selected based on filters
    selected_samples = NULL,

    # initialize pairwise_design based on yml input
    initialize_pairwise_design_yml = function(pairwise_design_file,
      in_src_dir) {
      logging::logdebug("Start yaml pairwise_design parsing")
      sample_details <- data.frame(
        sample = character(),
        batch = character(),
        group = character()
      )
      src_dir <- character()
      yamlsheet <- yaml::read_yaml(pairwise_design_file)

      # Parse optional variables dict
      if ("variables" %in% names(yamlsheet)) {
        additional_vars <- yamlsheet$variables
      }else {
        additional_vars <- NULL
      }

      # parse each run in design
      run_id <- 0
      for (run in yamlsheet$design) {
        run_id <- run_id + 1
        srun_id <- paste("run", run_id, sep = "")
        if ("path" %in% names(run)) {
          if (startsWith(run$path, "/")) {
            tmp_src_dir <- c(run$path)
          }else {
            logging::logerror(
              "paths in pairwise_design file should bo absolutes")
            stop()
          }
        } else if (! is.null(in_src_dir)) {
          tmp_src_dir <- in_src_dir
        } else {
          logging::logerror("paths are required!")
          stop()
        }
        names(tmp_src_dir) <- c(srun_id)
        src_dir <- c(src_dir, tmp_src_dir)
        for (iter_batch in names(run$import)) {
          batch <- run$import[[iter_batch]]
          for (iter_group in names(batch)) {
            group <- batch[[iter_group]]
            if (is.null(additional_vars)) {
              if (is.list(group)) {
                tmp <- data.frame(sample = names(group))
              }else {
                tmp <- data.frame(sample = group)
              }
            }else {
              if (is.list(group)) {
                tmp <- data.frame(group)
                tmp <- as.data.frame(t(tmp))
                names(tmp) <- additional_vars
                tmp$sample <- names(group)
              }else {
                tmp <- data.frame(sample = group)
              }
            }
            tmp$batch <- iter_batch
            tmp$b_label <- yamlsheet$b_labels[[iter_batch]]
            tmp$group <- iter_group
            tmp$g_label <- yamlsheet$g_labels[[iter_group]]
            tmp$run_id <- srun_id
            tmp$ctrl <- (iter_group == yamlsheet$ctrlGroups[[iter_batch]])
            row.names(tmp) <- NULL
            sample_details <- dplyr::bind_rows(sample_details, tmp)
          }
        }
      }

      ## parse label annotation for groups and batch

      if ("b_labels" %in% names(yamlsheet)) {
        result <- as.character(yamlsheet$b_labels)
        names(result) <- names(yamlsheet$b_labels)
        private$b_labels <- result
        sample_details$b_label <- private$b_labels[sample_details$batch]
      }
      if ("g_labels" %in% names(yamlsheet)) {
        result <- as.character(yamlsheet$g_labels)
        names(result) <- names(yamlsheet$g_labels)
        private$g_labels <- result
        sample_details$g_label <- private$g_labels[sample_details$group]
      }

      sample_details$replicate <- sample_details$sample # for retrocompatibility
      # with csv files, to be replaced by NAs

      ## populate the replicate factor according to the specifications

      if ("paired_by" %in% names(yamlsheet)) {
        private$paired <- names(yamlsheet$paired_by)
        for (iter_batch in names(yamlsheet$paired_by)) {
          sample_details[sample_details$batch == iter_batch, "replicate"] <- ""
          for (iter_var in yamlsheet$paired_by[[iter_batch]]) {
            sample_details[sample_details$batch == iter_batch, "replicate"] <-
              paste(
              sample_details[sample_details$batch == iter_batch, "replicate"],
              sample_details[sample_details$batch == iter_batch, iter_var])
          }
        }
      } else {
        private$paired <- c()
      }

      ## transform values as needed types
      sample_details$sample <- as.character(sample_details$sample)
      sample_details$replicate <- as.character(sample_details$replicate)
      sample_details$sample_base <- sample_details$sample
      rununique <- unique(sample_details[, c("sample", "run_id")])
      if (length(rununique$sample) > length(unique(rununique$sample))) {
        sample_details$sample <-
          paste(sample_details$run_id, sample_details$sample, sep = "_")
      }
      logging::logdebug("Ending yaml parsing")
      logging::logdebug(str(sample_details))
      private$src_dir <- src_dir
      private$pairwise_design <- sample_details
    },

    # Initialize the pairwise design from a CSV file.
    #
    # This method initializes the pairwise design from a CSV file.
    # The function takes in the following parameters:
    #
    # param pairwise_design_file Path to the pairwise design file in CSV
    # format.
    #
    # return None
    initialize_pairwise_design_csv = function(pairwise_design_file) {

      # Read in the CSV file using readr's read_csv function
      pairwise_design <- read.delim(pairwise_design_file, sep = ";")

      # Convert some columns to character
      for (str_col in c("batch", "group", "sample")) {
        pairwise_design[, str_col] <- as.character(pairwise_design[, str_col])
      }

      # Define mandatory fields for pairwise design CSV file
      mandatory_fields <- c("batch", "group", "ctrl", "sample")

      if (!all(mandatory_fields %in% names(pairwise_design))) {
        logging::logerror("CSV mandatory fields are not present")
        stop()
      }

      # Add default values for g_label and b_label if they are not present
      if (! "g_label" %in% names(pairwise_design)) {
        pairwise_design$g_label <- pairwise_design$group
      }
      if (! "b_label" %in% names(pairwise_design)) {
        pairwise_design$b_label <- pairwise_design$batch
      }

      # Make sure the sample column is unique
      # Create sample_base column and modify sample column for use in paired
      # analysis
      pairwise_design$sample_base <- pairwise_design$sample
      if (! "run_id" %in% names(pairwise_design)) {
        pairwise_design$run_id <- "run1"
      }
      pairwise_design$sample <-
        paste(pairwise_design$run_id, pairwise_design$sample_base, sep = "_")

      # Check that all batches have one and only one control group
      for (i_batch in unique(pairwise_design$batch)) {
        i_index <- (pairwise_design$batch == i_batch) & pairwise_design$ctrl
        group <- unique(pairwise_design$group[i_index])
        if (length(group) != 1) {
          logging::logerror("There is not on and only one ctrl group per batch")
          stop()
        }
        ij_index <- (pairwise_design$batch == i_batch) &
          (pairwise_design$group == group)
        if (! all(pairwise_design$ctrl[ij_index])) {
          logging::logerror("There are groups with mixed ctrl status")
          stop()
        }
      }

      # Save pairwise_design, b_labels, and g_labels to the object
      private$pairwise_design <- pairwise_design
      data_b_labels <- dplyr::distinct(pairwise_design[, c("batch", "b_label")])
      private$b_labels <- data_b_labels$b_label
      names(private$b_labels) <- data_b_labels$batch
      data_g_labels <- dplyr::distinct(pairwise_design[, c("group", "g_label")])
      private$g_labels <- data_g_labels$g_label
      names(private$g_labels) <- data_g_labels$group
    },

    # initialize annotation
    initialize_annotation = function(annotation_file, idmapping_file) {

      logging::logdebug("Begin annotation parsing")

      # Parse annotation
      annotations <- list(
        "gid" = list(),
        "txid" = list(),
        "tgid" = list()
      )

      if (!is.null(idmapping_file)) {
        annotations[["uniprot"]] <- list()
      }

      print("### DEBUG ###")
      print(annotation_file)

      annotation <- readr::read_delim(annotation_file, delim = " ",
        col_names = c("txid", "gid", "type", "symbol", "tax_id", "tax_name"),
        col_types = vroom::cols("c", "c", "c", "c", "c", "c"))

      # Remove the last digits in transcript ids that are version dependent.
      annotation$txid <- stringr::str_remove(annotation$txid, "-\\d+$")

      dest_vector <- c("type", "tax_id", "tax_name", "gid")

      for (dest in dest_vector) {
          annotations[["txid"]][[dest]] <- annotation[[dest]]
          names(annotations[["txid"]][[dest]]) <- annotation$txid
      }

      annotation_gene <- unique(annotation[, -1])
      dest_vector <- c("symbol", "tax_id", "tax_name")
      for (dest in dest_vector) {
        annotations[["gid"]][[dest]] <- annotation_gene[[dest]]
        names(annotations[["gid"]][[dest]]) <- annotation_gene$gid
      }

      if (! is.null(idmapping_file)) {
        logging::logdebug("run create_id_mappings function")
        idmapping <- create_id_mappings(idmapping_file, annotation_file)
        logging::logdebug("finished create_id_mappings function")
        gidmapping <- idmapping %>%
          group_by(gid) %>%
          summarise(uniprot = paste(unique(uniprot), collapse = ";"),
            protein_names = paste(unique(protein_names), collapse = ";"))

        dest_vector <- c("protein_names", "uniprot")
        for (dest in dest_vector) {
          annotations[["gid"]][[dest]] <- gidmapping[[dest]]
        names(annotations[["gid"]][[dest]]) <- gidmapping$gid
        }

        dest_vector <- c("protein_names", "gid")
        for (dest in dest_vector) {
          annotations[["uniprot"]][[dest]] <- (idmapping[[dest]])
          names(annotations[["uniprot"]][[dest]]) <- idmapping$uniprot
        }

        annotations[["uniprot"]][["symbol"]] <-
          annotations[["gid"]][["symbol"]][annotations[["uniprot"]][["gid"]]]
        names(annotations[["uniprot"]][["symbol"]]) <- idmapping$uniprot

      }


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

      logging::logdebug("End of annotation parsing")
    }
  )
)

# utils function for this class only
# scores min given for unique first and then for dups before merging
create_id_mappings <- function(
  uniprot_mapping_file, annotationfile,
  outfile = NULL, reviewed_min_scores = c(1, 4),
  unreviewed_min_scores = c(2, 5)) {
  # package dplyr required for %>% operator
  suppressWarnings(suppressMessages(require(dplyr, quietly = TRUE)))
  idmapping <- read.delim(uniprot_mapping_file, as.is = TRUE)
  annotation <- unique(read.table(
    annotationfile, header = FALSE, as.is = TRUE
  )[, c(2, 4)])
  names(annotation) <- c("gid", "symbol")
  if (ncol(idmapping) == 6) {
    names(idmapping) <-
      c("uniprot", "status", "names", "symbol", "gid", "score")
  }else if (ncol(idmapping) == 5) {
    names(idmapping) <- c("uniprot", "status", "names", "gid", "score")
  }else {
    stop("Uniprot mapping file not recognized")
  }

  idmapping$score <- as.integer(gsub(" out of 5", "", idmapping$score))
  # NB We do not consider protein spanning different genes that may be
  # fusion protein or readthrough transcripts
  idmapping <- idmapping[stringr::str_count(idmapping$gid, ";") == 1, ]
  uniprotgid <- idmapping[, c("uniprot", "gid", "status", "score")]
  uniprotgid <- subset(uniprotgid, uniprotgid$gid != "")
  status_levels <- c("unreviewed", "reviewed")
  uniprotgid$status <- factor(uniprotgid$status, levels = status_levels)

  # get unique
  filter_uniprotgid <- subset(uniprotgid,
    ((uniprotgid$status == "reviewed") &
    (uniprotgid$score >= reviewed_min_scores[1])) |
    ((uniprotgid$status == "unreviewed") &
    (uniprotgid$score >= unreviewed_min_scores[1]))
  )
  l_uniprotgid <-
    tidyr::separate_rows(filter_uniprotgid, "gid", sep = ";") %>%
    bind_rows() %>%
    filter(.data$gid != "")
  dups <- l_uniprotgid$gid[duplicated(l_uniprotgid$gid)] # dups definition here
  l_unique <- subset(l_uniprotgid, !(l_uniprotgid$gid %in% dups))

  # get dups
  filter_uniprotgid <- subset(uniprotgid,
    ((uniprotgid$status == "reviewed") &
    (uniprotgid$score >= reviewed_min_scores[2])) |
    ((uniprotgid$status == "unreviewed") &
    (uniprotgid$score >= unreviewed_min_scores[2]))
  )
  l_uniprotgid <-
    tidyr::separate_rows(filter_uniprotgid, "gid", sep = ";") %>%
    bind_rows() %>%
    filter(.data$gid != "")
  l_dups <- subset(l_uniprotgid, (l_uniprotgid$gid %in% dups))

  resolved <- l_dups %>%
    group_by(.data$gid) %>%
    filter(.data$score == max(c(0, .data$score)),
      .data$status == status_levels[min(c(2, as.integer(.data$status)))]
    ) %>%
    group_by(.data$gid) %>%
    mutate(uniprot = (function(x) paste0(x, collapse = ";"))(.data$uniprot))

  resolved <- subset(resolved, resolved$uniprot != "")
  gid2uniprot <- c(l_unique$uniprot, resolved$uniprot)
  names(gid2uniprot) <- c(l_unique$gid, resolved$gid)

  l_resolved <- tidyr::separate_rows(resolved, "uniprot", sep = ";")

  idmappings_final <- rbind(as.data.frame(l_unique), as.data.frame(l_resolved))
  idmappings_final$status <- NULL
  idmappings_final$score <- NULL

  uniprot2protname <- idmapping$names
  names(uniprot2protname) <- idmapping$uniprot
  idmappings_final$protein_names <- uniprot2protname[idmappings_final$uniprot]
  gid2symbol <- annotation$symbol
  names(gid2symbol) <- as.character(annotation$gid)
  idmappings_final$symbol <- gid2symbol[idmappings_final$gid]
  if (!is.null(outfile)) {
    write.table(idmappings_final, outfile, sep = "\t", row.names = FALSE)
  }
  return(idmappings_final)
}
