library(testthat)

# =============================================================================
# Integration Tests for Annotation Class with HTTPS GFF Sources
# =============================================================================
# These tests verify end-to-end functionality with live NCBI and Ensembl data.
# Note: R6 class checks, inheritance checks, and basic functionality are
# covered in test-unit-annotation-class.R
# =============================================================================

# Helper function to create snapshot data for headers
make_snapshot_data <- function(df, from_id) {
  knitr::kable(
    df[c(1:10, 100:110, 1000:1010, 3000:3010), ], format = "pipe"
  )
}

test_that("NCBI and Ensembl GFF integration tests", {
  skip_on_cran()
  skip_if_offline(host = "ftp.ncbi.nlm.nih.gov")
  skip_if_offline(host = "ftp.ensembl.org")
  skip_if_offline(host = "www.ebi.ac.uk")

  # =============================================================================
  # Setup: Create both Annotation objects once
  # =============================================================================

  ncbi_gff_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/",
    "GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz"
  )

  ensembl_gff_url <- c(
    "244447" = paste0(
      "https://ftp.ensembl.org/pub/release-115/gff3/cynoglossus_semilaevis/",
      "Cynoglossus_semilaevis.Cse_v1.0.115.gff3.gz"
    )
  )

  # Create NCBI annotation object
  annotation_ncbi <- NULL
  expect_no_error({
    annotation_ncbi <- nexodiff::Annotation$new(
      annotation = ncbi_gff_url,
      gff_source = "ncbi"
    )
  }, message = "Initialization from NCBI HTTPS GFF URL failed.")

  expect_false(
    is.null(annotation_ncbi),
    info = "Annotation object from NCBI HTTPS is NULL."
  )

  # Create Ensembl annotation object
  annotation_ensembl <- NULL
  expect_no_error({
    annotation_ensembl <- nexodiff::Annotation$new(
      annotation = ensembl_gff_url,
      gff_source = "ensembl"
    )
  }, message = "Initialization from Ensembl HTTPS GFF URL failed.")

  expect_false(
    is.null(annotation_ensembl),
    info = "Annotation object from Ensembl HTTPS is NULL."
  )

  # =============================================================================
  # Test: Minimum gid→uniprot mapping count (≥5000)
  # =============================================================================

  test_that("NCBI and Ensembl have minimum gid→uniprot mappings", {
    # NCBI check
    ncbi_tgid2uniprot <- annotation_ncbi$generate_translate_dict("tgid", "uniprot")
    ncbi_mapping_count <- length(na.omit(ncbi_tgid2uniprot))
    expect_gte(
      ncbi_mapping_count,
      5000
    )

    # Ensembl check
    ensembl_tgid2uniprot <- annotation_ensembl$generate_translate_dict("tgid", "uniprot")
    ensembl_mapping_count <- length(na.omit(ensembl_tgid2uniprot))
    expect_gte(
      ensembl_mapping_count,
      5000
    )
  })

  # =============================================================================
  # Test: Snapshot testing for table headers
  # =============================================================================

  test_that("NCBI snapshot: txid headers", {
    df_ids <- annotation_ncbi$export_to_df(from = "txid")
    snapshot_data <- make_snapshot_data(df_ids, "txid")
    expect_snapshot_value(snapshot_data, style = "deparse")
  })

  test_that("NCBI snapshot: gid headers", {
    df_ids <- annotation_ncbi$export_to_df(from = "gid")
    snapshot_data <- make_snapshot_data(df_ids, "gid")
    expect_snapshot_value(snapshot_data, style = "deparse")
  })

  test_that("NCBI snapshot: uniprot headers", {
    df_ids <- annotation_ncbi$export_to_df(from = "uniprot")
    snapshot_data <- make_snapshot_data(df_ids, "uniprot")
    expect_snapshot_value(snapshot_data, style = "deparse")
  })

  test_that("Ensembl snapshot: txid headers", {
    df_ids <- annotation_ensembl$export_to_df(from = "txid")
    snapshot_data <- make_snapshot_data(df_ids, "txid")
    expect_snapshot_value(snapshot_data, style = "deparse")
  })

  test_that("Ensembl snapshot: gid headers", {
    df_ids <- annotation_ensembl$export_to_df(from = "gid")
    snapshot_data <- make_snapshot_data(df_ids, "gid")
    expect_snapshot_value(snapshot_data, style = "deparse")
  })

  test_that("Ensembl snapshot: uniprot headers", {
    df_ids <- annotation_ensembl$export_to_df(from = "uniprot")
    snapshot_data <- make_snapshot_data(df_ids, "uniprot")
    expect_snapshot_value(snapshot_data, style = "deparse")
  })

  # =============================================================================
  # Test: Specific gene mappings (SEO1 as a known example)
  # =============================================================================

  test_that("NCBI and Ensembl have correct specific gene mappings", {
    target_txid <- "NM_001178208"
    target_symbol <- "SEO1"
    target_gid <- "851230"
    target_uniprot <- "P39709"

    # NCBI mappings
    ncbi_tx2sym <- annotation_ncbi$generate_translate_dict("txid", "symbol")
    expect_equal(
      ncbi_tx2sym[[target_txid]], target_symbol,
      info = paste("NCBI: txid to symbol mapping failed for", target_txid)
    )

    ncbi_tx2gid <- annotation_ncbi$generate_translate_dict("txid", "gid")
    expect_equal(
      ncbi_tx2gid[[target_txid]], target_gid,
      info = paste("NCBI: txid to gid mapping failed for", target_txid)
    )

    ncbi_gid2uniprot <- annotation_ncbi$generate_translate_dict("gid", "uniprot")
    expect_equal(
      ncbi_gid2uniprot[[target_gid]], target_uniprot,
      info = paste("NCBI: gid to uniprot mapping failed for", target_gid)
    )
  })

  # =============================================================================
  # Test: Round-trip (write to directory and reload)
  # =============================================================================

  test_that("NCBI round-trip test (write to directory and reload)", {
    temp_dir <- tempfile(pattern = "annot_ncbi_test_")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

    # Write to directory
    write_success <- FALSE
    expect_no_error(
      {
        write_success <- annotation_ncbi$write_to_directory(
          directory_path = temp_dir
        )
      },
      message = "Failed to write NCBI annotation to directory."
    )
    expect_true(write_success)
    expect_true(file.exists(file.path(temp_dir, "txid.csv")))

    # Reload from directory
    reloaded_obj <- NULL
    expect_no_error(
      {
        reloaded_obj <- nexodiff::Annotation$new(annotation_dir = temp_dir)
      },
      message = "Failed to reload NCBI annotation from directory."
    )
    expect_false(is.null(reloaded_obj))

    # Verify ID sets match
    from_ids_original <- sort(annotation_ncbi$get_from_ids())
    from_ids_reloaded <- sort(reloaded_obj$get_from_ids())
    expect_equal(
      from_ids_original, from_ids_reloaded,
      info = "NCBI: Original and reloaded objects have different 'from' ID sets."
    )

    # Verify snapshots match for all ID types
    for (from_id in from_ids_original) {
      df_original <- annotation_ncbi$export_to_df(from = from_id)
      df_reloaded <- reloaded_obj$export_to_df(from = from_id)

      expect_equal(
        names(df_original), names(df_reloaded),
        info = paste("NCBI: Column names differ for", from_id)
      )
      expect_equal(
        nrow(df_original), nrow(df_reloaded),
        info = paste("NCBI: Row counts differ for", from_id)
      )
    }
  })

  test_that("Ensembl round-trip test (write to directory and reload)", {
    temp_dir <- tempfile(pattern = "annot_ensembl_test_")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

    # Write to directory
    write_success <- FALSE
    expect_no_error(
      {
        write_success <- annotation_ensembl$write_to_directory(
          directory_path = temp_dir
        )
      },
      message = "Failed to write Ensembl annotation to directory."
    )
    expect_true(write_success)
    expect_true(file.exists(file.path(temp_dir, "txid.csv")))

    # Reload from directory
    reloaded_obj <- NULL
    expect_no_error(
      {
        reloaded_obj <- nexodiff::Annotation$new(annotation_dir = temp_dir)
      },
      message = "Failed to reload Ensembl annotation from directory."
    )
    expect_false(is.null(reloaded_obj))

    # Verify ID sets match
    from_ids_original <- sort(annotation_ensembl$get_from_ids())
    from_ids_reloaded <- sort(reloaded_obj$get_from_ids())
    expect_equal(
      from_ids_original, from_ids_reloaded,
      info = "Ensembl: Original and reloaded objects have different 'from' ID sets."
    )

    # Verify snapshots match for all ID types
    for (from_id in from_ids_original) {
      df_original <- annotation_ensembl$export_to_df(from = from_id)
      df_reloaded <- reloaded_obj$export_to_df(from = from_id)

      expect_equal(
        names(df_original), names(df_reloaded),
        info = paste("Ensembl: Column names differ for", from_id)
      )
      expect_equal(
        nrow(df_original), nrow(df_reloaded),
        info = paste("Ensembl: Row counts differ for", from_id)
      )
    }
  })
})
