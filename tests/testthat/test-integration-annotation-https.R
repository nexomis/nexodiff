library(testthat)
# Assuming nexodiff is loaded

test_that("Annotation object can be initialized from a live NCBI GFF URL", {
  skip_on_cran()
  skip_if_offline(host = "ftp.ncbi.nlm.nih.gov")
  skip_if_offline(host = "www.ebi.ac.uk")

  gff_url <- paste(
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/",
    "GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz",
    sep = ""
  )
  # gff_url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/" \

  annotation_obj_https <- NULL

  expect_no_error({
    annotation_obj_https <- nexodiff::Annotation$new(annotation = gff_url)
  }, message = "Initialization from HTTPS GFF URL failed.")

  expect_false(
    is.null(annotation_obj_https),
    info = "Annotation object from HTTPS is NULL."
  )
  expect_true(
    R6::is.R6(annotation_obj_https),
    info = "Annotation object is not an R6 object."
  )
  expect_true(
    inherits(annotation_obj_https, "Annotation"),
    info = "HTTPS Annotation object does not inherit from Annotation."
  )

  # Basic check: Ensure some common ID types are available
  expect_true(
    length(annotation_obj_https$get_from_ids()) > 0,
    info = "No 'from' ID types found in HTTPS annotation."
  )
  expect_true(
    "txid" %in% annotation_obj_https$get_from_ids(),
    info = "'txid' not found as a 'from' ID type."
  )

  # Test writing to directory and reloading
  temp_dir <- tempfile(pattern = "annot_https_test_")
  dir.create(temp_dir)

  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  write_success <- FALSE
  expect_no_error(
    {write_success <- annotation_obj_https$write_to_directory(
      directory_path = temp_dir
    )},
    message = "Failed to write HTTPS annotation to directory."
  )
  expect_true(
    write_success,
    info = "write_to_directory did not return TRUE or indicate success."
  )

  expect_true(
    file.exists(file.path(temp_dir, "txid.csv")),
    info = "txid.csv not found after writing to directory."
  )

  target_txid <- "NM_001178208"
  target_symbol <- "SEO1"
  target_gid <- "851230"
  target_uniprot <- "P39709"

  # Check mappings on the original HTTPS loaded object
  tx2sym <- annotation_obj_https$generate_translate_dict("txid", "symbol")
  expect_equal(
    tx2sym[[target_txid]], target_symbol,
    info = paste("HTTPS: txid to symbol mapping failed for", target_txid)
  )

  tx2gid <- annotation_obj_https$generate_translate_dict("txid", "gid")
  expect_equal(
    tx2gid[[target_txid]], target_gid,
    info = paste("HTTPS: txid to gid mapping failed for", target_txid)
  )
  
  # Check uniprot mapping (gid to uniprot, as txid to uniprot might not be direct)
  gid2uniprot <- annotation_obj_https$generate_translate_dict(
    "gid", "uniprot"
  )
  expect_equal(
    gid2uniprot[[target_gid]], target_uniprot,
    info = paste("HTTPS: gid to uniprot mapping failed for", target_gid)
  )

  # Row count and NA checks for original HTTPS object (yeast genome is much smaller)
  expected_counts <- list(
    txid = 5000, symbol = 5000, gid = 5000, uniprot = 3000
  )
  for (id_type in names(expected_counts)) {
    df <- annotation_obj_https$export_to_df(from = id_type)
    expect_false(any(is.na(df[[id_type]])),
                  info = paste("NAs found in primary ID column for", id_type))
    distinct_rows <- nrow(unique(df[, id_type, drop = FALSE]))
    expect_gt(
      distinct_rows,
      expected_counts[[id_type]]
    )
  }


  reloaded_annotation_obj <- NULL
  expect_no_error(
    {reloaded_annotation_obj <- nexodiff::Annotation$new(
      annotation_dir = temp_dir
    )},
    message = "Failed to reload annotation from directory."
  )

  expect_false(
    is.null(reloaded_annotation_obj),
    info = "Reloaded annotation object is NULL."
  )

  expect_true(
    inherits(reloaded_annotation_obj, "Annotation"),
    info = "Reloaded object does not inherit from Annotation."
  )

  from_ids_original <- sort(annotation_obj_https$get_from_ids())
  from_ids_reloaded <- sort(reloaded_annotation_obj$get_from_ids())
  expect_equal(
    from_ids_original, from_ids_reloaded,
    info = "Original and reloaded objects have different 'from' ID sets."
  )

  for (from_id in from_ids_original) {
    df_original <- annotation_obj_https$export_to_df(from = from_id)
    df_reloaded <- reloaded_annotation_obj$export_to_df(from = from_id)

    if (nrow(df_original) > 0 && from_id %in% names(df_original)) {
      df_original <- df_original[order(df_original[[from_id]]), ]
      df_original <- df_original[, sort(names(df_original))]
    }
    if (nrow(df_reloaded) > 0 && from_id %in% names(df_reloaded)) {
      df_reloaded <- df_reloaded[order(df_reloaded[[from_id]]), ]
      df_reloaded <- df_reloaded[, sort(names(df_reloaded))]
    }

    rownames(df_original) <- NULL # Remove rownames for cleaner comparison
    rownames(df_reloaded) <- NULL

    expect_equal(
      df_original, df_reloaded,
      info = paste(
        "Exported df for 'from' ID type '", from_id,
        "' differs between original and reloaded objects."
      )
    )
  }

})
