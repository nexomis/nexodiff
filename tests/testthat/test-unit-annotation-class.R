library(testthat)

context("Annotation Class Unit Tests")

test_that("Annotation object can be initialized from a local GFF file", {
  gff_file_path <- "fixtures/unit/dummy.gff"

  annotation_obj <- NULL
  expect_no_error({
    annotation_obj <- nexodiff::Annotation$new(annotation = gff_file_path)
  })

  expect_false(
    is.null(annotation_obj),
    info = "Annotation object is NULL after initialization."
  )
  expect_true(
    R6::is.R6(annotation_obj),
    info = "Object is not an R6 object."
  )
  expect_true(
    inherits(annotation_obj, "Annotation"),
    info = "Object does not inherit from Annotation."
  )

  # Test txid to symbol translation
  tx2symbol <- annotation_obj$generate_translate_dict(
    from = "txid",
    to = "symbol"
  )
  expect_equal(
    tx2symbol[["NM_001_1"]], "test_gene1",
    info = "tx2symbol mapping for NM_001_1 failed"
  )
  expect_equal(
    tx2symbol[["NM_002_2"]], "test_gene2",
    info = "tx2symbol mapping for NM_002_2 failed"
  )

  # Test txid to gid translation
  tx2gid <- annotation_obj$generate_translate_dict(from = "txid", to = "gid")
  expect_equal(
    tx2gid[["NM_001_1"]], "1",
    info = "tx2gid mapping for NM_001_1 failed"
  )
  expect_equal(
    tx2gid[["NM_002_2"]], "2",
    info = "tx2gid mapping for NM_002_2 failed"
  )

  # Test export_to_df for "txid"
  annotation_df_txid <- annotation_obj$export_to_df(from = "txid")
  expect_equal(
    nrow(annotation_df_txid), 15,
    info = "export_to_df(from = 'txid') did not return 15 rows."
  )

  # Check gid column properties
  expect_true(
    "gid" %in% names(annotation_df_txid),
    info = "'gid' column missing in exported txid data frame."
  )
  expect_type(
    annotation_df_txid$gid, "character"
  )
  expect_equal(
    length(unique(annotation_df_txid$gid)), 3,
    info = "'gid' column does not have 3 distinct values."
  )
  expect_setequal(
    unique(annotation_df_txid$gid), c("1", "2", "3")
  )

  # Check symbol column properties
  expect_true(
    "symbol" %in% names(annotation_df_txid),
    info = "'symbol' column missing in exported txid data frame."
  )
  expect_type(
    annotation_df_txid$symbol, "character"
  )
  expect_equal(
    length(unique(annotation_df_txid$symbol)), 3,
    info = "'symbol' column does not have 3 distinct values."
  )
  expect_setequal(
    unique(annotation_df_txid$symbol),
    c("test_gene1", "test_gene2", "test_gene3")
  )

  # Test initialization with gff_type_filter
  annotation_obj_filtered <- NULL
  expect_no_error({
    annotation_obj_filtered <- nexodiff::Annotation$new(
      annotation = gff_file_path,
      gff_type_filter = c(".*RNA")
    )
  })
  expect_false(
    is.null(annotation_obj_filtered),
    info = "Filtered Annotation object is NULL."
  )
  annotation_df_filtered <- annotation_obj_filtered$export_to_df(from = "txid")
  expect_equal(
    nrow(annotation_df_filtered), 10,
    info = "Filtered (.*RNA) did not result in 10 rows for txid export."
  )

})
