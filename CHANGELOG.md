# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.2] - 2026-03-25

### Changed

- **Plot ID Refactoring**: Refactored plotting functions to use unique IDs for internal operations while displaying potentially non-unique IDs as labels only. This improves the robustness of plot generation and clustering operations.
   - Modified [`extract_data_for_plot()`](R/class_PairwiseComp.r:490) to preserve original column names and handle ID selection more efficiently
   - Updated [`plot()`](R/class_PairwiseComp.r:602) to use `tag_id_unique` (based on `get_main_etag()`) for plot grouping and operations, while using `tag_id_show` for display labels
   - Enhanced [`plot_heatmap()`](R/class_PairwiseComp.r:1059) to use unique IDs for clustering operations while maintaining user-friendly labels
   - Added translation dictionary functionality to safely map between unique IDs and display labels
   - Improved error handling for invalid `tag_id_show` parameters

## [1.2.1] - 2026-03-25

### Fixed 

- **Heatmap Clustering with Duplicate Display IDs**: Fixed an error in [`PairwiseComp$plot_heatmap()`](R/class_PairwiseComp.r:986) where clustering was performed using `tag_id_show` (display IDs like gene symbols) which could contain duplicates. This caused the error `'list' object cannot be coerced to type 'double'` when computing distances. The fix now uses the actual unique `tag_id` (via [`get_main_etag()`](R/class_ExprData.r:1)) for clustering operations, while still displaying the user-friendly `tag_id_show` labels on the plot axes.

## [1.2.0] - 2026-03-24

### Added

- **Salmon Quantification Support**: Added support for Salmon quantification files (quant.sf) as an alternative to Kallisto (.h5 files).
   - New `quant_source` parameter in [`PairwiseDesign`](R/class_PairwiseDesign.r:1) class (default: "kallisto", alternative: "salmon")
   - Updated [`ExprDataTranscript`](R/class_ExprDataTranscript.r:1) to read Salmon quant.sf format (Name, Length, EffectiveLength, TPM, NumReads columns)
   - [`get_counts_file_path()`](R/utils.r:1) now supports multiple Salmon file patterns: quant.sf, quant.sf.gz, <sample>.sf, <sample>.sf.gz
   - [`make_test_data()`](R/make_test_data.r:1) can generate simulated data in both Kallisto (.h5) and Salmon (quant.sf) formats
   - Removed `format` parameter from [`ExprDataTranscript`](R/class_ExprDataTranscript.r:1) (now determined automatically from `quant_source`)

- **Ensembl GFF3 Support**: Added support for parsing Ensembl-style GFF3 files alongside existing NCBI format.
   - New `gff_source` parameter in [`Annotation`](R/class_Annotation.r:1) class and [`parse_gff_to_annotation()`](R/parse_gff_to_annotation.r:1) function
   - Supports both `"ncbi"` (default) and `"ensembl"` formats
   - Ensembl parsing handles `transcript:` and `gene:` prefixes in ID/Parent attributes
   - Extracts `transcript_id`, `gene_id`, and `Name` attributes for Ensembl format
   - C++ implementation updated in [`src/parse_gff_attributes.cpp`](src/parse_gff_attributes.cpp:1) with new `source_db` parameter
   - Ensembl GFF files require explicit `tax_id` (via named vector in annotation parameter) since they lack taxon info in region entries

- **UniProt ID Mapping Enhancement**: [`fetch_id_mapping()`](R/fetch_id_mapping.r:1) now supports `source_db` parameter to fetch mappings from either NCBI GeneID (`"ncbi"`) or Ensembl (`"ensembl"`) via UniProt API

- **Testing Documentation**: Added comprehensive [TESTING.md](TESTING.md) guide covering:
  - Test structure and organization
  - Running tests with `devtools::test()` and custom test runner
  - Visual regression testing with `vdiffr`
  - Coverage reporting and Allure visualization

- **New Test Fixtures**: Added separate GFF test files for NCBI and Ensembl formats:
  - [`tests/testthat/fixtures/unit/ncbi_dummy.gff`](tests/testthat/fixtures/unit/ncbi_dummy.gff:1)
  - [`tests/testthat/fixtures/unit/ensembl_dummy.gff`](tests/testthat/fixtures/unit/ensembl_dummy.gff:1)

### Changed

- **Documentation**:
  - Updated README.md with improved formatting and links to TESTING.md
  - Added Continuous Integration section documenting GitHub Actions triggers

- **Integration Test Reorganization**: Updated integration test suite to include Ensembl GFF3 support and improved snapshot testing

### Removed

### Fixed

### Security

## [1.1.1] - 2026-01-12

### Added
- Started maintaining a changelog from version 1.1.2 onwards.

[Unreleased]: https://github.com/nexomis/nexodiff/compare/v1.1.1...HEAD
[1.1.1]: https://github.com/nexomis/nexodiff/releases/tag/v1.1.1
