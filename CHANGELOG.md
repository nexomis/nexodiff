# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Ensembl GFF3 Support**: Added support for parsing Ensembl-style GFF3 files alongside existing NCBI format.
  - New `gff_style` parameter in [`Annotation`](R/class_Annotation.r:1) class and [`parse_gff_to_annotation()`](R/parse_gff_to_annotation.r:1) function
  - Supports both `"ncbi"` (default) and `"ensembl"` formats
  - Ensembl parsing handles `transcript:` and `gene:` prefixes in ID/Parent attributes
  - Extracts `transcript_id`, `gene_id`, and `Name` attributes for Ensembl format
  - C++ implementation updated in [`src/parse_gff_attributes.cpp`](src/parse_gff_attributes.cpp:1) with new `style` parameter

- **Testing Documentation**: Added comprehensive [TESTING.md](TESTING.md) guide covering:
  - Test structure and organization
  - Running tests with `devtools::test()` and custom test runner
  - Visual regression testing with `vdiffr`
  - Coverage reporting and Allure visualization

- **New Test Fixtures**: Added separate GFF test files for NCBI and Ensembl formats:
  - [`tests/testthat/fixtures/unit/ncbi_dummy.gff`](tests/testthat/fixtures/unit/ncbi_dummy.gff:1)
  - [`tests/testthat/fixtures/unit/ensembl_dummy.gff`](tests/testthat/fixtures/unit/ensembl_dummy.gff:1)

### Changed

- **Devcontainer Configuration**:
  - Changed VS Code extension from `saoudrizwan.claude-dev` to `kilocode.kilo-code`
  - Removed nvm, Node.js v24, and Allure CLI installation from Dockerfile

- **Documentation**:
  - Updated README.md with improved formatting and links to TESTING.md
  - Added Continuous Integration section documenting GitHub Actions triggers

- **Gitignore**: Added `/allure-report` and `/*gff3.gz` patterns


### Removed

### Fixed

### Security

## [1.1.1] - 2026-01-12

### Added
- Started maintaining a changelog from version 1.1.2 onwards.

[Unreleased]: https://github.com/nexomis/nexodiff/compare/v1.1.1...HEAD
[1.1.1]: https://github.com/nexomis/nexodiff/releases/tag/v1.1.1
