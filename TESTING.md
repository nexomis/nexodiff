# Testing Guide for nexodiff

This document covers testing for the `nexodiff` R package.

## Test Structure

```
tests/
├── testthat.R              # Main test entry point
└── testthat/
    ├── _snaps/             # Visual regression snapshots
    ├── test-comp-04.R      # Pairwise comparison tests
    ├── test-exprdata-*.R   # Expression data tests
    ├── test-integration-*   # Integration tests
    └── test-unit-*          # Unit tests
```

## Running Tests

### Quick Start

```R
devtools::test()
```

### Custom Test Runner

```R
# Run all tests (generates test-out.xml)
Rscript run_tests.R

# Verbose output
Rscript run_tests.R --reporter summary

# Filter tests
Rscript run_tests.R exprdata

# Force reinstall
Rscript run_tests.R --force
```

| Option | Description |
|--------|-------------|
| `-f`, `--force` | Force package reinstallation |
| `--reporter <name>` | Specify reporter (default: `junit`) |
| `<filter>` | Run tests matching filter |

## Visual Regression Testing

Uses [`vdiffr::expect_doppelganger()`](https://vdiffr.r-lib.org/reference/expect_doppelganger.html) to compare plots against snapshots in [`tests/testthat/_snaps/`](tests/testthat/_snaps/).

Update snapshots when plots change intentionally:

```R
vdiffr::manage_files()
```

## Code Coverage

```R
Rscript run_coverage.R
```

Creates an HTML coverage report in the current directory.

## Visualizing Test Results with Allure

[Allure](https://allurest.info/) provides enhanced test visualization with interactive reports, historical trends, and execution logs.

### Installing Allure

```bash
# Download and install nvm:
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.4/install.sh | bash

# in lieu of restarting the shell
\. "$HOME/.nvm/nvm.sh"

# Download and install Node.js:
nvm install 24

npm install -g allure
```

### Generating Reports

```bash
# Generate report from test results
allure generate tests

# Open report in browser
allure open tests/testthat
```

## Test Data

Test datasets are in [`inst/extdata/`](inst/extdata/):

| Dataset | Purpose |
|---------|---------|
| `sim_inputs/test01/` | Basic functionality |
| `sim_inputs/test04/` | Large-scale simulation |
| `test_data1/` | Integration tests |

Generate test data with:

```R
test <- nexodiff::make_test_data("test01", tempdir(TRUE))
```

## Continuous Integration

Tests run automatically on GitHub Actions for pushes to `main`/`master`, pull requests, and releases.

## Resources

- [testthat](https://testthat.r-lib.org/)
- [vdiffr](https://vdiffr.r-lib.org/)
- [covr](https://covr.r-lib.org/)
