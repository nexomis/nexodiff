# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for --force or -f flag
force_flags <- c("-f", "--force")
do_force <- any(args %in% force_flags)

if (do_force) {
  message("Force re-installation enabled.")
  # Remove the force flags from the arguments
  args <- args[!args %in% force_flags]
} else {
  do_force <- FALSE
}

temp_lib_path <- "tmp/test_local_lib"
dir.create(temp_lib_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(temp_lib_path, .libPaths()))

devtools::install_local(
  path = ".",
  force = do_force,
  quiet = FALSE,
  build_manual = FALSE,
  build_vignettes = FALSE,
  upgrade = "never"
)

options(testthat.output_file = "test-out.xml")

# Check if a specific test script is provided
if (length(args) > 0) {
  message("Running tests with filter")
  testthat::test_local(
    ".",
    filter = args[1],
    load_package = "installed"
  )
} else {
  message("Running all tests.")
  testthat::test_local(
    ".",
    reporter = "junit",
    load_package = "installed"
  )
}

