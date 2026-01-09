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

# Check for --reporter flag
reporter_flag <- "--reporter"
reporter_idx <- which(args == reporter_flag)
if (length(reporter_idx) > 0) {
  # Get the reporter value (next argument after the flag)
  if (reporter_idx[1] + 1 <= length(args)) {
    reporter <- args[reporter_idx[1] + 1]
    # Remove the flag and its value from arguments
    args <- args[-c(reporter_idx, reporter_idx[1] + 1)]
  } else {
    args <- args[-c(reporter_idx)]
    reporter <- "junit"
  }
} else {
  reporter <- "junit"
}

message("Using reporter: ", reporter)

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

if (reporter == "junit") {
  options(testthat.output_file = "test-out.xml")
}

# Check if a specific test script is provided
if (length(args) > 0) {
  message("Running tests with filter")
  testthat::test_local(
    ".",
    filter = args[1],
    reporter = reporter,
    load_package = "installed"
  )
} else {
  message("Running all tests.")
  testthat::test_local(
    ".",
    reporter = reporter,
    load_package = "installed"
  )
}

