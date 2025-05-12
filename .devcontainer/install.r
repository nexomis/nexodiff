args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No arguments provided. Please specify a file name.")
}

input_file <- args[1]

listpackages = read.table(input_file, header = F, stringsAsFactors = FALSE)$V1

rspm::enable()

# To keep track
installed_cran <- c()
installed_bioc <- c()
installed_github <- c()
not_installed <- c()

for (pkg_name in listpackages) {
  # Skip comments and blanks
  if ((!startsWith(pkg_name, "#")) & pkg_name != "") {

    # Check if already installed
    if (length(find.package(pkg_name, quiet = TRUE)) > 0) {
      next
    }

    cat("\n######################################################\n")
    cat("#######          START PACKAGE INSTALL         #######\n")
    cat(pkg_name, "\n")
    cat("######################################################\n\n")

    success <- FALSE

    # Try CRAN first
    tryCatch({
      install.packages(pkg_name, clean = TRUE)
      if (length(find.package(pkg_name, quiet = TRUE)) > 0) {
        installed_cran <- c(installed_cran, pkg_name)
        success <<- TRUE
      }
    }, error = function(e) {})

    # Try Bioconductor if not successful
    if (!success) {
      tryCatch({
        BiocManager::install(pkg_name, clean = TRUE, ask = FALSE, update = FALSE)
        if (length(find.package(pkg_name, quiet = TRUE)) > 0) {
          installed_bioc <- c(installed_bioc, pkg_name)
          success <<- TRUE
        }
      }, error = function(e) {})
    }

    # Try GitHub if not successful and pkg_name looks like github ("user/repo")
    if (!success && grepl(".+/.+", pkg_name)) {
      tryCatch({
        devtools::install_github(pkg_name, clean = TRUE)
        # github packages might have different installed name (usually repo)
        # Try to get the repo part to look for
        repo <- sub(".+/(.+)", "\\1", pkg_name)
        if (length(find.package(repo, quiet = TRUE)) > 0) {
          installed_github <- c(installed_github, pkg_name)
          success <<- TRUE
        }
      }, error = function(e) {})
    }

    # Record not installed
    if (!success) {
      not_installed <- c(not_installed, pkg_name)
    }
  }
}

# Print results
cat("\n=================== SUMMARY =====================\n\n")
cat("Packages installed from CRAN:\n")
print(installed_cran)
cat("\nPackages installed from Bioconductor:\n")
print(installed_bioc)
cat("\nPackages installed from Github:\n")
print(installed_github)
cat("\nPackages NOT successfully installed:\n")
print(not_installed)

# Optionally: install system requirements
rspm::install_sysreqs()