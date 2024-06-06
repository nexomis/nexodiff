
path <- getwd()
setwd(path)

message("### Install depedencies ###")

devtools::install_deps(pkgdir = path, upgrade = "never")

message("### roxygenize load classic ###")

roxygen2::roxygenize(package.dir = path)

message("### install ###")

devtools::install_local(path = path,
  force = TRUE,
  build = FALSE,
  upgrade = "never",
  keep_outputs = TRUE,
  build_vignettes = FALSE)
