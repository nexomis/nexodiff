
path <- getwd()
setwd(path)

temp_lib_path <- tempfile(pattern="Rlib")
dir.create(temp_lib_path)

.libPaths(c(temp_lib_path, .libPaths()))

#system("python3 tests/sim_data_functions.py")

message("### Install depedencies ###")

devtools::install_deps(pkgdir = path, upgrade = "never")

install.packages("covr",
  upgrade = "never",
  build_vignettes = FALSE)

devtools::install_local(path = path,
  force = TRUE,
  build = FALSE,
  upgrade = "never",
  keep_outputs = TRUE,
  build_vignettes = FALSE)

message("### Make tests ###")

devtools::test()

pc <- covr::package_coverage(
  path = ".",
  type = c("tests"),
)

l_cov <- as.integer(covr::percent_coverage(pc, by = "line"))
e_cov <- as.integer(covr::percent_coverage(pc, by = "expression"))

write(
  jinjar::render(
    fs::path(".ci/badge.svg.j2"),
    value = paste(as.character(l_cov), "%", sep = ""),
    type = "line cov",
    color = "blue"),
  file = "badge_linecov.svg"
)

write(
  jinjar::render(
    fs::path(".ci/badge.svg.j2"),
    value = paste(as.character(e_cov), "%", sep = ""),
    type = "expr cov",
    color = "blue"),
  file = "badge_exprcov.svg"
)

covr::report(
  x = pc,
  file = "tests.html",
  browse = FALSE
)
