
temp_lib_path <- "_trash/test_lib"
if (!dir.exists(temp_lib_path)) {
  dir.create(temp_lib_path, recursive = TRUE)
}
.libPaths(c(temp_lib_path, .libPaths()))

devtools::install_local(
  path = ".",
  force = TRUE,
  quiet = FALSE,
  build_manual = FALSE,
  build_vignettes = FALSE,
  upgrade = "never",
  lib = temp_lib_path
)

testthat::test_local(
  ".",
  load_package = "installed"
)