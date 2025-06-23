
local({
  temp_lib_path <- withr::local_tempdir()
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
})
