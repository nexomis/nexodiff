
local({
  temp_lib_path <- withr::local_tempdir()
  .libPaths(c(temp_lib_path, .libPaths()))

  covr::report(
    x = covr::package_coverage(
      type = c("tests"),
      install_path = temp_lib_path
    ),
    file = paste0(
      format(as.POSIXct(Sys.time(), tz = "UTC"), "%d-%m-%y_%Hh%Mm%S"),
      "-cov-report.html"
    ),
    browse = FALSE
  )

  testthat::test_local(
    ".",
    load_package = "installed"
  )
})
