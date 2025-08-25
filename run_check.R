# Use devtools to run R CMD check
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::check()