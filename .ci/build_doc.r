path <- getwd()

setwd(path)

t <- tryCatch(
  {
    roxygen2::roxygenize(package.dir = path)
    b_value <- "pass"
    b_color <- "green"
    devtools::build_manual(path = path)
    system("mv nexodiff_*.pdf documentation.pdf")
    return(c(b_value, b_color))
  },
  warning = function(w) {
    print("WARNING")
    print(w)
    b_value <- "warn"
    b_color <- "orange"
    return(c(b_value, b_color))
  },
  error = function(e) {
    print("WARNING")
    print(e)
    b_value <- "fail"
    b_color <- "red"
    system("touch documentation.pdf")
    return(c(b_value, b_color))
  }
)

write(
      jinjar::render(
        fs::path(".ci/badge.svg.j2"),
        value = t[1],
        type = "doc",
        color = t[2]),
      file = "badge_doc.svg"
    )
