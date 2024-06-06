# script file to read
# usage generate_data.r folder
# folder whould contains those files:
# - design.csv
# - expr.csv
# - annotation.csv
# - id_mapping.tab

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("only one argument expected: folder")
}

#folder <- args[1]
folder <- "./tests/sim_inputs/set1"
if (! dir.exists(folder)){
  stop("given arguments does not match an existing directory")
}

annot <- read.delim(paste(folder, "annotation.txt", sep = "/"), sep = " ",
  header = FALSE)

transcript_ids <- annot[, 1]

expr <- read.delim(paste(folder, "expr.csv", sep = "/"), sep = ";")

for (sample_name in names(expr)) {
  file_name <- paste(folder, "/", sample_name, ".h5", sep = "")
  if (file.exists(file_name)) {
    file.remove(file_name)
  }
  rhdf5::h5createFile(file_name)
  rhdf5::h5createGroup(file_name,"aux")
  rhdf5::h5write(expr[, sample_name], file_name, "est_counts")
  rhdf5::h5write(rep(1000, length(transcript_ids)), file_name,
    "aux/eff_lengths")
  rhdf5::h5write(transcript_ids, file_name, "aux/ids")
}
