% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{remove_id_version_suffix}
\alias{remove_id_version_suffix}
\title{Remove version suffixes from a vector of IDs}
\usage{
remove_id_version_suffix(ids, suffix_pattern = "\\\\.\\\\d+$")
}
\arguments{
\item{ids}{A character vector of identifiers.}

\item{suffix_pattern}{A string representing the regular expression pattern
for the suffix to remove. Defaults to \code{"\\\\.\\\\d+$"} (a hyphen followed by one
or more digits at the end of the string). If \code{NULL} or \code{FALSE} or an empty
string, the original IDs are returned unmodified.}
}
\value{
A character vector of IDs with the specified suffixes removed.
Issues a warning if new duplicates are created from previously unique
suffixed IDs.
}
\description{
This function takes a character vector of IDs and removes a specified suffix
pattern. It also checks if removing the suffix results in new duplicate IDs
that were not duplicates before the suffix removal.
}
\examples{
ids_to_clean <- c("ID-1", "ID-2.1", "ID-2.2", "ID-3", "ID-4.01", "ID-4.001")
remove_id_version_suffix(ids_to_clean, suffix_pattern = "\\\\.\\\\d+$")
# c("ID-1", "ID-2", "ID-2", "ID-3", "ID-4", "ID-4")
# Warning: Removing ID version suffixes created new duplicate IDs...

remove_id_version_suffix(ids_to_clean, suffix_pattern = "-\\\\d+$")
# c("ID", "ID-2.1", "ID-2.2", "ID", "ID-4.01", "ID-4.001")
# Warning: Removing ID version suffixes created new duplicate IDs...

remove_id_version_suffix(c("TX.1", "TX.2", "TX.3"), suffix_pattern = NULL)
# c("TX.1", "TX.2", "TX.3")
}
