% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_gff_to_annotation.r
\name{parse_gff_to_annotation}
\alias{parse_gff_to_annotation}
\title{Parse GFF file to annotation format}
\usage{
parse_gff_to_annotation(gff_file, type_filter = c(".*RNA", ".*transcript"))
}
\arguments{
\item{gff_file}{Path to the GFF file or URL}

\item{type_filter}{Vector of regex patterns to filter the type column
(3rd column)}
}
\value{
A data frame with columns: txid, gid, type, symbol, tax_id
}
\description{
This function parses a GFF file and extracts gene annotation information.
It processes the attributes column (9th column) to extract transcript ID,
gene ID, gene symbol, and taxonomic ID.
}
