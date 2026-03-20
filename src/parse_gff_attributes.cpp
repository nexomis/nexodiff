#include <Rcpp.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <utility>
#include <regex>
using namespace Rcpp;

// Helper function to strip prefix from an identifier
std::string strip_prefix(const std::string& id, const std::string& prefix) {
  if (id.substr(0, prefix.size()) == prefix) {
    return id.substr(prefix.size());
  }
  return id;
}

// Helper function to strip numerical suffix (e.g., -202, -888) from a string
// Returns the string without the trailing "-[0-9]+" pattern if present,
// otherwise returns the original string unchanged
std::string strip_numerical_suffix(const std::string& str) {
  std::regex suffix_pattern("-[0-9]+$");
  return std::regex_replace(str, suffix_pattern, "");
}

// [[Rcpp::export]]
DataFrame parse_gff_attributes(
    CharacterVector gff_attributes,
    CharacterVector feature_types,
    CharacterVector annotation_style = "ncbi"
) {
  int num_transcripts = gff_attributes.size();
  CharacterVector transcript_ids(num_transcripts),
                  gene_ids(num_transcripts),
                  gene_symbols(num_transcripts);

  std::string style_str = Rcpp::as<std::string>(annotation_style[0]);
  bool is_ensembl = (style_str == "ensembl");

  for (int idx = 0; idx < num_transcripts; ++idx) {
    std::string attributes = Rcpp::as<std::string>(gff_attributes[idx]);

    std::string transcript_id_value, gene_id_value, gene_symbol_value;

    std::stringstream attributes_stream(attributes);
    std::string attribute;

    while (getline(attributes_stream, attribute, ';')) {
      std::stringstream attribute_stream(attribute);
      std::string attribute_key;
      if (getline(attribute_stream, attribute_key, '=')) {
        std::string attribute_value;
        getline(attribute_stream, attribute_value, ';');

        if (is_ensembl) {
          // Ensembl style parsing
          if (attribute_key == "ID") {
            // Strip "transcript:" prefix if present
            transcript_id_value = strip_prefix(attribute_value, "transcript:");
          } else if (attribute_key == "transcript_id") {
            transcript_id_value = attribute_value;
          } else if (attribute_key == "gene_id") {
            gene_id_value = attribute_value;
          } else if (attribute_key == "Parent") {
            // Strip "gene:" prefix if present
            gene_id_value = strip_prefix(attribute_value, "gene:");
          } else if (attribute_key == "Name") {
            // Strip numerical suffix (e.g., -201) from gene symbol
            gene_symbol_value = strip_numerical_suffix(attribute_value);
          }
        } else {
          // NCBI style parsing (default)
          if (attribute_key == "ID") {
            transcript_ids[idx] = attribute_value;
          } else if (attribute_key == "Parent") {
            gene_ids[idx] = attribute_value;
          } else if (attribute_key == "gene") {
            gene_symbols[idx] = attribute_value;
          } else if (attribute_key == "Dbxref") {
            std::stringstream dbxref_stream(attribute_value);
            std::string dbxref_token;
            while (getline(dbxref_stream, dbxref_token, ',')) {
              std::stringstream dbxref_attribute_stream(dbxref_token);
              std::string dbxref_key;
              std::string dbxref_value;
              if (getline(dbxref_attribute_stream, dbxref_key, ':')) {
                if (getline(dbxref_attribute_stream, dbxref_value, ',')) {
                  if (dbxref_key == "GenBank") {
                    transcript_ids[idx] = dbxref_value;
                  } else if (dbxref_key == "GeneID") {
                    gene_ids[idx] = dbxref_value;
                  }
                }
              }
            }
          } else if (attribute_key == "gbkey") {
            feature_types[idx] = attribute_value;
          }
        }
      }
    }

    // Set values for Ensembl style at end of parsing
    if (is_ensembl) {
      transcript_ids[idx] = transcript_id_value;
      gene_ids[idx] = gene_id_value;
      gene_symbols[idx] = gene_symbol_value;
    }
  }

  return DataFrame::create(
      Named("transcript_id") = transcript_ids,
      Named("gene_id") = gene_ids,
      Named("feature_type") = feature_types,
      Named("gene_symbol") = gene_symbols
  );
}
