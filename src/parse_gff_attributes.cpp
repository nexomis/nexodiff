#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame parse_gff_attributes(CharacterVector attr_strings, CharacterVector types) {
  int n = attr_strings.size();
  CharacterVector tx_id(n), gene_id(n), gene_symbol(n);

  for (int i = 0; i < n; ++i) {
    std::string attr_string = Rcpp::as<std::string>(attr_strings[i]);
    std::string current_key, current_value;

    std::stringstream ss(attr_string);
    std::string token;

    while (getline(ss, token, ';')) {
      std::stringstream token_ss(token);
      std::string key;
      if (getline(token_ss, key, '=')) {
        std::string value;
        if (key == "ID") {
          getline(token_ss, value, ';');
          tx_id[i] = value;
        } else if (key == "Parent") {
          getline(token_ss, value, ';');
          gene_id[i] = value;
        } else if (key == "gene") {
          getline(token_ss, value, ';');
          gene_symbol[i] = value;
        } else if (key == "Dbxref") {
          std::string dbxref_value;
          getline(token_ss, dbxref_value, ';');
          std::stringstream ss2(dbxref_value);
          std::string token2;
          while (getline(ss2, token2, ',')) {
            std::stringstream token_ss2(token2);
            std::string key2;
            std::string value2;
            if (getline(token_ss2, key2, ':')) {
              if (getline(token_ss2, value2, ',')) {
                if (key2 == "GenBank") {
                  tx_id[i] = value2;
                } else if (key2 == "GeneID") {
                  gene_id[i] = value2;
                }
              }
            }
          }
        } else if (key == "gbkey") {
          getline(token_ss, value, ';');
          types[i] = value;
        }
      }
    }
  }

  return DataFrame::create(Named("txid") = tx_id,
                           Named("gid") = gene_id,
                           Named("type") = types,
                           Named("symbol") = gene_symbol
                           );
}
