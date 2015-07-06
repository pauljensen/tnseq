#include <Rcpp.h>
using namespace Rcpp;

#include <unordered_map>
#include <iostream>
#include <fstream>

typedef std::unordered_map<std::string, int> inner_t;

std::unordered_map<std::string, inner_t> hashes;

// [[Rcpp::export]]
void build_hashes(CharacterVector outer_keys) {
  for (int i=0; i<outer_keys.size(); i++) {
    inner_t hash;
    hashes.insert(std::make_pair(outer_keys[i], hash));
  }
}

// [[Rcpp::export]]
void erase_hashes(CharacterVector outer_keys) {
  for (int i=0; i<outer_keys.size(); i++) {
    hashes.erase(Rcpp::as<std::string>(outer_keys[i]));
  }
}

// [[Rcpp::export]]
void add_keys(std::string outer, CharacterVector keys) {
  int n = keys.size();
  for (int i=0; i<n; i++) {
    hashes[outer][Rcpp::as<std::string>(keys[i])] += 1;

  }
}

// [[Rcpp::export]]
NumericVector get_values(std::string outer, CharacterVector keys) {
  int n = keys.size();
  NumericVector values(n);
  for (int i=0; i<n; i++) {
    values[i] = hashes[outer][Rcpp::as<std::string>(keys[i])];
  }
  return values;
}

// Dumps a hash full of reads into a fasta file, with format
//
//    >i-n
//    actgacgt...
//
// where i is a unique (increasing) integer, and n is the number of reads.
// Returns the total number of reads in the file (sum(n) forall i).
// [[Rcpp::export]]
int write_hash_to_file(std::string outer_key, std::string filename) {
  std::ofstream myfile;
  myfile.open(filename);
  inner_t inner = hashes[outer_key];
  int count = 0;
  int readcount = 0;
  for (inner_t::iterator it = inner.begin(); it != inner.end(); ++it) {
    myfile << ">" << ++count << "-" << it->second << std::endl;
    myfile << it->first << std::endl;
    readcount += it->second;
  } 
  myfile.close();
  return readcount;
}
