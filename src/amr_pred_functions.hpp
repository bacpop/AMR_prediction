#pragma once
#include <string>
#include <utility>
#include <sdsl/suffix_arrays.hpp>
#include "./Matrix.h"

std::pair<sdsl::csa_wt<>,bool> create_index(std::string filepath); 
char complement(char n);
std::string invert(std::string seq);
bool check_species(const sdsl::csa_wt<>& fm_index);
Numeric_lib::Matrix<double,1> lookup_unitigs(const sdsl::csa_wt<>& fm_index, const std::vector<std::string>& unitigs);
std::string make_prediction(std::string assembly_filename);
std::string make_prediction_csv(std::string assembly_filename, std::vector<std::string> antibiotics);
std::string make_prediction_json(std::string assembly_filename);
