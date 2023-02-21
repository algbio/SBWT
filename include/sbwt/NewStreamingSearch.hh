#pragma once

#include <vector>
#include <utility>
#include <string>
#include "SBWT.hh"
#include "variants.hh"
#include "sdsl/vectors.hpp"

using namespace std;
using namespace sbwt;

vector<pair<int64_t,int64_t> > new_streaming_search(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& query){
    // Returns pairs (a_1, b_1), (a_2, b_2)..., such that 
    // - a_i is the length of the longest match ending at query[i]
    // - b_i is the colex rank of one arbitrary k-mer having the longest match.
}