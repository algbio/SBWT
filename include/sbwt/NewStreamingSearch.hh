#pragma once

#include <vector>
#include <utility>
#include <string>
#include "SBWT.hh"
#include "variants.hh"
#include "sdsl/vectors.hpp"

using namespace std;
using namespace sbwt;

// Returns pairs (a_1, b_1), (a_2, b_2)..., such that 
// - a_i is the length of the longest match ending at query[i]
// - b_i is the colex rank of one arbitrary k-mer having the longest match.
vector<pair<int64_t,int64_t> > new_streaming_search(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& query){
    const sdsl::bit_vector& A_bits = sbwt.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = sbwt.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = sbwt.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = sbwt.get_subset_rank_structure().T_bits;

    const sdsl::rank_support_v5<>& A_bits_rs = sbwt.get_subset_rank_structure().A_bits_rs;
    const sdsl::rank_support_v5<>& C_bits_rs = sbwt.get_subset_rank_structure().C_bits_rs;
    const sdsl::rank_support_v5<>& G_bits_rs = sbwt.get_subset_rank_structure().G_bits_rs;
    const sdsl::rank_support_v5<>& T_bits_rs = sbwt.get_subset_rank_structure().T_bits_rs;

    const vector<int64_t>& C = sbwt.get_C_array();

    const int64_t k = sbwt.get_k();
}