#pragma once

#include "sdsl/bit_vectors.hpp"
#include <iostream>

namespace sbwt{

using namespace std;

// Entropy of distribution P
double entropy(const vector<double>& P);

// Pushes the bits to the left end of the suffix group
void push_bits_left(sdsl::bit_vector& A_bits, 
                    sdsl::bit_vector& C_bits, 
                    sdsl::bit_vector& G_bits, 
                    sdsl::bit_vector& T_bits, 
                    const sdsl::bit_vector& suffix_group_marks);

// Maximally spreads the bits inside a suffix group. Assumes the bits have already been pushed to the left
void spread_bits_after_push_left(sdsl::bit_vector& A_bits,
                                 sdsl::bit_vector& C_bits, 
                                 sdsl::bit_vector& G_bits, 
                                 sdsl::bit_vector& T_bits, 
                                 const sdsl::bit_vector& suffix_group_marks);

sdsl::bit_vector mark_suffix_groups(const sdsl::bit_vector& A_bits,
                                    const sdsl::bit_vector& C_bits, 
                                    const sdsl::bit_vector& G_bits, 
                                    const sdsl::bit_vector& T_bits,
                                    int64_t k);

double compute_column_entropy(const sdsl::bit_vector& A_bits,
                              const sdsl::bit_vector& C_bits, 
                              const sdsl::bit_vector& G_bits, 
                              const sdsl::bit_vector& T_bits);

}