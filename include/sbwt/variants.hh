#pragma once

#include <string>
#include <cstring>
#include "cxxopts.hpp"
#include "globals.hh"
#include "SBWT.hh"
#include "SubsetWT.hh"
#include "stdlib_printing.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include <filesystem>
#include "MEF.hpp"

namespace sbwt{

// matrices
typedef SBWT<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> plain_matrix_sbwt_t;
typedef SBWT<SubsetMatrixRank<sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type>> rrr_matrix_sbwt_t;
typedef SBWT<SubsetMatrixRank<mod_ef_vector<>, mod_ef_vector<>::rank_1_type>> mef_matrix_sbwt_t;

// splits
typedef SBWT<SubsetSplitRank<sdsl::bit_vector, sdsl::rank_support_v5<>,
                            sdsl::bit_vector, sdsl::rank_support_v5<>>> plain_split_sbwt_t;

typedef SBWT<SubsetSplitRank<sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type,
                            sdsl::bit_vector, sdsl::rank_support_v5<>>> rrr_split_sbwt_t;


typedef SBWT<SubsetSplitRank<mod_ef_vector<>, mod_ef_vector<>::rank_1_type,
                            sdsl::bit_vector, sdsl::rank_support_v5<>>> mef_split_sbwt_t;

// concats
typedef SBWT<SubsetConcatRank<sdsl::bit_vector,
                            sdsl::bit_vector::select_0_type,
                            sdsl::wt_blcd<sdsl::bit_vector,
                                        sdsl::rank_support_v5<>,
                                        sdsl::select_support_scan<1>,
                                        sdsl::select_support_scan<0>>>
            > plain_concat_sbwt_t;

typedef SBWT<SubsetConcatRank<sd_vector<>,
                            sd_vector<>::select_0_type,
                            sdsl::wt_blcd<rrr_vector<63>,
                                        rrr_vector<>::rank_1_type,
                                        rrr_vector<>::select_1_type,
                                        rrr_vector<>::select_0_type>>
            > mef_concat_sbwt_t;

// wavelet trees
typedef SBWT<SubsetWT<sdsl::wt_blcd<sdsl::bit_vector,
                                sdsl::rank_support_v5<>,
                                sdsl::select_support_scan<1>,
                                sdsl::select_support_scan<0>>>
            > plain_sswt_sbwt_t;


typedef SBWT<SubsetWT<sdsl::wt_blcd<sdsl::rrr_vector<>,
                                sdsl::rrr_vector<>::rank_1_type,
                                rrr_vector<>::select_1_type,
                                rrr_vector<>::select_0_type>>
            > rrr_sswt_sbwt_t;

}
