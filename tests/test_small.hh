#pragma once

#include "setup_tests.hh"
#include "globals.hh"
#include "Kmer.hh"
#include "NodeBOSS.hh"
#include "SubsetMatrixRank.hh"
#include "libwheeler/BOSS_builder.hh"
#include "suffix_group_optimization.hh"
#include <gtest/gtest.h>


typedef long long LL;
typedef Kmer<32> kmer_t;

class TEST_SMALL : public ::testing::Test {
    protected:

    vector<string> input;
    set<kmer_t> kmers;
    LL k = 3;
    NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss;
    sdsl::bit_vector suffix_group_starts;
    sdsl::bit_vector A_bits, C_bits, G_bits, T_bits;

    void SetUp() override {
        input = {"TAGCAAGCACAGCATACAGA"};
        k = 3;
        BOSS_builder<BOSS<sdsl::bit_vector>, Kmer_stream_in_memory> bb;
        Kmer_stream_in_memory stream(input, k+1);
        BOSS<sdsl::bit_vector> wheelerBOSS = bb.build(stream, 1e6, 1);
        matrixboss.build_from_WheelerBOSS(wheelerBOSS);
        A_bits = matrixboss.subset_rank.A_bits;
        C_bits = matrixboss.subset_rank.C_bits;
        G_bits = matrixboss.subset_rank.G_bits;
        T_bits = matrixboss.subset_rank.T_bits;
        suffix_group_starts = mark_suffix_groups(A_bits, C_bits, G_bits, T_bits, matrixboss.C, matrixboss.k);
        push_bits_left(A_bits, C_bits, G_bits, T_bits, suffix_group_starts);
    }

};

TEST_F(TEST_SMALL, check_construction){
    // Matrix rows from the example in the paper
    sdsl::bit_vector true_A_bits = {0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1};
    sdsl::bit_vector true_C_bits = {0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0};
    sdsl::bit_vector true_G_bits = {0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    sdsl::bit_vector true_T_bits = {1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    ASSERT_EQ(true_A_bits, A_bits);
    ASSERT_EQ(true_C_bits, C_bits);
    ASSERT_EQ(true_G_bits, G_bits);
    ASSERT_EQ(true_T_bits, T_bits);
}