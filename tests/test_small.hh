#pragma once

#include "setup_tests.hh"
#include "globals.hh"
#include "Kmer.hh"
#include "NodeBOSS.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include "SubsetWT.hh"
#include "libwheeler/BOSS_builder.hh"
#include "suffix_group_optimization.hh"
#include <gtest/gtest.h>


typedef long long LL;
typedef Kmer<MAX_KMER_LENGTH> kmer_t;

typedef NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss_t;

typedef NodeBOSS<SubsetSplitRank<sdsl::bit_vector, sdsl::rank_support_v5<>,
                  sdsl::bit_vector, sdsl::rank_support_v5<>>> splitboss_t;

typedef NodeBOSS<SubsetWT<sdsl::wt_blcd<sdsl::bit_vector,
                                    sdsl::rank_support_v5<>,
                                    sdsl::select_support_scan<1>,
                                    sdsl::select_support_scan<0>>>
             > subsetwtboss_t;

typedef NodeBOSS<SubsetConcatRank<sdsl::bit_vector,
                        sdsl::bit_vector::select_0_type,
                        sdsl::wt_blcd<sdsl::bit_vector,
                                    sdsl::rank_support_v5<>,
                                    sdsl::select_support_scan<1>,
                                    sdsl::select_support_scan<0>>>
        > concatboss_t;


class TEST_SMALL : public ::testing::Test {
    protected:

    vector<string> input;
    set<string> kmers;
    LL k = 3;
    matrixboss_t matrixboss;
    splitboss_t splitboss;
    subsetwtboss_t subsetWTboss;
    concatboss_t concatboss;
    sdsl::bit_vector suffix_group_starts;
    sdsl::bit_vector A_bits, C_bits, G_bits, T_bits;

    void SetUp() override {
        input = {"TAGCAAGCACAGCATACAGA"};
        k = 3;
        kmers = get_all_kmers(input, k);

        BOSS_builder<BOSS<sdsl::bit_vector>, Kmer_stream_in_memory> bb;
        Kmer_stream_in_memory stream(input, k+1);
        BOSS<sdsl::bit_vector> wheelerBOSS = bb.build(stream);
        matrixboss.build_from_WheelerBOSS(wheelerBOSS);
        A_bits = matrixboss.subset_rank.A_bits;
        C_bits = matrixboss.subset_rank.C_bits;
        G_bits = matrixboss.subset_rank.G_bits;
        T_bits = matrixboss.subset_rank.T_bits;

        // Push bits to the left end of a suffix group
        suffix_group_starts = mark_suffix_groups(A_bits, C_bits, G_bits, T_bits, matrixboss.C, matrixboss.k);
        push_bits_left(A_bits, C_bits, G_bits, T_bits, suffix_group_starts);

        matrixboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        splitboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        subsetWTboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        concatboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
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

// Queries all 4^k k-mers and checks that the membership queries give the right answers
template<typename nodeboss_t>
void check_all_queries(const nodeboss_t& nodeboss, const set<string>& true_kmers){
    for(uint64_t mask = 0; mask < (1 << (2*nodeboss.k)); mask++){
        string kmer;
        for(int64_t i = 0; i < nodeboss.k; i++){
            if(((mask >> 2*i) & 0x3) == 0) kmer += 'A';
            if(((mask >> 2*i) & 0x3) == 1) kmer += 'C';
            if(((mask >> 2*i) & 0x3) == 2) kmer += 'G';
            if(((mask >> 2*i) & 0x3) == 3) kmer += 'T';
        }
        bool is_found = true_kmers.count(kmer); // Truth
        int64_t colex = nodeboss.search(kmer);
        if(is_found) ASSERT_GE(colex, 0); else ASSERT_EQ(colex, -1);
        logger << kmer << " " << colex << endl;
    }
}

TEST(TEST_IM_CONSTRUCTION, small_hand_crafted){
    matrixboss_t X;
    vector<string> strings = {"CCCGTGATGGCTA", "TAATGCTGTAGC", "TGGCTCGTGTAGTCGA"};
    X.build_from_strings(strings, 4);
    set<string> true_kmers = get_all_kmers(strings, 4);
    logger << "Queries on in-memory constructed matrixboss" << endl;
    check_all_queries(X, true_kmers);
}

TEST(TEST_IM_CONSTRUCTION, redundant_dummies){
    matrixboss_t X;
    vector<string> strings = {"AAAA", "ACCC", "ACCG", "CCCG", "TTTT"};
    X.build_from_strings(strings, 4);
    set<string> true_kmers = get_all_kmers(strings, 4);
    logger << "Queries on in-memory constructed matrixboss" << endl;
    check_all_queries(X, true_kmers);
    ASSERT_EQ(X.n_nodes, 9); // Dummies C, CC and CCC should not be there.
}

TEST_F(TEST_SMALL, matrixboss){
    check_all_queries(matrixboss, kmers);
}

TEST_F(TEST_SMALL, splitboss){
    check_all_queries(splitboss, kmers);
}

TEST_F(TEST_SMALL, subsetwtboss){
    check_all_queries(subsetWTboss, kmers);
}

TEST_F(TEST_SMALL, concatboss){
    check_all_queries(concatboss, kmers);
}

