#pragma once

#include "setup_tests.hh"
#include "globals.hh"
#include "Kmer.hh"
#include "variants.hh"
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

TEST(TEST_IM_CONSTRUCTION, test_all_variants){
    vector<string> strings = {"CCCGTGATGGCTA", "TAATGCTGTAGC", "TGGCTCGTGTAGTCGA"};
    LL k = 4;
    set<string> true_kmers = get_all_kmers(strings, k);

    vector<string> filenames;
    for(LL i = 0; i < 10; i++){ // Create temp file for each of the 10 variants
        filenames.push_back(get_temp_file_manager().create_filename());
    }

    // Build
    {
        plain_matrix_sbwt_t v1;
        rrr_matrix_sbwt_t v2;
        mef_matrix_sbwt_t v3;
        plain_split_sbwt_t v4;
        rrr_split_sbwt_t v5;
        mef_split_sbwt_t v6;
        plain_concat_sbwt_t v7;
        mef_concat_sbwt_t v8;
        plain_sswt_sbwt_t v9;
        rrr_sswt_sbwt_t v10;

        v1.build_from_strings(strings, k);
        v2.build_from_strings(strings, k);
        v3.build_from_strings(strings, k);
        v4.build_from_strings(strings, k);
        v5.build_from_strings(strings, k);
        v6.build_from_strings(strings, k);
        v7.build_from_strings(strings, k);
        v8.build_from_strings(strings, k);
        v9.build_from_strings(strings, k);
        v10.build_from_strings(strings, k);

        v1.serialize(filenames[0]);
        v2.serialize(filenames[1]);
        v3.serialize(filenames[2]);
        v4.serialize(filenames[3]);
        v5.serialize(filenames[4]);
        v6.serialize(filenames[5]);
        v7.serialize(filenames[6]);
        v8.serialize(filenames[7]);
        v9.serialize(filenames[8]);
        v10.serialize(filenames[9]);
    }

    // Load and query
    {
        plain_matrix_sbwt_t v1;
        rrr_matrix_sbwt_t v2;
        mef_matrix_sbwt_t v3;
        plain_split_sbwt_t v4;
        rrr_split_sbwt_t v5;
        mef_split_sbwt_t v6;
        plain_concat_sbwt_t v7;
        mef_concat_sbwt_t v8;
        plain_sswt_sbwt_t v9;
        rrr_sswt_sbwt_t v10;

        v1.load(filenames[0]);
        v2.load(filenames[1]);
        v3.load(filenames[2]);
        v4.load(filenames[3]);
        v5.load(filenames[4]);
        v6.load(filenames[5]);
        v7.load(filenames[6]);
        v8.load(filenames[7]);
        v9.load(filenames[8]);
        v10.load(filenames[9]);

        check_all_queries(v1, true_kmers);
        check_all_queries(v2, true_kmers);
        check_all_queries(v3, true_kmers);
        check_all_queries(v4, true_kmers);
        check_all_queries(v5, true_kmers);
        check_all_queries(v6, true_kmers);
        check_all_queries(v7, true_kmers);
        check_all_queries(v8, true_kmers);
        check_all_queries(v9, true_kmers);
        check_all_queries(v10, true_kmers);
    }
}

TEST(TEST_IM_CONSTRUCTION, redundant_dummies){
    plain_matrix_sbwt_t X;
    vector<string> strings = {"AAAA", "ACCC", "ACCG", "CCCG", "TTTT"};
    X.build_from_strings(strings, 4);
    set<string> true_kmers = get_all_kmers(strings, 4);
    logger << "Queries on in-memory constructed matrixboss" << endl;
    check_all_queries(X, true_kmers);
    ASSERT_EQ(X.n_nodes, 9); // Dummies C, CC and CCC should not be there.
}

TEST(TEST_IM_CONSTRUCTION, not_full_alphabet){
    plain_matrix_sbwt_t X;
    vector<string> strings = {"AAAA", "ACCC", "ACCG", "CCCG"}; // No 'T' exists
    X.build_from_strings(strings, 4);
    set<string> true_kmers = get_all_kmers(strings, 4);
    logger << "Queries on in-memory constructed matrixboss" << endl;
    check_all_queries(X, true_kmers);
}
