#pragma once

#include "setup_tests.hh"
#include "kmc_construct.hh"
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
        //logger << kmer << " " << colex << endl;
    }
}

// Queries all 4^k k-mers and checks that the membership queries give the right answers
template<typename nodeboss_t>
void check_streaming_queries(const nodeboss_t& nodeboss, const set<string>& true_kmers, const string& input){
    vector<int64_t> result = nodeboss.streaming_search(input.c_str(), input.size());

    // Check
    for(int64_t i = 0; i < (LL)result.size(); i++){
        string kmer = input.substr(i, nodeboss.k);
        bool is_found = true_kmers.count(kmer); // Truth
        if(is_found) ASSERT_GE(result[i], 0); else ASSERT_EQ(result[i], -1);
        //logger << kmer << " " << result[i] << endl;
    }
}

void run_small_testcase(const vector<string>& strings, LL k){
    vector<string> reverse_strings;
    for(const string& S : strings){
        string R(S.rbegin(), S.rend());
        reverse_strings.push_back(R);
    }

    plain_matrix_sbwt_t index_im;
    plain_matrix_sbwt_t index_kmc;

    index_im.build_from_strings(strings, k, false);

    string temp_filename = get_temp_file_manager().create_filename("", ".fna");
    write_seqs_to_fasta_file(reverse_strings, temp_filename);

    NodeBOSSKMCConstructor<plain_matrix_sbwt_t> X;
    X.build(temp_filename, index_kmc, k, 1, 2, false, 1);

    logger << index_im.subset_rank.A_bits << endl;
    logger << index_im.subset_rank.C_bits << endl;
    logger << index_im.subset_rank.G_bits << endl;
    logger << index_im.subset_rank.T_bits << endl;
    logger << "--" << endl;
    logger << index_kmc.subset_rank.A_bits << endl;
    logger << index_kmc.subset_rank.C_bits << endl;
    logger << index_kmc.subset_rank.G_bits << endl;
    logger << index_kmc.subset_rank.T_bits << endl;

    ASSERT_EQ(index_im.subset_rank.A_bits, index_kmc.subset_rank.A_bits);
    ASSERT_EQ(index_im.subset_rank.C_bits, index_kmc.subset_rank.C_bits);
    ASSERT_EQ(index_im.subset_rank.G_bits, index_kmc.subset_rank.G_bits);
    ASSERT_EQ(index_im.subset_rank.T_bits, index_kmc.subset_rank.T_bits);

    set<string> true_kmers = get_all_kmers(strings, k);
    check_all_queries(index_im, true_kmers);
    check_all_queries(index_kmc, true_kmers);
}

TEST(TEST_KMC_CONSTRUCT, not_all_dummies_needed){
    vector<string> strings = {"CCCGTGATGGCTA", "TAATGCTGTAGC", "TGGCTCGTGTAGTCGA"};
    LL k = 4;
    run_small_testcase(strings, k);
}


TEST(TEST_IM_CONSTRUCTION, redundant_dummies){
    vector<string> strings = {"AAAA", "ACCC", "ACCG", "CCCG", "TTTT"};

    run_small_testcase(strings, 4);

    logger << "Checking that there are no extra dummies" << endl;

    plain_matrix_sbwt_t X;
    X.build_from_strings(strings, 4, false);
    ASSERT_EQ(X.n_nodes, 9); // Dummies C, CC and CCC should not be there.

    plain_matrix_sbwt_t X2;
    string filename = get_temp_file_manager().create_filename("", ".fna");
    write_seqs_to_fasta_file(strings, filename);
    X2.build_using_KMC(filename, 4, false, 1, 2, 1);
    ASSERT_EQ(X2.n_nodes, 9); // Dummies C, CC and CCC should not be there.
}

TEST(TEST_IM_CONSTRUCTION, not_full_alphabet){
    vector<string> strings = {"AAAA", "ACCC", "ACCG", "CCCG"}; // No 'T' exists
    run_small_testcase(strings, 3);
}

TEST(TEST_IM_CONSTRUCTION, lots_of_dummies){
    vector<string> strings;
    for(int64_t i = 0; i < 20; i++){
        strings.push_back(generate_random_kmer(6));
    }
    run_small_testcase(strings, 6);
}



TEST(TEST_IM_CONSTRUCTION, test_serialization){
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

        v1.build_from_strings(strings, k, true);
        v2.build_from_strings(strings, k, true);
        v3.build_from_strings(strings, k, true);
        v4.build_from_strings(strings, k, true);
        v5.build_from_strings(strings, k, true);
        v6.build_from_strings(strings, k, true);
        v7.build_from_strings(strings, k, true);
        v8.build_from_strings(strings, k, true);
        v9.build_from_strings(strings, k, true);
        v10.build_from_strings(strings, k, true);

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

        vector<string> streaming_query_inputs = strings; // input strings
        streaming_query_inputs.push_back(generate_random_kmer(100));

        for(const string& S : streaming_query_inputs){
            check_streaming_queries(v1, true_kmers, S);
            check_streaming_queries(v2, true_kmers, S);
            check_streaming_queries(v3, true_kmers, S);
            check_streaming_queries(v4, true_kmers, S);
            check_streaming_queries(v5, true_kmers, S);
            check_streaming_queries(v6, true_kmers, S);
            check_streaming_queries(v7, true_kmers, S);
            check_streaming_queries(v8, true_kmers, S);
            check_streaming_queries(v9, true_kmers, S);
            check_streaming_queries(v10, true_kmers, S);
        }
    }
}