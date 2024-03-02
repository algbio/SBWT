#pragma once

#include "setup_tests.hh"
#include "kmc_construct.hh"
#include "globals.hh"
#include "Kmer.hh"
#include "variants.hh"
#include "SBWT.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include "SubsetMatrixSelectSupport.hh"
#include "SubsetWT.hh"
#include "suffix_group_optimization.hh"
#include <gtest/gtest.h>
#include <set>

using namespace sbwt;


typedef Kmer<MAX_KMER_LENGTH> kmer_t;

// Queries all 4^k k-mers and checks that the membership queries give the right answers
template<typename nodeboss_t>
void check_all_queries(const nodeboss_t& nodeboss, const set<string>& true_kmers){
    for(uint64_t mask = 0; mask < (1 << (2*nodeboss.get_k())); mask++){
        string kmer;
        for(int64_t i = 0; i < nodeboss.get_k(); i++){
            if(((mask >> 2*i) & 0x3) == 0) kmer += 'A';
            if(((mask >> 2*i) & 0x3) == 1) kmer += 'C';
            if(((mask >> 2*i) & 0x3) == 2) kmer += 'G';
            if(((mask >> 2*i) & 0x3) == 3) kmer += 'T';
        }
        bool is_found = true_kmers.count(kmer); // Truth
        int64_t column = nodeboss.search(kmer);
        if(is_found) ASSERT_GE(column, 0); else ASSERT_EQ(column, -1);
        //logger << kmer << " " << colex << endl;
    }

    // Check NN...N
    string NNN(nodeboss.get_k(), 'N');
    ASSERT_EQ(nodeboss.search(NNN), -1);
}

// Queries all 4^k k-mers and checks that the membership queries give the right answers
template<typename nodeboss_t>
void check_streaming_queries(const nodeboss_t& nodeboss, const set<string>& true_kmers, const string& input){
    vector<int64_t> result = nodeboss.streaming_search(input.c_str(), input.size());

    // Check
    for(int64_t i = 0; i < (int64_t)result.size(); i++){
        string kmer = input.substr(i, nodeboss.get_k());
        bool is_found = true_kmers.count(kmer); // Truth
        if(is_found) ASSERT_GE(result[i], 0); else ASSERT_EQ(result[i], -1);
        //logger << kmer << " " << result[i] << endl;
    }

    // Check NN...N
    string NNN(100, 'N');
    for(int64_t x : nodeboss.streaming_search(NNN)){
        ASSERT_EQ(x, -1);
    }
}

void run_small_testcase(const vector<string>& strings, int64_t k){
    plain_matrix_sbwt_t index_im;
    build_nodeboss_in_memory(strings, index_im, k, false); 

    string temp_filename = get_temp_file_manager().create_filename("", ".fna");
    write_seqs_to_fasta_file(strings, temp_filename);

    plain_matrix_sbwt_t::BuildConfig config;
    config.input_files = {temp_filename};
    config.k = k;
    config.build_streaming_support = false;
    config.ram_gigas = 2;
    config.n_threads = 2;
    config.min_abundance = 1;
    plain_matrix_sbwt_t index_kmc(config);

    logger << index_im.get_subset_rank_structure().A_bits << endl;
    logger << index_im.get_subset_rank_structure().C_bits << endl;
    logger << index_im.get_subset_rank_structure().G_bits << endl;
    logger << index_im.get_subset_rank_structure().T_bits << endl;
    logger << "--" << endl;
    logger << index_kmc.get_subset_rank_structure().A_bits << endl;
    logger << index_kmc.get_subset_rank_structure().C_bits << endl;
    logger << index_kmc.get_subset_rank_structure().G_bits << endl;
    logger << index_kmc.get_subset_rank_structure().T_bits << endl;

    ASSERT_EQ(index_im.get_subset_rank_structure().A_bits, index_kmc.get_subset_rank_structure().A_bits);
    ASSERT_EQ(index_im.get_subset_rank_structure().C_bits, index_kmc.get_subset_rank_structure().C_bits);
    ASSERT_EQ(index_im.get_subset_rank_structure().G_bits, index_kmc.get_subset_rank_structure().G_bits);
    ASSERT_EQ(index_im.get_subset_rank_structure().T_bits, index_kmc.get_subset_rank_structure().T_bits);

    set<string> true_kmers = get_all_kmers(strings, k);
    check_all_queries(index_im, true_kmers);
    check_all_queries(index_kmc, true_kmers);
}

template<typename nodeboss_t>
void test_partial_search(){
    nodeboss_t sbwt;
    vector<string> strings = {"CCCGTGATGGCTA", "TAATGCTGTAGC", "TGGCTCGTGTAGTCGA"};
    int64_t k = 6;
    build_nodeboss_in_memory(strings, sbwt, k, false); 
    string kmers_concat = sbwt.reconstruct_all_kmers();

    string x = "GCAAAA"; // GCA is the longest matching prefix
    pair<int64_t, int64_t> I; int64_t len;
    std::tie(I, len) = sbwt.partial_search(x);
    ASSERT_EQ(len, 2);

    cerr << I << " " << len << endl;

    // Check that the k-mers in the interval are exactly those that have GC as a suffix.
    for(int64_t i = 0; i < sbwt.number_of_subsets(); i++){
        string kmer = kmers_concat.substr(i*k, k);
        // cerr << kmer << endl;
        if(i >= I.first && i <= I.second) 
            ASSERT_EQ(kmer.substr(k-len), "GC");
        else
            ASSERT_NE(kmer.substr(k-len), "GC");
    }

}

struct TrivialSubsetSelect{

    // Precalculates and stores all answers

    map<char, vector<int64_t>> select_answers;

    template<typename sbwt_t>
    TrivialSubsetSelect(sbwt_t& sbwt) {
        for(char c : string("ACGT")){
            select_answers[c].push_back(-1); // select(0) is undefined
            for(int64_t i = 0; i < sbwt.number_of_subsets(); i++){
                if(sbwt.get_subset_rank_structure().contains(i,c)){
                    select_answers[c].push_back(i);
                }
            }
        }
    }

    int64_t select(int64_t i, char c) const{
        return select_answers.at(c)[i];
    }
};

template<typename nodeboss_t>
void test_get_kmer(){
    nodeboss_t sbwt;
    vector<string> strings = {"CCCGTGATGGCTA", "TAATGCTGTAGC", "TGGCTCGTGTAGTCGA"};
    int64_t k = 6;
    build_nodeboss_in_memory(strings, sbwt, k, false); 
    string kmers_concat = sbwt.reconstruct_all_kmers();

    TrivialSubsetSelect ss(sbwt);
    vector<char> buf(k);
    for(int64_t i = 0; i < sbwt.number_of_subsets(); i++){
        string true_kmer = kmers_concat.substr(i*k, k);

        sbwt.get_kmer(i, buf.data());
        string test_kmer = string(buf.data(), buf.data()+k);

        cerr << true_kmer << " " << test_kmer << endl;
        ASSERT_EQ(true_kmer, test_kmer);
    }

}

TEST(TEST_GET_KMER, fast){
    plain_matrix_sbwt_t sbwt;
    vector<string> strings = {"CCCGTGATGGCTA", "TAATGCTGTAGC", "TGGCTCGTGTAGTCGA"};
    int64_t k = 6;
    build_nodeboss_in_memory(strings, sbwt, k, false); 
    string kmers_concat = sbwt.reconstruct_all_kmers();

    TrivialSubsetSelect trivial_ss(sbwt);
    string ss_filename = get_temp_file_manager().create_filename();
    {
        SubsetMatrixSelectSupport<sdsl::bit_vector> nontrivial_ss(sbwt);
        throwing_ofstream out(ss_filename, ios::binary);
        nontrivial_ss.serialize(out.stream);
    } // Flushes stream and frees the select support

    // Load the select support back from disk
    SubsetMatrixSelectSupport<sdsl::bit_vector> nontrivial_ss;
    throwing_ifstream ss_in(ss_filename);
    nontrivial_ss.load(ss_in.stream, ios::binary);

    vector<char> buf_trivial(k);
    vector<char> buf_nontrivial(k);
    for(int64_t i = 0; i < sbwt.number_of_subsets(); i++){
        string true_kmer = kmers_concat.substr(i*k, k);

        sbwt.get_kmer_fast<TrivialSubsetSelect>(i, buf_trivial.data(), trivial_ss);
        string trivial_kmer = string(buf_trivial.data(), buf_trivial.data()+k);

        sbwt.get_kmer_fast<SubsetMatrixSelectSupport<sdsl::bit_vector>>(i, buf_nontrivial.data(), nontrivial_ss);
        string nontrivial_kmer = string(buf_nontrivial.data(), buf_nontrivial.data()+k);

        cerr << true_kmer << " " << trivial_kmer << " " << nontrivial_kmer << endl;
        ASSERT_EQ(true_kmer, trivial_kmer);
        ASSERT_EQ(true_kmer, nontrivial_kmer);
    }

}

TEST(TEST_PARTIAL_SEARCH, all){
    // mef variants are commented out because they don't compile because the mef bit vector
    // does not support access currently.

    test_partial_search<plain_matrix_sbwt_t>();
    test_partial_search<rrr_matrix_sbwt_t>();
    //test_partial_search<mef_matrix_sbwt_t>();
    test_partial_search<plain_split_sbwt_t>();
    test_partial_search<rrr_split_sbwt_t>();
    //test_partial_search<mef_split_sbwt_t>();
    test_partial_search<plain_concat_sbwt_t>();
    //test_partial_search<mef_concat_sbwt_t>();
    test_partial_search<plain_sswt_sbwt_t>();
    test_partial_search<rrr_sswt_sbwt_t>();
}


TEST(TEST_GET_KMER, all){
    // mef variants are commented out because they don't compile because the mef bit vector
    // does not support access currently.

    test_get_kmer<plain_matrix_sbwt_t>();
    test_get_kmer<rrr_matrix_sbwt_t>();
    //test_get_kmer<mef_matrix_sbwt_t>();
    test_get_kmer<plain_split_sbwt_t>();
    test_get_kmer<rrr_split_sbwt_t>();
    //test_get_kmer<mef_split_sbwt_t>();
    test_get_kmer<plain_concat_sbwt_t>();
    //test_get_kmer<mef_concat_sbwt_t>();
    test_get_kmer<plain_sswt_sbwt_t>();
    test_get_kmer<rrr_sswt_sbwt_t>();
}

TEST(TEST_KMC_CONSTRUCT, not_all_dummies_needed){
    vector<string> strings = {"CCCGTGATGGCTA", "TAATGCTGTAGC", "TGGCTCGTGTAGTCGA"};
    int64_t k = 4;
    run_small_testcase(strings, k);
}

TEST(TEST_KMC_CONSTRUCT, multiple_input_files){
    vector<string> strings = {"CCCGTGATGGCTA", "TAATGCTGTAGC", "TGGCTCGTGTAGTCGA"};
    int64_t k = 4;
    string f1 = get_temp_file_manager().create_filename("",".fna");
    string f2 = get_temp_file_manager().create_filename("",".fna");
    string f3 = get_temp_file_manager().create_filename("",".fna");
    string f123 = get_temp_file_manager().create_filename("",".fna");
    write_seqs_to_fasta_file({strings[0]}, f1);
    write_seqs_to_fasta_file({strings[1]}, f2);
    write_seqs_to_fasta_file({strings[2]}, f3);
    write_seqs_to_fasta_file(strings, f123);

    NodeBOSSKMCConstructor<plain_matrix_sbwt_t> X;
    plain_matrix_sbwt_t index1, index2;
    X.build({f123}, index1, k, 1, 2, true, 1, 1e9, 2);
    X.build({f1, f2, f3}, index2, k, 1, 2, true, 1, 1e9, 2);

    ASSERT_EQ(index1.get_subset_rank_structure().A_bits, index2.get_subset_rank_structure().A_bits);
    ASSERT_EQ(index1.get_subset_rank_structure().C_bits, index2.get_subset_rank_structure().C_bits);
    ASSERT_EQ(index1.get_subset_rank_structure().G_bits, index2.get_subset_rank_structure().G_bits);
    ASSERT_EQ(index1.get_subset_rank_structure().T_bits, index2.get_subset_rank_structure().T_bits);
    ASSERT_EQ(index1.get_streaming_support(), index2.get_streaming_support());
}


TEST(TEST_IM_CONSTRUCTION, redundant_dummies){
    vector<string> strings = {"AAAA", "ACCC", "ACCG", "CCCG", "TTTT"};

    run_small_testcase(strings, 4);

    logger << "Checking that there are no extra dummies" << endl;

    plain_matrix_sbwt_t X;
    build_nodeboss_in_memory(strings, X, 4, false);
    ASSERT_EQ(X.number_of_subsets(), 9); // Dummies C, CC and CCC should not be there.

    string filename = get_temp_file_manager().create_filename("", ".fna");
    write_seqs_to_fasta_file(strings, filename);
    plain_matrix_sbwt_t::BuildConfig config;
    config.input_files = {filename};
    config.k = 4;
    config.build_streaming_support = false;
    config.n_threads = 1;
    config.ram_gigas = 2;
    config.min_abundance = 1;
    plain_matrix_sbwt_t X2(config);
    ASSERT_EQ(X2.number_of_subsets(), 9); // Dummies C, CC and CCC should not be there.
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

TEST(TEST_IM_CONSTRUCTION, cyclic){
    vector<string> strings = {"ACGTACGTACGT"}; // String is cyclic
    run_small_testcase(strings, 3);
}


TEST(TEST_IM_CONSTRUCTION, test_serialization){
    vector<string> strings = {"CCCGTGATGGCTA", "TAATGCTGTAGC", "TGGCTCGTGTAGTCGA", "NNAAAAAAAAAAAA"}; // The last string tests N's with k-mer precalc of length 2
    int64_t k = 4;
    set<string> true_kmers = get_all_kmers(strings, k);
    set<string> true_kmers2; // Remove k-mers with N's
    for(const string& S : true_kmers){
        if(S.find('N',0) == string::npos)
            true_kmers2.insert(S);
    }
    true_kmers = true_kmers2;

    vector<string> filenames;
    for(int64_t i = 0; i < 10; i++){ // Create temp file for each of the 10 variants
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

        build_nodeboss_in_memory(strings, v1, k, true);
        build_nodeboss_in_memory(strings, v2, k, true);
        build_nodeboss_in_memory(strings, v3, k, true);
        build_nodeboss_in_memory(strings, v4, k, true);
        build_nodeboss_in_memory(strings, v5, k, true);
        build_nodeboss_in_memory(strings, v6, k, true);
        build_nodeboss_in_memory(strings, v7, k, true);
        build_nodeboss_in_memory(strings, v8, k, true);
        build_nodeboss_in_memory(strings, v9, k, true);
        build_nodeboss_in_memory(strings, v10, k, true);

        v1.do_kmer_prefix_precalc(2);

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