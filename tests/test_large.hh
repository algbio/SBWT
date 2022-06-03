#pragma once

#include "setup_tests.hh"
#include "globals.hh"
#include "Kmer.hh"
#include "SBWT.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include "SubsetWT.hh"
#include "SeqIO.hh"
#include "suffix_group_optimization.hh"
#include "variants.hh"
#include <gtest/gtest.h>
#include <unordered_set>

using namespace sbwt;

typedef long long LL;
typedef Kmer<MAX_KMER_LENGTH> kmer_t;

typedef SBWT<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss_t;

class TEST_LARGE : public ::testing::Test {
    protected:

    static matrixboss_t matrixboss;
    static matrixboss_t matrixboss_reference;
    static LL k;
    static vector<string> seqs;

    static void SetUpTestSuite(){
        string filename = "example_data/coli3.fna";

        vector<string> seqs;
        SeqIO::Unbuffered_Reader sr(filename);
        while(!sr.done()){
            string S = sr.get_next_query_stream().get_all();
            seqs.push_back(S);
        }

        k = 30;
        logger << "Building E. coli in memory..." << endl;
        NodeBOSSInMemoryConstructor<plain_matrix_sbwt_t> builder;
        builder.build(seqs, matrixboss_reference, k, true);

        logger << "Building E. coli with external memory..." << endl;
        plain_matrix_sbwt_t::BuildConfig config;
        config.input_files = {filename};
        config.k = k;
        config.build_streaming_support = true;
        config.ram_gigas = 2;
        config.n_threads = 2;
        config.min_abundance = 1;
        matrixboss = plain_matrix_sbwt_t(config);
    }

    // static void SetUpTestCase(){
    //     // SetUpTestCase was renamed to SetUpTestSuite at some point.
    //     // This funtion is for legacy support.
    //     SetUpTestSuite();
    // }

};

matrixboss_t TEST_LARGE::matrixboss;
matrixboss_t TEST_LARGE::matrixboss_reference;
LL TEST_LARGE::k;
vector<string> TEST_LARGE::seqs;

TEST_F(TEST_LARGE, verify_suffix_group_starts){
    sdsl::bit_vector reference = mark_suffix_groups(matrixboss.get_subset_rank_structure().A_bits, matrixboss.get_subset_rank_structure().C_bits, matrixboss.get_subset_rank_structure().G_bits, matrixboss.get_subset_rank_structure().T_bits, k);
    ASSERT_EQ(reference, matrixboss.get_streaming_support());
}

TEST_F(TEST_LARGE, check_bit_vectors){
    logger << matrixboss_reference.get_subset_rank_structure().A_bits.size() << " " << matrixboss.get_subset_rank_structure().A_bits.size() << endl;

    ASSERT_EQ(matrixboss.get_subset_rank_structure().A_bits, matrixboss_reference.get_subset_rank_structure().A_bits);
    ASSERT_EQ(matrixboss.get_subset_rank_structure().C_bits, matrixboss_reference.get_subset_rank_structure().C_bits);
    ASSERT_EQ(matrixboss.get_subset_rank_structure().G_bits, matrixboss_reference.get_subset_rank_structure().G_bits);
    ASSERT_EQ(matrixboss.get_subset_rank_structure().T_bits, matrixboss_reference.get_subset_rank_structure().T_bits);
    ASSERT_EQ(matrixboss.get_streaming_support(), matrixboss_reference.get_streaming_support());
}

TEST_F(TEST_LARGE, streaming_queries){
    SeqIO::Unbuffered_Reader sr("example_data/queries.fastq");
    while(!sr.done()){
        string S = sr.get_next_query_stream().get_all();
        vector<int64_t> result = matrixboss.streaming_search(S);
        for(LL i = 0; i < (LL)S.size()-k+1; i++){
            LL x = matrixboss.search(S.c_str() + i);
            ASSERT_EQ(result[i], x);
        }
    }
}

TEST_F(TEST_LARGE, query_lots_of_kmers){
    unordered_set<kmer_t> all_kmers; // Also collect a set of all k-mers in the input for later
    LL search_count = 0;

    logger << "Querying all input k-mers..." << endl;
    for(const string& S : seqs){
        for(LL i = 0; i < (LL)S.size() - k + 1; i++){
            string kmer = S.substr(i,k);
            bool is_valid = true;
            for(char c : kmer) if(c != 'A' && c != 'C' && c != 'G' && c != 'T') is_valid = false;
            if(is_valid){
                LL colex = matrixboss.search(kmer);
                ASSERT_GE(colex, 0); // Should be sound
                all_kmers.insert(kmer);
                search_count++;

                // Print verbose output
                if(search_count % 100000 == 0)  logger << kmer << " " << colex << endl;

            }
        }
    }
    logger << "Querying random k-mers that are not in the input..." << endl;
    srand(12514);
    for(LL rep = 0; rep < 1e6; rep++){
        string S = generate_random_kmer(k);
        if(all_kmers.count(S) == 0){
            LL colex = matrixboss.search(S);
            ASSERT_EQ(colex, -1); // Should not be found

            // Print verbose output
            if(rep % 100000 == 0)
                logger << S << " " << colex << endl;
        }
    }
}
