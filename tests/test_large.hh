#pragma once

#include "setup_tests.hh"
#include "globals.hh"
#include "Kmer.hh"
#include "SBWT.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include "SubsetWT.hh"
#include "SeqIO/SeqIO.hh"
#include "suffix_group_optimization.hh"
#include "variants.hh"
#include <gtest/gtest.h>
#include <unordered_set>

using namespace sbwt;


typedef Kmer<MAX_KMER_LENGTH> kmer_t;

typedef SBWT<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss_t;

class TEST_LARGE : public ::testing::Test {
    protected:

    static matrixboss_t matrixboss;
    static matrixboss_t matrixboss_reference;
    static int64_t k;
    static vector<string> seqs;
    static unordered_set<Kmer<MAX_KMER_LENGTH>> all_kmers;

    static void SetUpTestSuite(){
        string filename = "example_data/coli3.fna";

        k = 30;
        int64_t precalc_k = 5;

        seq_io::Reader sr(filename);
        string S;
        logger << "Reading sequences and hashing all k-mers" << endl;
        while(true){
            string S = sr.get_next_read();
            if(S == "") break;
            seqs.push_back(S);
            for(int64_t i = 0; i < (int64_t)S.size()-k+1; i++){
                string kmer = S.substr(i,k);
                bool is_valid = true;
                for(char c : kmer) if(c != 'A' && c != 'C' && c != 'G' && c != 'T') is_valid = false;
                if(is_valid) all_kmers.insert(kmer_t(kmer));
            }
        }

        
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

        // Add precalc
        matrixboss.do_kmer_prefix_precalc(precalc_k);
        matrixboss_reference.do_kmer_prefix_precalc(precalc_k);

    }

    // static void SetUpTestCase(){
    //     // SetUpTestCase was renamed to SetUpTestSuite at some point.
    //     // This funtion is for legacy support.
    //     SetUpTestSuite();
    // }

};

matrixboss_t TEST_LARGE::matrixboss;
matrixboss_t TEST_LARGE::matrixboss_reference;
int64_t TEST_LARGE::k;
vector<string> TEST_LARGE::seqs;
unordered_set<Kmer<MAX_KMER_LENGTH>> TEST_LARGE::all_kmers;

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
    seq_io::Reader sr("example_data/queries.fastq");
    while(true){
        string S = sr.get_next_read();
        if(S == "") break;
        vector<int64_t> result = matrixboss.streaming_search(S);
        for(int64_t i = 0; i < (int64_t)S.size()-k+1; i++){
            int64_t x = matrixboss.search(S.c_str() + i);
            ASSERT_EQ(result[i], x);
        }
    }
}

TEST_F(TEST_LARGE, dummy_node_marks){
    // This test is really weak but it's something
    sdsl::bit_vector marks = matrixboss.compute_dummy_node_marks();
    int64_t n_marks = 0;
    for(bool x : marks) n_marks += x;
    logger << "dummy node marks test: " << matrixboss.number_of_subsets() << " " << matrixboss.number_of_kmers() + n_marks << endl;
    ASSERT_EQ(matrixboss.number_of_subsets(), matrixboss.number_of_kmers() + n_marks);
}

TEST_F(TEST_LARGE, query_lots_of_kmers){
    int64_t search_count = 0;

    string ACGT = "ACGT";
    logger << "Querying all input k-mers..." << endl;
    for(const string& S : seqs){
        for(int64_t i = 0; i < (int64_t)S.size() - k + 1; i++){
            string kmer = S.substr(i,k);
            bool is_valid = true;
            for(char c : kmer) if(c != 'A' && c != 'C' && c != 'G' && c != 'T') is_valid = false;
            if(is_valid){
                int64_t colex = matrixboss.search(kmer);
                ASSERT_GE(colex, 0); // Should be found
                search_count++;

                // Try all forward moves
                for(char c : ACGT){
                    kmer_t next_kmer(kmer.substr(1) + c);
                    int64_t r = matrixboss.forward(colex, c);
                    if(all_kmers.count(next_kmer) == 0)
                        ASSERT_EQ(r, -1); // Should not be found
                    else
                        ASSERT_EQ(r, matrixboss.search(next_kmer.to_string()));
                }

                // Print verbose output
                if(search_count % 100000 == 0)  logger << kmer << " " << colex << endl;

            }
        }
    }
    logger << "Querying random k-mers that are not in the input..." << endl;
    srand(12514);
    for(int64_t rep = 0; rep < 1e6; rep++){
        string S = generate_random_kmer(k);
        if(all_kmers.count(S) == 0){
            int64_t colex = matrixboss.search(S);
            ASSERT_EQ(colex, -1); // Should not be found

            // Print verbose output
            if(rep % 100000 == 0)
                logger << S << " " << colex << endl;
        }
    }
}
