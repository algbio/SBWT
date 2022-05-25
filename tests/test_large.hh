#pragma once

#include "setup_tests.hh"
#include "globals.hh"
#include "Kmer.hh"
#include "NodeBOSS.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include "SubsetWT.hh"
#include "input_reading.hh"
#include "libwheeler/BOSS_builder.hh"
#include "suffix_group_optimization.hh"
#include "variants.hh"
#include <gtest/gtest.h>
#include <unordered_set>

typedef long long LL;
typedef Kmer<MAX_KMER_LENGTH> kmer_t;

typedef NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss_t;

class TEST_LARGE : public ::testing::Test {
    protected:

    matrixboss_t matrixboss;
    matrixboss_t matrixboss_reference;
    LL k;
    vector<string> seqs;

    void reverse_seqs_in_fasta(std::string infile, std::string outfile){
        Sequence_Reader sr(infile, FASTA_MODE);
        throwing_ofstream out(outfile);
        while(!sr.done()){
            Read_stream rs = sr.get_next_query_stream();
            string seq = rs.get_all();
            std::reverse(seq.begin(), seq.end());
            out << ">" + rs.header << "\n" << seq << "\n";
        }
    }

    void SetUp() override {
        string filename = "example_data/coli3.fna";

        string rev_file = get_temp_file_manager().create_filename("",".fna");
        reverse_seqs_in_fasta(filename, rev_file);

        Sequence_Reader sr(filename);
        while(!sr.done())
            seqs.push_back(sr.get_next_query_stream().get_all());

        k = 30;
        logger << "Building E. coli in memory..." << endl;
        matrixboss_reference.build_from_strings(seqs, k, true);

        logger << "Building E. coli with KMC..." << endl;
        matrixboss.build_using_KMC({rev_file}, k, true, 2, 2, 1);
    }

};

TEST_F(TEST_LARGE, check_matrix_bits){
    logger << matrixboss_reference.subset_rank.A_bits.size() << " " << matrixboss.subset_rank.A_bits.size() << endl;

    ASSERT_EQ(matrixboss.subset_rank.A_bits, matrixboss_reference.subset_rank.A_bits);
    ASSERT_EQ(matrixboss.subset_rank.C_bits, matrixboss_reference.subset_rank.C_bits);
    ASSERT_EQ(matrixboss.subset_rank.G_bits, matrixboss_reference.subset_rank.G_bits);
    ASSERT_EQ(matrixboss.subset_rank.T_bits, matrixboss_reference.subset_rank.T_bits);
}

TEST_F(TEST_LARGE, streaming_queries){
    Sequence_Reader sr("example_data/queries.fastq");
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