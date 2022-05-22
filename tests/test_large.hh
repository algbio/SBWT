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
#include <gtest/gtest.h>
#include <unordered_set>

typedef long long LL;
typedef Kmer<MAX_KMER_LENGTH> kmer_t;

typedef NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss_t;

TEST(TEST_LARGE, e_coli){
    Sequence_Reader sr("example_data/coli3.fna", FASTA_MODE);
    vector<string> seqs;
    while(!sr.done())
        seqs.push_back(sr.get_next_query_stream().get_all());
    matrixboss_t matrixboss;
    logger << "Building E. coli..." << endl;
    LL k = 30;
    matrixboss.build_from_strings(seqs, k);
    
    logger << "Querying all k-mers in the input..." << endl;
    unordered_set<kmer_t> all_kmers; // Also collect a set of all k-mers in the input
    LL search_count = 0;
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