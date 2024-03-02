#include <iostream>
#include "throwing_streams.hh"
#include "globals.hh"
#include "variants.hh"
#include "SeqIO/SeqIO.hh"
#include "SubsetMatrixSelectSupport.hh"

using namespace std;
using namespace sbwt;

int main(int argc, char** argv){

    int64_t k = 6;

    // Build the index
    plain_matrix_sbwt_t::BuildConfig config;
    config.k = k;
    config.build_streaming_support = true; // One extra bit vector for speeding up positive streaming queries
    config.precalc_k = 4; // Speed up search by precalculating all 4^p p-mer intervals
    config.input_files = {"sequences.fna"};
    config.n_threads = 4;
    config.ram_gigas = 4;
    plain_matrix_sbwt_t sbwt(config);

    // Search for k-mer GATGGC
    cout << sbwt.search("GATGGC") << endl;

    // Search for all k-mers of TAATGCTGTAGC
    for(int64_t colex_rank : sbwt.streaming_search("TAATGCTGTAGC")){
        cout << colex_rank << endl;
    }

    // Dump all k-mers out of the data structure at once (fast)
    string kmer_dump = sbwt.reconstruct_all_kmers();
    for(int64_t i = 0; i < kmer_dump.size(); i += k){
        string kmer = kmer_dump.substr(i, k);

        // If the k-mer is not a dummy k-mer, print it
        if(kmer[0] != '$') cout << kmer << endl;
    }
    cout << "--" << endl;

    // List k-mers one by one
    SubsetMatrixSelectSupport<sdsl::bit_vector> select_support(sbwt.get_subset_rank_structure());

    vector<char> buf(k+1); // The k-mer will be written here
    for(int64_t i = 0; i < sbwt.number_of_subsets(); i++){
        sbwt.get_kmer_fast(i, buf.data(), select_support);

        // If the k-mer is not a dummy k-mer, print it
        if(buf[0] != '$') cout << buf.data() << endl;
    }

}