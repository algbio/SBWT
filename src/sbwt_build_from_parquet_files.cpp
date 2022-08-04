#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "globals.hh"
#include "SBWT.hh"
#include "SubsetWT.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include "Parquet_Kmer_Stream.hh"
#include "throwing_streams.hh"
#include "variants.hh"
#include <filesystem>

//#include "MEF.hpp"

using namespace sdsl;

using namespace std;
using namespace sbwt;

typedef long long LL;

vector<string> split(string s, char delimiter){
    stringstream test(s);
    string segment;
    vector<string> seglist;

    while(getline(test, segment, delimiter))
    {
        seglist.push_back(segment);
    }
    return seglist;
}

// Format is like this:
// ATGAGGGCGCGCTCAATGGGCACGCCGTTT,1,AG
// ATTAGGGCGCGCTCAATGGGCACGCCGTTT,1,G
// ATGAGGGCGTGCTCAATGGGCACGCCGTTT,1,G

vector<string> read_lines(string filename){
    string line;
    throwing_ifstream in(filename);
    vector<string> lines;
    while(getline(in.stream, line)){
        lines.push_back(line);
    }
    return lines;
}

int main(int argc, char** argv){

    set_log_level(LogLevel::MINOR);

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "This program creates a plain matrix SBWT out of output of the Spark-Themisto code of Jaakko Vuohtoniemi.");

    options.add_options()
        ("i,in-directory", "Parquet file directory", cxxopts::value<string>())
        ("o,out-file", "Output file path.", cxxopts::value<string>())
        ("k", "The k-mer k.", cxxopts::value<LL>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string in_directory = opts["in-directory"].as<string>();
    string outfile = opts["out-file"].as<string>();
    LL k = opts["k"].as<LL>();

    check_writable(outfile);

    Parquet_Kmer_Stream stream(in_directory, k);

    string prev_kmer;
    LL start_of_suffix_group = 0;
    LL kmer_idx = 0;
    vector<vector<bool> > rows(4); // DNA alphabet
    bool root_kmer_exists = false;
    string kmer, kmer_outlabels;
    LL n_kmers_without_dummies = 0;
    while(stream.next(kmer, kmer_outlabels)){
        if(kmer.back() == '$') root_kmer_exists = true; // K-mer is all dollars if the last k-mer is a dollar
        if(std::find(kmer.begin(), kmer.end(), '$') == kmer.end()) n_kmers_without_dummies++;
        if(prev_kmer == "" || prev_kmer.substr(1) != kmer.substr(1)){
            start_of_suffix_group = kmer_idx;
        }

        for(LL i = 0; i < 4; i++) rows[i].push_back(0);
        for(char c : kmer_outlabels){
            if(c == 'A') rows[0][start_of_suffix_group] = 1;
            if(c == 'C') rows[1][start_of_suffix_group] = 1;
            if(c == 'G') rows[2][start_of_suffix_group] = 1;
            if(c == 'T') rows[3][start_of_suffix_group] = 1;
        }

        kmer_idx++;
        prev_kmer = kmer;
    
    }

    cerr << "Root k-mer exists: " << (root_kmer_exists ? "true" : "false") << endl;

    LL n = rows[0].size();
    vector<sdsl::bit_vector> sdsl_vectors;
    for(LL i = 0; i < 4; i++){
        sdsl::bit_vector vec;
        if(root_kmer_exists) vec.resize(n);
        else{
            vec.resize(n+1);
            vec[0] = 0;
        }
        for(LL j = 0; j < n; j++) vec[j + (!root_kmer_exists)] = rows[i][j];
        sdsl_vectors.push_back(vec);
    }

    sdsl::bit_vector empty;
    plain_matrix_sbwt_t matrixboss(sdsl_vectors[0], sdsl_vectors[1], sdsl_vectors[2], sdsl_vectors[3], empty, k, n_kmers_without_dummies);
    throwing_ofstream out(outfile, ios::binary);
    matrixboss.serialize(out.stream);

}