#include <string>
#include <cstring>
#include "cxxopts.hpp"
#include "globals.hh"
#include "SBWT.hh"
#include "SubsetWT.hh"
#include "stdlib_printing.hh"
#include "SeqIO.hh"
#include "SubsetMatrixRank.hh"
#include "buffered_streams.hh"
#include "variants.hh"
#include "commands.hh"
#include <filesystem>
#include <cstdio>

using namespace std;
typedef long long LL;
using namespace sbwt;

int stats_main(int argc, char** argv){

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "SBWT stats.");

    options.add_options()
        ("i,index-file", "Index input file. On works for plain matrix.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string indexfile = opts["index-file"].as<string>();
    check_readable(indexfile);

    throwing_ifstream in(indexfile, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: only works for plain-matrix" << endl;
        return 1;
    }

    plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);

    map<vector<bool>, int64_t> column_counts;

    vector<bool> col(4);
    for(int64_t i = 0; i < sbwt.number_of_subsets(); i++){
        col[0] = sbwt.get_subset_rank_structure().A_bits[i];
        col[1] = sbwt.get_subset_rank_structure().C_bits[i];
        col[2] = sbwt.get_subset_rank_structure().G_bits[i];
        col[3] = sbwt.get_subset_rank_structure().T_bits[i];
        column_counts[col]++;
    }

    vector<double> distribution;
    for(auto [v, count] : column_counts){
        for(bool b : v) cout << b;
        cout << " ";
        cout << count << endl;
        distribution.push_back((double)count / sbwt.number_of_subsets());
    }

    cout << "Number of columns: " << sbwt.number_of_subsets() << endl;
    cout << "Entropy: " << entropy(distribution) << endl;

    return 0;

}