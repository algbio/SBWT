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

void print_basic_stats(const plain_matrix_sbwt_t& sbwt){
    cout << "k = " << sbwt.get_k() << endl;
    cout << "Number of kmers: " << sbwt.number_of_kmers() << endl;
    cout << "Number of columns: " << sbwt.number_of_subsets() << endl;
}

void print_column_distribution(const plain_matrix_sbwt_t& sbwt){
    map<vector<bool>, LL> counts; // Column bits, count of that column
    vector<bool> col(4);
    for(LL i = 0; i < sbwt.number_of_subsets(); i++){
        col[0] = sbwt.get_subset_rank_structure().A_bits[i];
        col[1] = sbwt.get_subset_rank_structure().C_bits[i];
        col[2] = sbwt.get_subset_rank_structure().G_bits[i];
        col[3] = sbwt.get_subset_rank_structure().T_bits[i];
        counts[col]++;
    }
    vector<pair<LL, vector<bool>>> counts_vec;
    for(auto [column, count] : counts)
        counts_vec.push_back({count, column});
    sort(counts_vec.begin(), counts_vec.end());

    cout << "Column distribution: " << endl;
    for(auto [count, column] : counts_vec){
        cout << column[0] << column[1] << column[2] << column[3] << " " << count << endl;
    }
}

int stats_main(int argc, char** argv){

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Compute statistics on the index (only plain-matrix variant supported)");

    options.add_options()
        ("i,index-file", "Index input file.", cxxopts::value<string>())
        ("d,column-distribution", "Also print column distribution.", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string indexfile = opts["index-file"].as<string>();
    bool do_column_distribution = opts["column-distribution"].as<bool>();
    check_readable(indexfile);

    vector<string> variants = get_available_variants();
    throwing_ifstream in(indexfile, ios::binary);
    string variant = load_string(in.stream); // read variant type

    if(variant != "plain-matrix"){
        throw std::runtime_error("Error: index is not in plain-matrix format. Only plain-matrix is supported at the moment");
    }

    plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);

    print_basic_stats(sbwt);
    if(do_column_distribution) print_column_distribution(sbwt);
    
    return 0;

}
