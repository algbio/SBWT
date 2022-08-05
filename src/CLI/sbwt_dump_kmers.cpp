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

char incoming_label(const plain_matrix_sbwt_t& sbwt, int64_t node){
    if(node < sbwt.get_C_array()[0]) return '$';
    else if(node < sbwt.get_C_array()[1]) return 'A';
    else if(node < sbwt.get_C_array()[2]) return 'C';
    else if(node < sbwt.get_C_array()[3]) return 'G';
    else return 'T';
}


int dump_main(int argc, char** argv){

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "Dump all k-mers out of the data structure (only plain matrix supported).");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,index-file", "Index input file.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string outfile = opts["out-file"].as<string>();
    string indexfile = opts["index-file"].as<string>();

    check_writable(outfile);
    check_readable(indexfile);

    write_log("Loading the DBG.", LogLevel::MAJOR);
    plain_matrix_sbwt_t sbwt;
    throwing_ifstream in(indexfile, ios::binary);
    if(load_string(in.stream) != "plain-matrix"){
        cerr << "Erro: only plain-matrix variant is supported" << endl;
        return 1;
    }    
    sbwt.load(in.stream);

    write_log("Building select supports", LogLevel::MAJOR);
    vector<sdsl::select_support_mcl<>> select_supports(4);
    sdsl::util::init_support(select_supports[0], &sbwt.get_subset_rank_structure().A_bits);
    sdsl::util::init_support(select_supports[1], &sbwt.get_subset_rank_structure().C_bits);
    sdsl::util::init_support(select_supports[2], &sbwt.get_subset_rank_structure().G_bits);
    sdsl::util::init_support(select_supports[3], &sbwt.get_subset_rank_structure().T_bits);

    write_log("Dumping k-mers to " + outfile, LogLevel::MAJOR);
    throwing_ofstream out(outfile);
    LL k = sbwt.get_k();
    for(LL i = 0; i < sbwt.number_of_subsets(); i++){
        vector<char> label(k);
        LL node = i;
        for(LL j = 0; j < k; j++){
            char c = incoming_label(sbwt, node);
            label[k - 1 - j] = c;
            if(c == '$') 
                continue;
            if(c == 'A')
                node = select_supports[0].select(node - sbwt.get_C_array()[0] + 1);
            if(c == 'C')
                node = select_supports[1].select(node - sbwt.get_C_array()[1] + 1);
            if(c == 'G')
                node = select_supports[2].select(node - sbwt.get_C_array()[2] + 1);
            if(c == 'T')
                node = select_supports[3].select(node - sbwt.get_C_array()[3] + 1);
        }

        if(label[0] != '$'){
            // Write out the label
            for(char c : label) out.stream << c;
            out.stream << "\n";
        }
    }

}
