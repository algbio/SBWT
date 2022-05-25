#include "libwheeler/BOSS.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "globals.hh"
#include "NodeBOSS.hh"
#include "SubsetWT.hh"
#include "stdlib_printing.hh"
#include "SeqIO.hh"
#include "SubsetMatrixRank.hh"
#include <filesystem>

//#include "MEF.hpp"

using namespace std;
typedef long long LL;

void check_readable(string filename){
    throwing_ifstream F(filename); // Throws on failure
}

// Also clears the file
void check_writable(string filename){
    throwing_ofstream F(filename, std::ofstream::out | std::ofstream::app); // Throws on failure
}

int main(int argc, char** argv){

    set_log_level(LogLevel::MINOR);

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "Query all k-mers of all input reads. Assumes all reads only contain characters A,C,G and T.");

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
    NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss;
    throwing_ifstream in(indexfile, ios::binary);
    matrixboss.load(in.stream);

    write_log("Building select supports", LogLevel::MAJOR);
    vector<sdsl::select_support_mcl<>> select_supports(4);
    sdsl::util::init_support(select_supports[0], &matrixboss.subset_rank.A_bits);
    sdsl::util::init_support(select_supports[1], &matrixboss.subset_rank.C_bits);
    sdsl::util::init_support(select_supports[2], &matrixboss.subset_rank.G_bits);
    sdsl::util::init_support(select_supports[3], &matrixboss.subset_rank.T_bits);

    write_log("Dumping k-mers to " + outfile, LogLevel::MAJOR);
    throwing_ofstream out(outfile);
    LL k = matrixboss.k;
    for(LL i = 0; i < matrixboss.n_nodes; i++){
        vector<char> label(k);
        LL node = i;
        for(LL j = 0; j < k; j++){
            char c = matrixboss.incoming_label(node);
            label[k - 1 - j] = c;
            if(c == '$') 
                continue;
            if(c == 'A')
                node = select_supports[0].select(node - matrixboss.C[0] + 1);
            if(c == 'C')
                node = select_supports[1].select(node - matrixboss.C[1] + 1);
            if(c == 'G')
                node = select_supports[2].select(node - matrixboss.C[2] + 1);
            if(c == 'T')
                node = select_supports[3].select(node - matrixboss.C[3] + 1);
        }

        // Write out the label
        for(char c : label) out.stream << c;
        out.stream << "\n";
        
    }

}

