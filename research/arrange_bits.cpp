#include "libwheeler/BOSS.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "globals.hh"
#include "NodeBOSS.hh"
#include "SubsetWT.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include "suffix_group_optimization.hh"
#include <filesystem>

//#include "MEF.hpp"

using namespace sdsl;

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

    cxxopts::Options options(argv[0], "This program prints a NodeBOSS representation of the index. Not production quality code.");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,input-matrixboss", "The input matrixBOSS", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("s,strategy", "Bit shifting strategy to use. Possible options: \"pushleft\" or \"spread\"", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage examples: " << argv[0] << " -i in.matrixboss -o out.matrixboss --temp-dir temp --strategy pushleft" << endl;
        exit(1);
    }

    string outfile = opts["out-file"].as<string>();
    string infile = opts["input-matrixboss"].as<string>();
    string temp_dir = opts["temp-dir"].as<string>();
    string strategy = opts["strategy"].as<string>();

    if(strategy != "pushleft" && strategy != "spread"){
        cerr << "Invalid strategy: " << strategy << endl;
        return 1;
    }

    check_writable(outfile);
    check_readable(infile);
    std::filesystem::create_directory(temp_dir);

    write_log("Loading the matrixBOSS", LogLevel::MAJOR);
    SBWT<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss;
    throwing_ifstream in(infile, ios::binary);
    matrixboss.load(in.stream);

    LL n_nodes = matrixboss.n_nodes;
    LL k = matrixboss.k;

    // Mark suffix group starts
    write_log("Computing suffix group boundaries", LogLevel::MAJOR);
    sdsl::bit_vector suffix_group_starts = mark_suffix_groups(matrixboss.subset_rank.A_bits,
                                                              matrixboss.subset_rank.C_bits,
                                                              matrixboss.subset_rank.G_bits,
                                                              matrixboss.subset_rank.T_bits,
                                                              matrixboss.C,
                                                              matrixboss.k);

    sdsl::bit_vector A_copy = matrixboss.subset_rank.A_bits;
    sdsl::bit_vector C_copy = matrixboss.subset_rank.C_bits;
    sdsl::bit_vector G_copy = matrixboss.subset_rank.G_bits;
    sdsl::bit_vector T_copy = matrixboss.subset_rank.T_bits;

    write_log("Arranging bits", LogLevel::MAJOR);
    if(strategy == "pushleft"){
        push_bits_left(A_copy, C_copy, G_copy, T_copy, suffix_group_starts);
    }

    if(strategy == "spread"){
        push_bits_left(A_copy, C_copy, G_copy, T_copy, suffix_group_starts);
        spread_bits_after_push_left(A_copy, C_copy, G_copy, T_copy, suffix_group_starts);
    }


    write_log("Building new MatrixBOSS", LogLevel::MAJOR);
    SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>> new_matrixrank(A_copy, C_copy, G_copy, T_copy);
    SBWT<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> new_matrixboss;
    new_matrixboss.build_from_bit_matrix(A_copy, C_copy, G_copy, T_copy, matrixboss.k);

    throwing_ofstream out(outfile);
    LL n_bytes = new_matrixboss.serialize(out.stream);
    cerr << n_bytes << " bytes written to " << outfile << endl;



}

