#include "libwheeler/BOSS.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "globals.hh"
#include "NodeBOSS.hh"
#include "SubsetWT.hh"
#include "stdlib_printing.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include <filesystem>
#include "MEF.hpp"

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
        ("o,out-file", "Output filename prefix.", cxxopts::value<string>())
        ("i,matrixboss-file", "The .matrixboss file of an index.", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("variant", "MatrixBOSS variant to build. Possible values are: plain-matrix, rrr-matrix, mef-matrix, plain-split, rrr-split, mef-split, plain-concat, vari-concat, plain-subset and rrr-subset", cxxopts::value<string>())
        ("h,help", "Print usage")
        ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage examples: " << argv[0] << " -i in_prefix -o out_prefix --temp-dir temp" << endl;
        exit(1);
    }

    string out_prefix = opts["out-file"].as<string>();
    string index_dbg_file = opts["matrixboss-file"].as<string>();
    string temp_dir = opts["temp-dir"].as<string>();
    string variant_opt = opts["variant"].as<string>();

    enum boss_variant {
        PLAIN_MATRIX, RRR_MATRIX, MEF_MATRIX,
        PLAIN_SPLIT, RRR_SPLIT, MEF_SPLIT,
        PLAIN_CONCAT, VARI_CONCAT,
        PLAIN_SUBSET, RRR_SUBSET
        //, MEF_SUBSET
    };

    boss_variant selected_var = PLAIN_MATRIX;

    if (variant_opt == "plain-matrix")
        selected_var = PLAIN_MATRIX;
    if (variant_opt == "rrr-matrix")
        selected_var = RRR_MATRIX;
    if (variant_opt == "mef-matrix")
        selected_var = MEF_MATRIX;
    if (variant_opt == "plain-split")
        selected_var = PLAIN_SPLIT;
    if (variant_opt == "rrr-split")
        selected_var = RRR_SPLIT;
    if (variant_opt == "mef-split")
        selected_var = MEF_SPLIT;
    if (variant_opt == "plain-concat")
        selected_var = PLAIN_CONCAT;
    if (variant_opt == "vari-concat")
        selected_var = VARI_CONCAT;
    if (variant_opt == "plain-subset")
        selected_var = PLAIN_SUBSET;
    if (variant_opt == "rrr-subset")
        selected_var = RRR_SUBSET;
    // if (variant_opt == "mef-subset")
    //     selected_var = MEF_SUBSET;

    check_readable(index_dbg_file);
    std::filesystem::create_directory(temp_dir);

    write_log("Loading the DBG", LogLevel::MAJOR);

    throwing_ifstream in(index_dbg_file, ios::binary);

    NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss;
    write_log("Building MatrixBOSS representation", LogLevel::MAJOR);
    matrixboss.load(in.stream);
    sdsl::bit_vector A_bits = matrixboss.subset_rank.A_bits;
    sdsl::bit_vector C_bits = matrixboss.subset_rank.C_bits;
    sdsl::bit_vector G_bits = matrixboss.subset_rank.G_bits;
    sdsl::bit_vector T_bits = matrixboss.subset_rank.T_bits;
    const int64_t k = matrixboss.k;
    const int64_t n_nodes = matrixboss.n_nodes;

    // matrices
    NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> plain_matrixboss;
    NodeBOSS<SubsetMatrixRank<sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type>> rrr_matrixboss;
    NodeBOSS<SubsetMatrixRank<mod_ef_vector<>, mod_ef_vector<>::rank_1_type>> mef_matrixboss;

    // splits
    NodeBOSS<SubsetSplitRank<sdsl::bit_vector, sdsl::rank_support_v5<>,
                             sdsl::bit_vector, sdsl::rank_support_v5<>>> plain_splitboss;

    NodeBOSS<SubsetSplitRank<sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type,
                             sdsl::bit_vector, sdsl::rank_support_v5<>>> rrr_splitboss;


    NodeBOSS<SubsetSplitRank<mod_ef_vector<>, mod_ef_vector<>::rank_1_type,
                             sdsl::bit_vector, sdsl::rank_support_v5<>>> mef_splitboss;

    // concats
    NodeBOSS<SubsetConcatRank<sdsl::bit_vector,
                              sdsl::bit_vector::select_0_type,
                              sdsl::wt_blcd<sdsl::bit_vector,
                                            sdsl::rank_support_v5<>,
                                            sdsl::select_support_scan<1>,
                                            sdsl::select_support_scan<0>>>
             > plain_concatboss;

    NodeBOSS<SubsetConcatRank<sd_vector<>,
                              sd_vector<>::select_0_type,
                              sdsl::wt_blcd<rrr_vector<63>,
                                            rrr_vector<>::rank_1_type,
                                            rrr_vector<>::select_1_type,
                                            rrr_vector<>::select_0_type>>
             > vari_concatboss;

    // subsets
    NodeBOSS<SubsetWT<sdsl::wt_blcd<sdsl::bit_vector,
                                    sdsl::rank_support_v5<>,
                                    sdsl::select_support_scan<1>,
                                    sdsl::select_support_scan<0>>>
             > plain_sswtboss;


    NodeBOSS<SubsetWT<sdsl::wt_blcd<sdsl::rrr_vector<>,
                                    sdsl::rrr_vector<>::rank_1_type,
                                    rrr_vector<>::select_1_type,
                                    rrr_vector<>::select_0_type>>
             > rrr_sswtboss;

    // NodeBOSS<SubsetWT<sdsl::wt_blcd<mod_ef_vector<>,
    //                                 mod_ef_vector<>::rank_1_type,
    //                                 sdsl::select_support_scan<1>,
    //                                 sdsl::select_support_scan<0>>>
    //          > mef_sswtboss;

    throwing_ofstream boss_out(out_prefix, ios::binary);

    switch (selected_var) {
    case PLAIN_MATRIX:
        write_log("Writing plain MatrixBOSS to disk", LogLevel::MAJOR);
        write_log("plain MatrixBOSS: " +
                  to_string(matrixboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
                  " bits per k-mer", LogLevel::MAJOR);
        break;
    case RRR_MATRIX:
        write_log("Building rrr MatrixBOSS representation", LogLevel::MAJOR);
        rrr_matrixboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        write_log("Writing rrr matrixBOSS to disk", LogLevel::MAJOR);
        write_log("rrr MatrixBOSS: " +
                  to_string(rrr_matrixboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
                  " bits per k-mer", LogLevel::MAJOR);
        break;

    case MEF_MATRIX:
        write_log("Building mef MatrixBOSS representation", LogLevel::MAJOR);
        mef_matrixboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        write_log("Writing mef matrixBOSS to disk", LogLevel::MAJOR);
        write_log("mef MatrixBOSS: " +
                  to_string(mef_matrixboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
                  " bits per k-mer", LogLevel::MAJOR);
        break;
    case PLAIN_SPLIT:
        write_log("Writing plain SplitBOSS to disk", LogLevel::MAJOR);
        plain_splitboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        write_log("plain SplitBOSS: " +
                  to_string(plain_splitboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
                  " bits per k-mer", LogLevel::MAJOR);
        break;
    case RRR_SPLIT:
        write_log("Building rrr SplitBOSS representation", LogLevel::MAJOR);
        rrr_splitboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        write_log("Writing rrr splitBOSS to disk", LogLevel::MAJOR);
        write_log("rrr SplitBOSS: " +
                  to_string(rrr_splitboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
                  " bits per k-mer", LogLevel::MAJOR);
        break;
    case MEF_SPLIT:
        write_log("Building mef SplitBOSS representation", LogLevel::MAJOR);
        mef_splitboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        write_log("Writing mef splitBOSS to disk", LogLevel::MAJOR);
        write_log("mef SplitBOSS: " +
                  to_string(mef_splitboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
                  " bits per k-mer", LogLevel::MAJOR);
        break;
    case PLAIN_CONCAT:
        write_log("Writing plain ConcatBOSS to disk", LogLevel::MAJOR);
        plain_concatboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        write_log("plain ConcatBOSS: " +
                  to_string(plain_concatboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
                  " bits per k-mer", LogLevel::MAJOR);
        break;
    case VARI_CONCAT:
        write_log("Building vari ConcatBOSS representation", LogLevel::MAJOR);
        vari_concatboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        write_log("Writing vari ConcatBOSS to disk", LogLevel::MAJOR);
        write_log("vari ConcatBOSS: " +
                  to_string(vari_concatboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
                  " bits per k-mer", LogLevel::MAJOR);
        break;
    case PLAIN_SUBSET:
        write_log("Writing plain SswtBOSS to disk", LogLevel::MAJOR);
        plain_sswtboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        write_log("plain SswtBOSS: " +
                  to_string(plain_sswtboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
                  " bits per k-mer", LogLevel::MAJOR);
        break;

    case RRR_SUBSET:
        write_log("Building rrr SswtBOSS representation", LogLevel::MAJOR);
        rrr_sswtboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        write_log("Writing rrr sswtBOSS to disk", LogLevel::MAJOR);
        write_log("rrr SswtBOSS: " +
                  to_string(rrr_sswtboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
                  " bits per k-mer", LogLevel::MAJOR);
        break;
        // case MEF_SUBSET:
        //     write_log("Building mef SswtBOSS representation", LogLevel::MAJOR);
        //     mef_sswtboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        //     write_log("Writing mef sswtBOSS to disk", LogLevel::MAJOR);
        //     write_log("mef SswtBOSS: " +
        //               to_string(mef_sswtboss.serialize(boss_out.stream) * 8.0 / n_nodes) +
        //               " bits per k-mer", LogLevel::MAJOR);
        //     break;
    default:
        break;
    }
}
