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
#include <chrono>

using namespace sdsl;

using namespace std;
typedef long long LL;

// Returns pair (mean, stddev) in microseconds pre query
template<typename boss_t>
pair<double, double> time_queries(const boss_t& boss, const vector<string>& queries, LL repeats){

    vector<double> trials;

    for(LL repeat = 0; repeat < repeats; repeat++){
        LL total_micros = 0; // microseconds

        for(const string& kmer : queries){
            // Search nodeboss
            LL t0 = cur_time_micros();
            boss.search(kmer);
            total_micros += cur_time_micros() - t0;
        }
        trials.push_back((double) total_micros / queries.size());
    }

    double mean = 0;
    for(double x : trials) mean += x / trials.size();

    double variance = 0;
    for(double x : trials) variance += (mean - x) * (mean - x) / trials.size();

    return {mean, sqrt(variance)};
}

void check_readable(string filename){
    throwing_ifstream F(filename); // Throws on failure
}

// Also clears the file
void check_writable(string filename){
    throwing_ofstream F(filename, std::ofstream::out | std::ofstream::app); // Throws on failure
}

int main(int argc, char** argv){

    set_log_level(LogLevel::OFF);

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "This program prints a NodeBOSS representation of the index. Not production quality code.");

    options.add_options()
        ("i,index-file", "The index file.", cxxopts::value<string>())
        ("q,query-file", "The query file to get kmers from.", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("variant", "MatrixBOSS variant in the index file. The value should be the same as the value that was used to build the variant.", cxxopts::value<string>())
        ("repeats", "Number of times to run the query set", cxxopts::value<LL>()->default_value("1"))
        ("h,help", "Print usage")
        ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage examples: " << argv[0] << " -i in_prefix -o out_prefix --temp-dir temp" << endl;
        exit(1);
    }

    string index_dbg_file = opts["index-file"].as<string>();
    string temp_dir = opts["temp-dir"].as<string>();
    string variant_opt = opts["variant"].as<string>();
    string query_file = opts["query-file"].as<string>();
    LL repeats = opts["repeats"].as<LL>();


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

    // write_log("Loading the DBG", LogLevel::MAJOR);

    throwing_ifstream in(index_dbg_file, ios::binary);

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

    plain_matrixboss.load(in.stream);
    sdsl::bit_vector A_bits = plain_matrixboss.subset_rank.A_bits;
    sdsl::bit_vector C_bits = plain_matrixboss.subset_rank.C_bits;
    sdsl::bit_vector G_bits = plain_matrixboss.subset_rank.G_bits;
    sdsl::bit_vector T_bits = plain_matrixboss.subset_rank.T_bits;
    const int64_t k = plain_matrixboss.k;
    const int64_t n_nodes = plain_matrixboss.n_nodes;

    // write_log("Loading the queries", LogLevel::MAJOR);

    // ifstream qin(query_file, ios::binary);

    vector<string> queries;
    string s;
    ifstream qin(query_file);
    while (std::getline(qin, s)) {
        queries.push_back(s);
    }

    qin.close();

    double mean = 0;
    double stddev = 0;
    switch (selected_var) {
    case PLAIN_MATRIX:
        // plain_matrixboss.load(in.stream);
        std::tie(mean, stddev) = time_queries(plain_matrixboss, queries, repeats);
        cout << mean << " " << stddev << endl;
        break;
    case RRR_MATRIX:
        // rrr_matrixboss.load(in.stream);
        rrr_matrixboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        std::tie(mean, stddev) = time_queries(rrr_matrixboss, queries, repeats);
        cout << mean << " " << stddev << endl;
        break;
    case MEF_MATRIX:
        // mef_matrixboss.load(in.stream);
        mef_matrixboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        std::tie(mean, stddev) = time_queries(mef_matrixboss, queries, repeats);
        cout << mean << " " << stddev << endl;
        break;
    case PLAIN_SPLIT:
        // plain_splitboss.load(in.stream);
        plain_splitboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        std::tie(mean, stddev) = time_queries(plain_splitboss, queries, repeats);
        cout << mean << " " << stddev << endl;
        break;
    case RRR_SPLIT:
        // rrr_splitboss.load(in.stream);
        rrr_splitboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        std::tie(mean, stddev) = time_queries(rrr_splitboss, queries, repeats);
        cout << mean << " " << stddev << endl;
        break;
    case MEF_SPLIT:
        // mef_splitboss.load(in.stream);
        mef_splitboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        std::tie(mean, stddev) = time_queries(mef_splitboss, queries, repeats);
        cout << mean << " " << stddev << endl;
        break;
    case PLAIN_CONCAT:
        // plain_concatboss.load(in.stream);
        plain_concatboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        std::tie(mean, stddev) = time_queries(plain_concatboss, queries, repeats);
        cout << mean << " " << stddev << endl;
        break;
    case VARI_CONCAT:
        // vari_concatboss.load(in.stream);
        vari_concatboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        std::tie(mean, stddev) = time_queries(vari_concatboss, queries, repeats);
        cout << mean << " " << stddev << endl;
        break;
    case PLAIN_SUBSET:
        // plain_sswtboss.load(in.stream);
        plain_sswtboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        std::tie(mean, stddev) = time_queries(plain_sswtboss, queries, repeats);
        cout << mean << " " << stddev << endl;
        break;
    case RRR_SUBSET:
        // rrr_sswtboss.load(in.stream);
        rrr_sswtboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k);
        std::tie(mean, stddev) = time_queries(rrr_sswtboss, queries, repeats);
        cout << mean << " " << stddev << endl;
        break;
        // case MEF_SUBSET:
        //     break;
    default:
        break;
    }
}
