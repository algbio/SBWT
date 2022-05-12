#include "libwheeler/BOSS.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "globals.hh"
#include "NodeBOSS.hh"
#include "SubsetWT.hh"
#include "stdlib_printing.hh"
#include "input_reading.hh"
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

template<typename otherboss_t>
void run_test(const BOSS<sdsl::bit_vector>& boss, const otherboss_t& otherboss){

    LL otherboss_total_millis = 0;
    LL wheleerboss_total_millis = 0;
    LL n_queries = 0;
    // Search for all nodes
    for(LL v = 0; v < boss.number_of_nodes(); v++){
        string kmer = boss.get_node_label(v);
        if(kmer.size() < boss.get_k()) continue;
        else n_queries++;

        // Search nodeboss
        LL t0_otherboss = cur_time_millis();
        LL search_result = otherboss.search(kmer);
        otherboss_total_millis += cur_time_millis() - t0_otherboss;
        
        LL t0_wheelerboss = cur_time_millis();
        LL search_result2 = boss.find_kmer(kmer);
        wheleerboss_total_millis += cur_time_millis() - t0_wheelerboss;

        if(search_result != search_result2){
            cout << "WRONG ANSWER " << v << " " << kmer << endl;
            exit(0);
        }

        if(n_queries == 1000000) break;
    }

    cout << "WheleerBOSS ms/query: " << wheleerboss_total_millis / (double)n_queries << endl;
    cout << "OtherBOSS ms/query: " << otherboss_total_millis / (double)n_queries << endl;
    cout << "Speedup: " << (double)wheleerboss_total_millis / otherboss_total_millis << endl;

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
        ("q,query-file", "The query in FASTA format.", cxxopts::value<string>())
        ("k,node-length", "The k of the k-mers.", cxxopts::value<LL>())
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
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
    string temp_dir = opts["temp-dir"].as<string>();
    string queryfile = opts["query-file"].as<string>();
    LL k = opts["k"].as<LL>();

    check_writable(outfile);
    check_readable(indexfile);
    check_readable(queryfile);
    std::filesystem::create_directory(temp_dir);

    write_log("Loading the DBG", LogLevel::MAJOR);
    NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss;
    throwing_ifstream in(indexfile, ios::binary);
    matrixboss.load(in.stream);

    throwing_ofstream out(outfile);
    Sequence_Reader sr(queryfile, FASTA_MODE);
    while(!sr.done()){
        string read = sr.get_next_query_stream().get_all();
        for(LL i = 0; i < (LL)read.size() - k + 1; i++){
            LL v = matrixboss.search(read.c_str() + i, k);
            out.stream << v << " ";
        }
        out << "\n";
    }

}

