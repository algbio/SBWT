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
#include "variants.hh"
#include <filesystem>

//#include "MEF.hpp"

using namespace std;
typedef long long LL;

template<typename sbwt_t>
void run_queries(const string& queryfile, const string& outfile, const sbwt_t& sbwt){
    write_log("Running queries", LogLevel::MAJOR);
    throwing_ofstream out(outfile);
    Sequence_Reader_Buffered sr(queryfile, FASTA_MODE);
    LL k = sbwt.k;
    while(true){ 
        LL len = sr.get_next_read_to_buffer();
        if(len == 0) break;
        for(LL i = 0; i < len - k + 1; i++){
            LL v = sbwt.search(sr.read_buf + i, k);
            out.stream << v << " ";
        }
        out << "\n";
    }
}

int search_main(int argc, char** argv){

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Query all k-mers of all input reads.");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,index-file", "Index input file.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA format.", cxxopts::value<string>())
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
    string queryfile = opts["query-file"].as<string>();

    check_writable(outfile);
    check_readable(indexfile);
    check_readable(queryfile);

    vector<string> variants = get_available_variants();

    throwing_ifstream in(indexfile, ios::binary);
    string variant = load_string(in.stream);
    if(std::find(variants.begin(), variants.end(), variant) == variants.end()){
        cerr << "Error loading index from file: unrecognized variant specified in the file" << endl;
        return 1;
    }

    write_log("Loading the index variant " + variant, LogLevel::MAJOR);

    if (variant == "plain-matrix"){
        plain_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        run_queries(queryfile, outfile, sbwt);
    }
    if (variant == "rrr-matrix"){
        rrr_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        run_queries(queryfile, outfile, sbwt);
    }
    if (variant == "mef-matrix"){
        mef_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        run_queries(queryfile, outfile, sbwt);
    }
    if (variant == "plain-split"){
        plain_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        run_queries(queryfile, outfile, sbwt);
    }
    if (variant == "rrr-split"){
        rrr_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        run_queries(queryfile, outfile, sbwt);
    }
    if (variant == "mef-split"){
        mef_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        run_queries(queryfile, outfile, sbwt);
    }
    if (variant == "plain-concat"){
        plain_concat_sbwt_t sbwt;
        sbwt.load(in.stream);
        run_queries(queryfile, outfile, sbwt);
    }
    if (variant == "mef-concat"){
        mef_concat_sbwt_t sbwt;
        sbwt.load(in.stream);
        run_queries(queryfile, outfile, sbwt);
    }
    if (variant == "plain-subsetwt"){
        plain_sswt_sbwt_t sbwt;
        sbwt.load(in.stream);
        run_queries(queryfile, outfile, sbwt);
    }
    if (variant == "rrr-subsetwt"){
        rrr_sswt_sbwt_t sbwt;
        sbwt.load(in.stream);
        run_queries(queryfile, outfile, sbwt);
    }

    return 0;

}

