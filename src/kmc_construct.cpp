#include "KMC/include/kmc_runner.h"
#include "KMC_code.hh"
#include "globals.hh"
#include "cxxopts.hpp"
#include "variants.hh"
#include "throwing_streams.hh"
#include "kmc_construct.hh"
#include "stdint.h"
#include <iostream>

typedef long long LL;

int main(int argc, char** argv){

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Construct the SBWT index");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,in-file", "Input filename in FASTA format.", cxxopts::value<string>())
        ("k,node-length", "The k of the k-mers.", cxxopts::value<LL>())
        ("t,n-threads", "Number of threads", cxxopts::value<LL>()->default_value("1"))
        ("m,ram-gigabytes", "Maximum number of RAm to use in gigabytes.", cxxopts::value<LL>()->default_value("4"))
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string outfile = opts["out-file"].as<string>();
    string infile = opts["in-file"].as<string>();
    string temp_dir = opts["temp-dir"].as<string>();
    LL k = opts["k"].as<LL>();
    LL n_threads = opts["n-threads"].as<LL>();
    LL ram_gigas = opts["ram-gigabytes"].as<LL>();

    check_writable(outfile);
    check_readable(infile);
    std::filesystem::create_directory(temp_dir);

    get_temp_file_manager().set_dir(temp_dir);

    NodeBOSSKMCConstructor<plain_matrix_sbwt_t> X;
    plain_matrix_sbwt_t index;
    X.build(infile, index, k, n_threads, ram_gigas, false);
    
}