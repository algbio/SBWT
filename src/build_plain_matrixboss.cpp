
#include "globals.hh"
#include "throwing_streams.hh"
#include "cxxopts.hpp"
#include "BOSS.hh"
#include "NodeBOSS.hh"
#include "SubsetMatrixRank.hh"

typedef long long LL;
using namespace std;

void check_readable(string filename){
    throwing_ifstream F(filename); // Throws on failure
}

// Also clears the file
void check_writable(string filename){
    throwing_ofstream F(filename, std::ofstream::out | std::ofstream::app); // Throws on failure
}


int main(int argc, char** argv){

    set_log_level(LogLevel::MAJOR);

    cxxopts::Options options(argv[0], "This constructs a matrixboss out of a .tdbg file. Not production quality code.");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,dbg-file", "The .tdbg file of an index.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage examples: " << argv[0] << " -i in.tdbg -o out.matrixboss" << endl;
        exit(1);
    }

    string out_file = opts["out-file"].as<string>();
    string tdbg_file = opts["dbg-file"].as<string>();

    check_writable(out_file);
    check_readable(tdbg_file);

    write_log("Loading the DBG", LogLevel::MAJOR);
    BOSS<sdsl::bit_vector> wheelerBOSS;
    throwing_ifstream in(tdbg_file, ios::binary);
    wheelerBOSS.load(in.stream);
    LL n_nodes = wheelerBOSS.number_of_nodes();

    write_log("WheelerBOSS has order k = " + to_string(wheelerBOSS.get_k()) + " (node label length)", LogLevel::MAJOR);
    
    NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss_plain;
    write_log("Building MatrixBOSS representation", LogLevel::MAJOR);
    matrixboss_plain.build_from_WheelerBOSS(wheelerBOSS);

    write_log("Writing plain matrixBOSS to disk", LogLevel::MAJOR);
    throwing_ofstream matrixboss_out(out_file, ios::binary);
    write_log("MatrixBOSS: " + to_string(matrixboss_plain.serialize(matrixboss_out.stream) * 8.0 / n_nodes) + " bits per column", LogLevel::MAJOR);
}