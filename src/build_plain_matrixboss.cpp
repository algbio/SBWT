
#include "globals.hh"
#include "throwing_streams.hh"
#include "cxxopts.hpp"
#include "BOSS.hh"
#include "NodeBOSS.hh"
#include "SubsetMatrixRank.hh"
#include "input_reading.hh"

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
        ("in-fasta", "Build in internal memory from a FASTA file (takes a lot of memory).", cxxopts::value<string>()->default_value(""))
        ("in-themisto", "Build from a Themisto .tdbg file.", cxxopts::value<string>()->default_value(""))
        ("k", "Value of k (must not be given if --in-themisto is given because themisto defines the k)", cxxopts::value<LL>()->default_value("0"))
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
    check_writable(out_file);

    string in_fasta = opts["in-fasta"].as<string>();
    string in_themisto = opts["in-themisto"].as<string>();

    if(in_fasta == "" && in_themisto == ""){
        cerr << "Error: No input file given" << endl;
        return 1;
    }
    
    if(in_fasta != "" && in_themisto != ""){
        cerr << "Error: Please give only one of the options --in-fasta and --in-themisto" << endl;
        return 1;
    }

    if(in_themisto != ""){
        check_readable(in_themisto);
        if(opts["k"].as<LL>() != 0){
            cerr << "Error: -k must not be given if building from Themisto because Themisto defines the k" << endl;
            return 1;
        }

        write_log("Loading the DBG", LogLevel::MAJOR);
        BOSS<sdsl::bit_vector> wheelerBOSS;
        throwing_ifstream in(in_themisto, ios::binary);
        wheelerBOSS.load(in.stream);
        LL n_nodes = wheelerBOSS.number_of_nodes();

        write_log("WheelerBOSS has order k = " + to_string(wheelerBOSS.get_k()) + " (node label length)", LogLevel::MAJOR);
        
        NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss_plain;
        write_log("Building MatrixBOSS representation", LogLevel::MAJOR);
        matrixboss_plain.build_from_WheelerBOSS(wheelerBOSS);

        write_log("Writing plain matrixBOSS to disk", LogLevel::MAJOR);
        throwing_ofstream matrixboss_out(out_file, ios::binary);
        write_log("MatrixBOSS: " + to_string(matrixboss_plain.serialize(matrixboss_out.stream) * 8.0 / n_nodes) + " bits per column", LogLevel::MAJOR);
    } else{
        if(opts["k"].as<LL>() == 0){
            cerr << "Error: -k not specified" << endl;
            return 1;
        }
        LL k = opts["k"].as<LL>();
        check_readable(in_fasta);
        NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss_plain;

        write_log("Reading the input", LogLevel::MAJOR);
        vector<string> input;
        Sequence_Reader_Buffered sr(in_fasta, FASTA_MODE);
        while(true) { 
           LL len = sr.get_next_read_to_buffer();
           if(len == 0) break;
           input.push_back(string(sr.read_buf));
        }

        write_log("Building plain MatrisBOSS in memory", LogLevel::MAJOR);
        matrixboss_plain.build_from_strings(input, k);

        write_log("Writing plain matrixBOSS to disk", LogLevel::MAJOR);
        throwing_ofstream matrixboss_out(out_file, ios::binary);
        write_log("MatrixBOSS: " + to_string(matrixboss_plain.serialize(matrixboss_out.stream) * 8.0 / matrixboss_plain.n_nodes) + " bits per column", LogLevel::MAJOR);
    }
}