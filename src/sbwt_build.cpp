
#include "globals.hh"
#include "throwing_streams.hh"
#include "cxxopts.hpp"
#include "BOSS.hh"
#include "NodeBOSS.hh"
#include "SubsetMatrixRank.hh"
#include "input_reading.hh"
#include "variants.hh"

typedef long long LL;
using namespace std;

int build_main(int argc, char** argv){

    set_log_level(LogLevel::MAJOR);

    cxxopts::Options options(argv[0], "Construct an SBWT variant.");

    vector<string> variants = get_available_variants();
    string all_variants_string;
    for(string variant : variants) all_variants_string += " " + variant;

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("variant", "The SBWT variant to build. Available variants:" + all_variants_string, cxxopts::value<string>())
        ("streaming-support", "Build the auxiliary bit vector for streaming query support.", cxxopts::value<bool>()->default_value("false"))
        ("in-fasta", "Build in internal memory from a FASTA file (takes a lot of memory).", cxxopts::value<string>()->default_value(""))
        ("in-themisto", "Build from a Themisto .tdbg file.", cxxopts::value<string>()->default_value(""))
        ("k", "Value of k (must not be given if --in-themisto is given because themisto defines the k)", cxxopts::value<LL>()->default_value("0"))
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage example: " << argv[0] << " -o temp/out.sbwt --variant plain-matrix --in-fasta example_data/coli3.fna -k 30" << endl;
        exit(1);
    }

    string variant = opts["variant"].as<string>();
    if(std::find(variants.begin(), variants.end(), variant) == variants.end()){
        cerr << "Error: unknown variant: " << variant << endl;
        cerr << "Available variants are:" << all_variants_string << endl;
        return 1;
    }

    string out_file = opts["out-file"].as<string>();
    check_writable(out_file);

    string in_fasta = opts["in-fasta"].as<string>();
    string in_themisto = opts["in-themisto"].as<string>();
    bool streaming_support = opts["streaming-support"].as<bool>();

    if(in_fasta == "" && in_themisto == ""){
        cerr << "Error: No input file given" << endl;
        return 1;
    }
    
    if(in_fasta != "" && in_themisto != ""){
        cerr << "Error: Please give only one of the options --in-fasta and --in-themisto" << endl;
        return 1;
    }

    plain_matrix_sbwt_t matrixboss_plain;

    LL k = 0;

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
        k = wheelerBOSS.get_k();

        write_log("Themisto WheelerBOSS has order k = " + to_string(k) + " (node label length)", LogLevel::MAJOR);   
        write_log("Building MatrixBOSS representation", LogLevel::MAJOR);
        matrixboss_plain.build_from_WheelerBOSS(wheelerBOSS, streaming_support);
    } else{
        if(opts["k"].as<LL>() == 0){
            cerr << "Error: -k not specified" << endl;
            return 1;
        }
        k = opts["k"].as<LL>();
        check_readable(in_fasta);

        write_log("Reading the input", LogLevel::MAJOR);
        vector<string> input;
        Sequence_Reader_Buffered sr(in_fasta, FASTA_MODE);
        while(true) { 
           LL len = sr.get_next_read_to_buffer();
           if(len == 0) break;
           input.push_back(string(sr.read_buf));
        }

        write_log("Building SBWT subset sequence in memory", LogLevel::MAJOR);
        matrixboss_plain.build_from_strings(input, k, streaming_support);
    }

    throwing_ofstream out(out_file, ios::binary);
    LL bytes_written = 0;
    bytes_written += serialize_string(variant, out.stream); // Write variant string to file
    write_log("Building subset rank support", LogLevel::MAJOR);

    
    sdsl::bit_vector& A_bits = matrixboss_plain.subset_rank.A_bits;
    sdsl::bit_vector& C_bits = matrixboss_plain.subset_rank.C_bits;
    sdsl::bit_vector& G_bits = matrixboss_plain.subset_rank.G_bits;
    sdsl::bit_vector& T_bits = matrixboss_plain.subset_rank.T_bits;
    if (variant == "plain-matrix"){
        bytes_written = matrixboss_plain.serialize(out.stream);
    }
    if (variant == "rrr-matrix"){
        rrr_matrix_sbwt_t sbwt;
        sbwt.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k, streaming_support);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "mef-matrix"){
        mef_matrix_sbwt_t sbwt;
        sbwt.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k, streaming_support);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "plain-split"){
        plain_split_sbwt_t sbwt;
        sbwt.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k, streaming_support);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "rrr-split"){
        rrr_split_sbwt_t sbwt;
        sbwt.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k, streaming_support);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "mef-split"){
        mef_split_sbwt_t sbwt;
        sbwt.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k, streaming_support);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "plain-concat"){
        plain_concat_sbwt_t sbwt;
        sbwt.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k, streaming_support);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "mef-concat"){
        mef_concat_sbwt_t sbwt;
        sbwt.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k, streaming_support);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "plain-subsetwt"){
        plain_sswt_sbwt_t sbwt;
        sbwt.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k, streaming_support);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "rrr-subsetwt"){
        rrr_sswt_sbwt_t sbwt;
        sbwt.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k, streaming_support);
        bytes_written = sbwt.serialize(out.stream);
    }

    write_log("Built variant " + variant + " to file " + out_file, LogLevel::MAJOR);
    write_log("Space on disk: " + to_string(bytes_written * 8.0 / matrixboss_plain.n_nodes) + " bits per column", LogLevel::MAJOR);

    return 0;
}