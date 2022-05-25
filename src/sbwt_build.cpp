
#include "globals.hh"
#include "throwing_streams.hh"
#include "cxxopts.hpp"
#include "BOSS.hh"
#include "NodeBOSS.hh"
#include "SubsetMatrixRank.hh"
#include "SeqIO.hh"
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
        ("i,in-file", "The input sequences as a FASTA or FASTQ file. If the file extension is .txt, the file is interpreted as a list of input files, one file on each line.", cxxopts::value<string>())
        ("o,out-file", "Output file for the constructed index.", cxxopts::value<string>())
        ("k,kmer-length", "The k-mer length.", cxxopts::value<LL>())
        ("variant", "The SBWT variant to build. Available variants:" + all_variants_string, cxxopts::value<string>()->default_value("plain-matrix"))
        ("add-reverse-complements", "Also add the reverse complement of every k-mer to the index (Warning: this creates a temporary reverse-complemented duplicate of each input file before construction. Make sure that the directory at --temp-dir can handle this amount of data).", cxxopts::value<bool>()->default_value("false"))
        ("no-streaming-support", "Save space by not building the streaming query support bit vector. This leads to slower queries.", cxxopts::value<bool>()->default_value("false"))
        ("t,n-threads", "Number of parallel threads.", cxxopts::value<LL>()->default_value("1"))
        ("a,min-abundance", "Discard all k-mers occurring fewer than this many times. By default we keep all k-mers. Note that we consider a k-mer distinct from its reverse complement.", cxxopts::value<LL>()->default_value("1"))
        ("m,ram-gigas", "RAM budget in gigabytes (not strictly enforced). Must be at least 2.", cxxopts::value<LL>()->default_value("2"))
        ("temp-dir", "Location for temporary files.", cxxopts::value<string>()->default_value("."))
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage example: " << argv[0] << " -i example_data/coli3.fna -o index.sbwt -k 30" << endl;
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

    string in_file = opts["in-file"].as<string>();
    vector<string> input_files;
    if(in_file.size() >= 4 && in_file.substr(in_file.size() - 4) == ".txt"){
        input_files = readlines(in_file);
    } else{
        input_files = {in_file};
    }
    for(string file : input_files) check_readable(file);

    bool streaming_support = !(opts["no-streaming-support"].as<bool>());
    bool revcomps = opts["add-reverse-complements"].as<bool>();
    LL n_threads = opts["n-threads"].as<LL>();
    LL ram_gigas = opts["ram-gigas"].as<LL>();
    LL k = opts["k"].as<LL>();
    LL min_abundance = opts["min-abundance"].as<LL>();
    string temp_dir = opts["temp-dir"].as<string>();
    get_temp_file_manager().set_dir(temp_dir);    

    if(revcomps){
        write_log("Creating a reverse-complemented version of each input file to " + temp_dir, LogLevel::MAJOR);
        vector<string> new_files = create_reverse_complement_files(input_files);
        for(string f : new_files) input_files.push_back(f);
    }

    plain_matrix_sbwt_t matrixboss_plain;

    write_log("Building SBWT subset sequence using KMC", LogLevel::MAJOR);
    matrixboss_plain.build_using_KMC(input_files, k, streaming_support, n_threads, ram_gigas, min_abundance);
    char colex = false; // Lexicographic or colexicographic index? KMC sorts in lexicographic order.

    throwing_ofstream out(out_file, ios::binary);
    LL bytes_written = 0;
    bytes_written += serialize_string(variant, out.stream); // Write variant string to file
    out.stream.write(&colex, 1); bytes_written++; // Write colex flag

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
