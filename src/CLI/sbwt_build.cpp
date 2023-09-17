#include "globals.hh"
#include "throwing_streams.hh"
#include "cxxopts.hpp"
#include "SBWT.hh"
#include "SubsetMatrixRank.hh"
#include "SeqIO/SeqIO.hh"
#include "variants.hh"
#include "commands.hh"


using namespace std;

std::vector<std::string> get_available_variants(){
    return {"plain-matrix", "rrr-matrix", "mef-matrix", "plain-split", "rrr-split", "mef-split", "plain-concat", "mef-concat", "plain-subsetwt", "rrr-subsetwt"};
}

// Return the format, or throws if not all files have the same format
seq_io::FileFormat check_that_all_files_have_the_same_format(const vector<string>& filenames){
    if(filenames.size() == 0) runtime_error("Error: empty input file list");
    seq_io::FileFormat f1 = seq_io::figure_out_file_format(filenames[0]);
    for(int64_t i = 1; i < filenames.size(); i++){
        seq_io::FileFormat f2 = seq_io::figure_out_file_format(filenames[i]);
        if(f1.format != f2.format || f1.gzipped != f2.gzipped){
            throw runtime_error("Error: not all input files have the same format (" + filenames[0] + " vs " + filenames[i] + ")");
        }
    }
    return f1;
}

int build_main(int argc, char** argv){

    sbwt::set_log_level(sbwt::LogLevel::MAJOR);

    cxxopts::Options options(argv[0], "Construct an SBWT variant.");

    vector<string> variants = get_available_variants();
    string all_variants_string;
    for(string variant : variants) all_variants_string += " " + variant;

    options.add_options()
        ("i,in-file", "The input sequences as a FASTA or FASTQ file, possibly gzipped. If the file extension is .txt, the file is interpreted as a list of input files, one file on each line. All input files must be in the same format.", cxxopts::value<string>())
        ("o,out-file", "Output file for the constructed index.", cxxopts::value<string>())
        ("k,kmer-length", "The k-mer length.", cxxopts::value<int64_t>())
        ("p,precalc-length", "Precalculate SBWT intervals of strings of this length. Speeds up query, but takes 4^(p+2) bytes of memory.", cxxopts::value<int64_t>()->default_value("8"))
        ("variant", "The SBWT variant to build. Available variants:" + all_variants_string, cxxopts::value<string>()->default_value("plain-matrix"))
        ("add-reverse-complements", "Also add the reverse complement of every k-mer to the index. Warning: this creates a temporary reverse-complemented duplicate of each input file before construction. Make sure that the directory at --temp-dir can handle this amount of data. If the input is gzipped, the duplicate will also be compressed, which might take a while.", cxxopts::value<bool>()->default_value("false"))
        ("no-streaming-support", "Save space by not building the streaming query support bit vector. This leads to slower queries.", cxxopts::value<bool>()->default_value("false"))
        ("t,n-threads", "Number of parallel threads.", cxxopts::value<int64_t>()->default_value("1"))
        ("a,min-abundance", "Discard all k-mers occurring fewer than this many times. By default we keep all k-mers. Note that we consider a k-mer distinct from its reverse complement.", cxxopts::value<int64_t>()->default_value("1"))
        ("b,max-abundance", "Discard all k-mers occurring more than this many times.", cxxopts::value<int64_t>()->default_value("1000000000"))
        ("m,ram-gigas", "RAM budget in gigabytes (not strictly enforced). Must be at least 2.", cxxopts::value<int64_t>()->default_value("2"))
        ("d,temp-dir", "Location for temporary files.", cxxopts::value<string>()->default_value("."))
        ("v,verbose", "Print more verbose output.", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
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
    sbwt::check_writable(out_file);

    string in_file = opts["in-file"].as<string>();
    vector<string> input_files;
    if(in_file.size() >= 4 && in_file.substr(in_file.size() - 4) == ".txt"){
        input_files = sbwt::readlines(in_file);
    } else{
        input_files = {in_file};
    }
    for(string file : input_files) sbwt::check_readable(file);

    bool streaming_support = !(opts["no-streaming-support"].as<bool>());
    bool revcomps = opts["add-reverse-complements"].as<bool>();
    bool verbose = opts["verbose"].as<bool>();
    int64_t n_threads = opts["n-threads"].as<int64_t>();
    int64_t ram_gigas = opts["ram-gigas"].as<int64_t>();
    int64_t k = opts["k"].as<int64_t>();
    int64_t min_abundance = opts["min-abundance"].as<int64_t>();
    int64_t max_abundance = opts["max-abundance"].as<int64_t>();
    int64_t precalc_length = opts["precalc-length"].as<int64_t>();
    string temp_dir = opts["temp-dir"].as<string>();
    sbwt::get_temp_file_manager().set_dir(temp_dir);

    if(verbose){
        sbwt::set_log_level(sbwt::LogLevel::MINOR);
    }

    if(precalc_length > k){
        write_log("Warning: precalc length " + to_string(precalc_length) + " is longer than k = " + to_string(k), sbwt::LogLevel::MAJOR);
        write_log("Setting precalc length to " + to_string(k), sbwt::LogLevel::MAJOR);
        precalc_length = k;
    }

    seq_io::FileFormat fileformat = check_that_all_files_have_the_same_format(input_files);
    if(revcomps){
        sbwt::write_log("Creating a reverse-complemented version of each input file to " + temp_dir, sbwt::LogLevel::MAJOR);
        vector<string> new_files;
        for(int64_t i = 0; i < input_files.size(); i++){
            new_files.push_back(sbwt::get_temp_file_manager().create_filename("", fileformat.extension));
        }
        if(fileformat.gzipped){
            seq_io::create_reverse_complement_files<
                seq_io::Reader<seq_io::Buffered_ifstream<seq_io::zstr::ifstream>>,
                seq_io::Writer<seq_io::Buffered_ofstream<seq_io::zstr::ofstream>>>(input_files, new_files);
        } else{
            seq_io::create_reverse_complement_files<
                seq_io::Reader<seq_io::Buffered_ifstream<std::ifstream>>,
                seq_io::Writer<seq_io::Buffered_ofstream<std::ofstream>>>(input_files, new_files);
        }
        for(string f : new_files) input_files.push_back(f);
    }

    write_log("Building SBWT subset sequence using KMC", sbwt::LogLevel::MAJOR);
    sbwt::plain_matrix_sbwt_t::BuildConfig config;
    config.input_files = input_files;
    config.k = k;
    config.build_streaming_support = streaming_support;
    config.n_threads = n_threads;
    config.min_abundance = min_abundance;
    config.max_abundance = max_abundance;
    config.ram_gigas = ram_gigas;
    config.temp_dir = temp_dir;
    config.precalc_k = 0; // No precalc yet at this point to avoid doing it twice

    sbwt::plain_matrix_sbwt_t matrixboss_plain(config);

    sbwt::throwing_ofstream out(out_file, ios::binary);
    int64_t bytes_written = 0;
    bytes_written += sbwt::serialize_string(variant, out.stream); // Write variant string to file

    write_log("Build SBWT for " + to_string(matrixboss_plain.number_of_kmers()) + " distinct k-mers", sbwt::LogLevel::MAJOR);
    write_log("SBWT has " + to_string(matrixboss_plain.number_of_subsets()) + " subsets", sbwt::LogLevel::MAJOR);

    sbwt::write_log("Building subset rank support", sbwt::LogLevel::MAJOR);
    
    const sdsl::bit_vector& A_bits = matrixboss_plain.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = matrixboss_plain.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = matrixboss_plain.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = matrixboss_plain.get_subset_rank_structure().T_bits;
    const sdsl::bit_vector& ssupport = matrixboss_plain.get_streaming_support();
    int64_t n_kmers = matrixboss_plain.number_of_kmers();

    if (variant == "plain-matrix"){
        matrixboss_plain.do_kmer_prefix_precalc(precalc_length);
        bytes_written = matrixboss_plain.serialize(out.stream);
    }
    if (variant == "rrr-matrix"){
        sbwt::rrr_matrix_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_length);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "mef-matrix"){
        sbwt::mef_matrix_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_length);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "plain-split"){
        sbwt::plain_split_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_length);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "rrr-split"){
        sbwt::rrr_split_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_length);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "mef-split"){
        sbwt::mef_split_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_length);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "plain-concat"){
        sbwt::plain_concat_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_length);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "mef-concat"){
        sbwt::mef_concat_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_length);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "plain-subsetwt"){
        sbwt::plain_sswt_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_length);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "rrr-subsetwt"){
        sbwt::rrr_sswt_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_length);
        bytes_written = sbwt.serialize(out.stream);
    }

    sbwt::write_log("Built variant " + variant + " to file " + out_file, sbwt::LogLevel::MAJOR);
    sbwt::write_log("Space on disk: " + 
                    to_string(bytes_written * 8.0 / matrixboss_plain.number_of_subsets()) + " bits per column, " +
                    to_string(bytes_written * 8.0 / matrixboss_plain.number_of_kmers()) + " bits per k-mer" , 
                    sbwt::LogLevel::MAJOR);

    return 0;
}
