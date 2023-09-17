#include "globals.hh"
#include "throwing_streams.hh"
#include "cxxopts.hpp"
#include "SBWT.hh"
#include "SubsetMatrixRank.hh"
#include "SeqIO/SeqIO.hh"
#include "variants.hh"
#include "commands.hh"


using namespace std;

int build_from_plain_main(int argc, char** argv){

    sbwt::set_log_level(sbwt::LogLevel::MAJOR);

    cxxopts::Options options(argv[0], "Construct an SBWT variant from a plain matrix SBWT.");

    vector<string> variants = get_available_variants();
    string all_variants_string;
    for(string variant : variants) all_variants_string += " " + variant;

    options.add_options()
        ("i,in-file", "Index file of a plain matrix SBWT.", cxxopts::value<string>())
        ("o,out-file", "Output file for the constructed variant.", cxxopts::value<string>())
        ("variant", "The SBWT variant to build. Available variants:" + all_variants_string, cxxopts::value<string>()->default_value("plain-matrix"))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
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
    sbwt::check_readable(in_file);

    sbwt::throwing_ifstream in(in_file, ios::binary);
    string variant_on_disk = sbwt::load_string(in.stream); // read variant type
    if(variant_on_disk != "plain-matrix"){
        cerr << "Error input is not a plain-matrix SBWT." << endl;
        return 1;
    }

    sbwt::plain_matrix_sbwt_t matrixboss_plain;
    write_log("Reading input.", sbwt::LogLevel::MAJOR);    
    matrixboss_plain.load(in.stream);

    sbwt::write_log("Building variant " + variant, sbwt::LogLevel::MAJOR);
    
    const sdsl::bit_vector& A_bits = matrixboss_plain.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = matrixboss_plain.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = matrixboss_plain.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = matrixboss_plain.get_subset_rank_structure().T_bits;
    const sdsl::bit_vector& ssupport = matrixboss_plain.get_streaming_support();
    int64_t n_kmers = matrixboss_plain.number_of_kmers();
    int64_t k = matrixboss_plain.get_k();
    int64_t precalc_k = matrixboss_plain.get_precalc_k();

    int64_t bytes_written = 0;
    sbwt::throwing_ofstream out(out_file, ios::binary);

    sbwt::serialize_string(variant, out.stream);
    if (variant == "plain-matrix"){
        bytes_written = matrixboss_plain.serialize(out.stream);
    }
    if (variant == "rrr-matrix"){
        sbwt::rrr_matrix_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_k);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "mef-matrix"){
        sbwt::mef_matrix_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_k);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "plain-split"){
        sbwt::plain_split_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_k);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "rrr-split"){
        sbwt::rrr_split_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_k);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "mef-split"){
        sbwt::mef_split_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_k);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "plain-concat"){
        sbwt::plain_concat_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_k);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "mef-concat"){
        sbwt::mef_concat_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_k);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "plain-subsetwt"){
        sbwt::plain_sswt_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_k);
        bytes_written = sbwt.serialize(out.stream);
    }
    if (variant == "rrr-subsetwt"){
        sbwt::rrr_sswt_sbwt_t sbwt(A_bits, C_bits, G_bits, T_bits, ssupport, k, n_kmers, precalc_k);
        bytes_written = sbwt.serialize(out.stream);
    }

    sbwt::write_log("Built variant " + variant + " to file " + out_file, sbwt::LogLevel::MAJOR);
    sbwt::write_log("Space on disk: " + 
                    to_string(bytes_written * 8.0 / matrixboss_plain.number_of_subsets()) + " bits per column, " +
                    to_string(bytes_written * 8.0 / matrixboss_plain.number_of_kmers()) + " bits per k-mer" , 
                    sbwt::LogLevel::MAJOR);

    return 0;
}
