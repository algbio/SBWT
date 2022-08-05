#include "globals.hh"
#include "throwing_streams.hh"
#include "cxxopts.hpp"
#include "SBWT.hh"
#include "SubsetMatrixRank.hh"
#include "SeqIO.hh"
#include "variants.hh"
#include "commands.hh"
#include "throwing_streams.hh"

typedef long long LL;
using namespace std;
using namespace sbwt;

int main(int argc, char** argv){

    sbwt::set_log_level(sbwt::LogLevel::MAJOR);

    cxxopts::Options options(argv[0], "Construct an SBWT variant from a plain matrix SBWT.");

    options.add_options()
        ("i,in-prefix", "Input file prefix that was given to the parquet-to-raw-bit-vectors program.", cxxopts::value<string>())
        ("o,out-file", "Output file for the constructed plain matrix SBWT.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string out_file = opts["out-file"].as<string>();
    sbwt::check_writable(out_file);

    string in_prefix = opts["in-prefix"].as<string>();

    throwing_ifstream A_in(in_prefix + "_A_bits.bin");
    throwing_ifstream C_in(in_prefix + "_C_bits.bin");
    throwing_ifstream G_in(in_prefix + "_G_bits.bin");
    throwing_ifstream T_in(in_prefix + "_T_bits.bin");
    throwing_ifstream n_columns_in(in_prefix + "_n_columns.txt");
    throwing_ifstream has_root_in(in_prefix + "_has_root.txt");

    LL n_columns; n_columns_in.stream >> n_columns;
    string has_root; has_root_in.stream >> has_root;

    cerr << "Number of columns: " << n_columns << endl;
    cerr << "Has root: " << has_root << endl;
    if(has_root != "yes"){
        throw std::runtime_error("Error: Todo: has_root == false not implemented");
    }

/*
    sbwt::plain_matrix_sbwt_t matrixboss_plain;
    write_log("Reading input.", sbwt::LogLevel::MAJOR);    
    matrixboss_plain.load(in.stream);

    sbwt::write_log("Building variant " + variant, sbwt::LogLevel::MAJOR);
    
    const sdsl::bit_vector& A_bits = matrixboss_plain.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = matrixboss_plain.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = matrixboss_plain.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = matrixboss_plain.get_subset_rank_structure().T_bits;
    const sdsl::bit_vector& ssupport = matrixboss_plain.get_streaming_support();
    LL n_kmers = matrixboss_plain.number_of_kmers();
    LL k = matrixboss_plain.get_k();

    LL bytes_written = 0;
    sbwt::throwing_ofstream out(out_file, ios::binary);

    sbwt::serialize_string(variant, out.stream);

    sbwt::write_log("Built variant " + variant + " to file " + out_file, sbwt::LogLevel::MAJOR);
    sbwt::write_log("Space on disk: " + 
                    to_string(bytes_written * 8.0 / matrixboss_plain.number_of_subsets()) + " bits per column, " +
                    to_string(bytes_written * 8.0 / matrixboss_plain.number_of_kmers()) + " bits per k-mer" , 
                    sbwt::LogLevel::MAJOR);

    return 0;
*/
}
