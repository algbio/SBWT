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

sdsl::bit_vector read_raw_bit_vector(string filename, LL n_bits){
    sdsl::bit_vector bv(n_bits, 0); // Reading the bits into here

    Buffered_ifstream<> in(filename, ios::binary);

    LL n_bytes = n_bits / 8 + (n_bits % 8 > 0); // ceil(n_bits/8)
    unsigned char byte;
    LL bits_read = 0;
    for(LL i = 0; i < n_bytes; i++){
        in.read((char*)(&byte), 1);
        for(LL j = 0; j < 8; j++){
            if(bits_read < n_bits){
                bv[bits_read++] = (byte & (1 << (7-j))) > 0; // Extract the j-th bit from left
            }
        }
    }

    return bv;
    
}

int main(int argc, char** argv){

    sbwt::set_log_level(sbwt::LogLevel::MAJOR);

    cxxopts::Options options(argv[0], "Construct an SBWT variant from a plain matrix SBWT.");

    options.add_options()
        ("i,in-prefix", "Input file prefix that was given to the parquet-to-raw-bit-vectors program.", cxxopts::value<string>())
        ("o,out-file", "Output file for the constructed plain matrix SBWT.", cxxopts::value<string>())
        ("k", "The k-mer k (node length).", cxxopts::value<LL>())
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

    LL k = opts["k"].as<LL>();

    throwing_ifstream n_columns_in(in_prefix + "_n_columns.txt");
    throwing_ifstream n_kmers_in(in_prefix + "_n_kmers.txt");
    throwing_ifstream has_root_in(in_prefix + "_has_root.txt");
    
    LL n_columns; n_columns_in.stream >> n_columns;
    LL n_kmers; n_kmers_in.stream >> n_kmers;
    string has_root; has_root_in.stream >> has_root;

    cerr << "Number of columns: " << n_columns << endl;
    cerr << "Has root: " << has_root << endl;
    if(has_root != "yes"){
        throw std::runtime_error("Error: Todo: has_root == false not implemented");
    }

    sdsl::bit_vector A_bits = read_raw_bit_vector(in_prefix + "_A_bits.bin", n_columns);
    sdsl::bit_vector C_bits = read_raw_bit_vector(in_prefix + "_C_bits.bin", n_columns);
    sdsl::bit_vector G_bits = read_raw_bit_vector(in_prefix + "_G_bits.bin", n_columns);
    sdsl::bit_vector T_bits = read_raw_bit_vector(in_prefix + "_T_bits.bin", n_columns);

    sdsl::bit_vector empty;

    sbwt::plain_matrix_sbwt_t matrixboss_plain(A_bits, C_bits, G_bits, T_bits, empty, k, n_kmers);

    throwing_ofstream out(out_file, ios::binary);
    
    sbwt::serialize_string("plain-matrix", out.stream);
    matrixboss_plain.serialize(out.stream);

    return 0;

}
