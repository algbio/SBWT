#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "globals.hh"
#include "SBWT.hh"
#include "SubsetWT.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include "Parquet_Kmer_Stream.hh"
#include "throwing_streams.hh"
#include "variants.hh"
#include <filesystem>

//#include "MEF.hpp"

using namespace sdsl;

using namespace std;
using namespace sbwt;

typedef long long LL;

class Output_Bit_Stream{

public:

    static const LL buf_cap = 1 << 27; // In BITS. MUST BE DIVISIBLE BY 8. 2^27 bits = 2^24 bytes = 16 MB

    ofstream out;
    unsigned char* buf;
    LL bits_in_buf = 0;

    Output_Bit_Stream(const string& filename) : out(filename, ios_base::binary){
        if(!out.good()) throw runtime_error("Error opening file: " + filename);
        buf = (unsigned char*)malloc(buf_cap/8); // ASSUMES BUFFER HAS A NUMBER OF BITS DIVISIBLE BY 8
        for(LL i = 0; i < buf_cap/8; i++) buf[i] = 0; // Zero-initialize
    }

    void add_bit(bool b){
        if(bits_in_buf == buf_cap) flush();
        LL byte_index = bits_in_buf / 8;
        LL byte_offset = bits_in_buf % 8;
        buf[byte_index] |= (b << (7 - byte_offset));
        bits_in_buf++;
    }

    void flush(){
        LL bytes = bits_in_buf / 8 + (bits_in_buf % 8 > 0); // ceil(bits_in_buf / 8)
        out.write((char*)buf, bytes);
        out.flush();
        for(LL i = 0; i < buf_cap/8; i++) buf[i] = 0;
        bits_in_buf = 0;
    }

    ~Output_Bit_Stream(){
        flush();
        free(buf);
    }

};

vector<string> split(string s, char delimiter){
    stringstream test(s);
    string segment;
    vector<string> seglist;

    while(getline(test, segment, delimiter))
    {
        seglist.push_back(segment);
    }
    return seglist;
}

// Format is like this:
// ATGAGGGCGCGCTCAATGGGCACGCCGTTT,1,AG
// ATTAGGGCGCGCTCAATGGGCACGCCGTTT,1,G
// ATGAGGGCGTGCTCAATGGGCACGCCGTTT,1,G

vector<string> read_lines(string filename){
    string line;
    throwing_ifstream in(filename);
    vector<string> lines;
    while(getline(in.stream, line)){
        lines.push_back(line);
    }
    return lines;
}

int main(int argc, char** argv){

    set_log_level(LogLevel::MINOR);

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "This program creates a plain matrix SBWT out of output of the Spark-Themisto code of Jaakko Vuohtoniemi.");

    options.add_options()
        ("i,in-directory", "Parquet file directory", cxxopts::value<string>())
        ("o,out-prefix", "Output file path prefix.", cxxopts::value<string>())
        ("k", "The k-mer k.", cxxopts::value<LL>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string in_directory = opts["in-directory"].as<string>();
    string output_prefix = opts["out-prefix"].as<string>();
    LL k = opts["k"].as<LL>();

    Parquet_Kmer_Stream stream(in_directory, k);

    string prev_kmer;
    LL start_of_suffix_group = 0;
    LL column_idx = 0;
    vector<vector<bool> > rows(4); // DNA alphabet
    bool root_kmer_exists = false;
    string kmer, kmer_outlabels;
    LL n_kmers_without_dummies = 0;

    bool A_bit = 0;
    bool C_bit = 0;
    bool G_bit = 0;
    bool T_bit = 0;

    Output_Bit_Stream A_out(output_prefix + "_A_bits.bin");
    Output_Bit_Stream C_out(output_prefix + "_C_bits.bin");
    Output_Bit_Stream G_out(output_prefix + "_G_bits.bin");
    Output_Bit_Stream T_out(output_prefix + "_T_bits.bin");
    throwing_ofstream n_columns_out(output_prefix + "_n_columns.txt");
    throwing_ofstream has_root_out(output_prefix + "_has_root.txt");

    // Writes a suffix group of given width bit bits set by A_bits, C_bits, G_bits and T_bits.
    auto write_suffix_group = [&](LL width){

        // Add the column for the first k-mer
        A_out.add_bit(A_bit);
        C_out.add_bit(C_bit);
        G_out.add_bit(G_bit);
        T_out.add_bit(T_bit);

        // Add zero-columns for the rest
        for(LL i = 0; i < width - 1; i++){
            A_out.add_bit(0);
            C_out.add_bit(0);
            G_out.add_bit(0);
            T_out.add_bit(0);
        }
    };

    while(stream.next(kmer, kmer_outlabels)){
        if(kmer.back() == '$') root_kmer_exists = true; // K-mer is all dollars if the last k-mer is a dollar
        if(std::find(kmer.begin(), kmer.end(), '$') == kmer.end()) n_kmers_without_dummies++;
        if(prev_kmer == "" || prev_kmer.substr(1) != kmer.substr(1)){
            // New suffix group starts

            if(column_idx != 0){
                // Write the previous suffix group
                write_suffix_group(column_idx - start_of_suffix_group);
                A_bit = 0; C_bit = 0; G_bit = 0; T_bit = 0;
            }

            start_of_suffix_group = column_idx;   
        }

        for(char c : kmer_outlabels){
            if(c == 'A') A_bit = 1;
            if(c == 'C') C_bit = 1;
            if(c == 'G') G_bit = 1;
            if(c == 'T') T_bit = 1;
        }

        column_idx++;
        prev_kmer = kmer;
    
    }

    // Write the last suffix group
    write_suffix_group(column_idx - start_of_suffix_group);
    A_bit = 0; C_bit = 0; G_bit = 0; T_bit = 0;

    n_columns_out.stream << column_idx << endl;
    has_root_out.stream << (root_kmer_exists ? "yes" : "no") << endl;

    std::cerr << "Root k-mer exists: " << (root_kmer_exists ? "true" : "false") << endl;

/*
    LL n = rows[0].size();
    vector<sdsl::bit_vector> sdsl_vectors;
    for(LL i = 0; i < 4; i++){
        sdsl::bit_vector vec;
        if(root_kmer_exists) vec.resize(n);
        else{
            vec.resize(n+1);
            vec[0] = 0;
        }
        for(LL j = 0; j < n; j++) vec[j + (!root_kmer_exists)] = rows[i][j];
        sdsl_vectors.push_back(vec);
    }

    sdsl::bit_vector empty;
    plain_matrix_sbwt_t matrixboss(sdsl_vectors[0], sdsl_vectors[1], sdsl_vectors[2], sdsl_vectors[3], empty, k, n_kmers_without_dummies);
    throwing_ofstream out(outfile, ios::binary);
    matrixboss.serialize(out.stream);
*/

}