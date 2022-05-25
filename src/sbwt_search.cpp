#include "libwheeler/BOSS.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "globals.hh"
#include "NodeBOSS.hh"
#include "SubsetWT.hh"
#include "stdlib_printing.hh"
#include "input_reading.hh"
#include "SubsetMatrixRank.hh"
#include "buffered_streams.hh"
#include "variants.hh"
#include <filesystem>
#include <cstdio>

using namespace std;
typedef long long LL;

// Assumes values of v are -1 or larger
inline void print_vector(const vector<int64_t>& v, Buffered_ofstream& out){
    // Fast manual integer-to-string conversion
    char buffer[32];
    char newline = '\n';
    for(int64_t x : v){
        LL i = 0;
        if(x == -1){
            buffer[0] = '1';
            buffer[1] = '-';
            i = 2;
        } else{
            while(x > 0){
                buffer[i++] = '0' + (x % 10);
                x /= 10;
            }
        }
        std::reverse(buffer, buffer + i);
        buffer[i] = ' ';
        out.write(buffer, i+1);
    }
    out.write(&newline, 1);
}

template<typename sbwt_t>
LL run_queries_streaming(Sequence_Reader_Buffered& sr, const string& outfile, const sbwt_t& sbwt, bool colex){
    write_log("Running streaming queries", LogLevel::MAJOR);
    Buffered_ofstream out(outfile);
    
    LL total_micros = 0;
    LL number_of_queries = 0;
    while(true){ 
        LL len = sr.get_next_read_to_buffer();
        if(len == 0) break;

        if(!colex) std::reverse(sr.read_buf, sr.read_buf + len);

        LL t0 = cur_time_micros();
        vector<int64_t> out_buffer = sbwt.streaming_search(sr.read_buf, len);
        total_micros += cur_time_micros() - t0;

        number_of_queries += out_buffer.size();

        // Write out
        if(!colex) std::reverse(out_buffer.begin(), out_buffer.end());
        print_vector(out_buffer, out);
    }
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    return number_of_queries;
}

// Returns number of queries executed
template<typename sbwt_t>
LL run_queries(Sequence_Reader_Buffered& sr, const string& outfile, const sbwt_t& sbwt, bool colex){
    vector<int64_t> out_buffer;
    if(sbwt.has_streaming_query_support()){
        return run_queries_streaming(sr, outfile, sbwt, colex);
    } else{
        LL total_micros = 0;
        LL number_of_queries = 0;
        write_log("Running queries", LogLevel::MAJOR);
        Buffered_ofstream out(outfile);
        LL k = sbwt.k;
        while(true){ 
            LL len = sr.get_next_read_to_buffer();
            if(len == 0) break;

            if(!colex) std::reverse(sr.read_buf, sr.read_buf + len);
            for(LL i = 0; i < len - k + 1; i++){
                LL t0 = cur_time_micros();
                LL ans = sbwt.search(sr.read_buf + i);
                total_micros += cur_time_micros() - t0;
                number_of_queries++;
                out_buffer.push_back(ans);
            }
            if(!colex) std::reverse(out_buffer.begin(), out_buffer.end());

            print_vector(out_buffer, out);
            out_buffer.clear();
        }
        write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
        return number_of_queries;
    }
    
}


int search_main(int argc, char** argv){

    LL micros_start = cur_time_micros();

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Query all k-mers of all input reads.");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,index-file", "Index input file.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA or FASTQ format. Multi-line FASTQ is not supported.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string outfile = opts["out-file"].as<string>();
    string indexfile = opts["index-file"].as<string>();
    string queryfile = opts["query-file"].as<string>();

    check_writable(outfile);
    check_readable(indexfile);
    check_readable(queryfile);
    
    Sequence_Reader_Buffered sr(queryfile);

    vector<string> variants = get_available_variants();

    throwing_ifstream in(indexfile, ios::binary);
    string variant = load_string(in.stream); // read variant type
    char colex; in.stream.read(&colex, 1); // Read colex flag
    if(std::find(variants.begin(), variants.end(), variant) == variants.end()){
        cerr << "Error loading index from file: unrecognized variant specified in the file" << endl;
        return 1;
    }

    write_log("Loading the index variant " + variant, LogLevel::MAJOR);
    LL number_of_queries = 0;

    if (variant == "plain-matrix"){
        plain_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(sr, outfile, sbwt, colex);
    }
    if (variant == "rrr-matrix"){
        rrr_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(sr, outfile, sbwt, colex);
    }
    if (variant == "mef-matrix"){
        mef_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(sr, outfile, sbwt, colex);
    }
    if (variant == "plain-split"){
        plain_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(sr, outfile, sbwt, colex);
    }
    if (variant == "rrr-split"){
        rrr_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(sr, outfile, sbwt, colex);
    }
    if (variant == "mef-split"){
        mef_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(sr, outfile, sbwt, colex);
    }
    if (variant == "plain-concat"){
        plain_concat_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(sr, outfile, sbwt, colex);
    }
    if (variant == "mef-concat"){
        mef_concat_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(sr, outfile, sbwt, colex);
    }
    if (variant == "plain-subsetwt"){
        plain_sswt_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(sr, outfile, sbwt, colex);
    }
    if (variant == "rrr-subsetwt"){
        rrr_sswt_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(sr, outfile, sbwt, colex);
    }

    LL total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);

    return 0;

}

