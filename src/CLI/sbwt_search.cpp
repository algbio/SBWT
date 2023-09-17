#include <string>
#include <cstring>
#include "cxxopts.hpp"
#include "globals.hh"
#include "SBWT.hh"
#include "SubsetWT.hh"
#include "stdlib_printing.hh"
#include "SeqIO/SeqIO.hh"
#include "SubsetMatrixRank.hh"
#include "SeqIO/buffered_streams.hh"
#include "variants.hh"
#include "commands.hh"
#include <filesystem>
#include <cstdio>

using namespace std;

using namespace sbwt;

// Assumes values of v are -1 or larger
template <typename writer_t>
inline void print_vector(const vector<int64_t>& v, writer_t& out){
    // Fast manual integer-to-string conversion
    char buffer[32];
    char newline = '\n';
    for(int64_t x : v){
        int64_t i = 0;
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

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_streaming(reader_t& reader, writer_t& writer, const sbwt_t& sbwt){

    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;

        int64_t t0 = cur_time_micros();
        vector<int64_t> out_buffer = sbwt.streaming_search(reader.read_buf, len);
        total_micros += cur_time_micros() - t0;

        number_of_queries += out_buffer.size();

        // Write out
        print_vector(out_buffer, writer);
    }
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_not_streaming(reader_t& reader, writer_t& writer, const sbwt_t& sbwt){

    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    int64_t k = sbwt.get_k();
    vector<int64_t> out_buffer;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;

        for(int64_t i = 0; i < len - k + 1; i++){
            int64_t t0 = cur_time_micros();
            int64_t ans = sbwt.search(reader.read_buf + i);
            total_micros += cur_time_micros() - t0;
            number_of_queries++;
            out_buffer.push_back(ans);
        }

        print_vector(out_buffer, writer);
        out_buffer.clear();
    }
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_file(const string& infile, const string& outfile, const sbwt_t& sbwt){
    reader_t reader(infile);
    writer_t writer(outfile);
    if(sbwt.has_streaming_query_support()){
        write_log("Running streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
    else{
        write_log("Running non-streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_not_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
}

// Returns number of queries executed
template<typename sbwt_t>
int64_t run_queries(const vector<string>& infiles, const vector<string>& outfiles, const sbwt_t& sbwt, bool gzip_output){

    if(infiles.size() != outfiles.size()){
        string count1 = to_string(infiles.size());
        string count2 = to_string(outfiles.size());
        throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
    }

    typedef seq_io::Reader<seq_io::Buffered_ifstream<seq_io::zstr::ifstream>> in_gzip;
    typedef seq_io::Reader<seq_io::Buffered_ifstream<std::ifstream>> in_no_gzip;

    typedef seq_io::Buffered_ofstream<seq_io::zstr::ofstream> out_gzip;
    typedef seq_io::Buffered_ofstream<std::ofstream> out_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = seq_io::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input && gzip_output){
            n_queries_run += run_file<sbwt_t, in_gzip, out_gzip>(infiles[i], outfiles[i], sbwt);
        }
        if(gzip_input && !gzip_output){
            n_queries_run += run_file<sbwt_t, in_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt);
        }
        if(!gzip_input && gzip_output){
            n_queries_run += run_file<sbwt_t, in_no_gzip, out_gzip>(infiles[i], outfiles[i], sbwt);
        }
        if(!gzip_input && !gzip_output){
            n_queries_run += run_file<sbwt_t, in_no_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt);
        }
    }
    return n_queries_run;

}

int search_main(int argc, char** argv){

    int64_t micros_start = cur_time_micros();

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Query all k-mers of all input reads.");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,index-file", "Index input file.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
        ("z,gzip-output", "Writes output in gzipped form. This can shrink the output files by an order of magnitude.", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string indexfile = opts["index-file"].as<string>();
    check_readable(indexfile);

    // Interpret input file
    string queryfile = opts["query-file"].as<string>();
    vector<string> input_files;
    bool multi_file = queryfile.size() >= 4 && queryfile.substr(queryfile.size() - 4) == ".txt";
    if(multi_file){
        input_files = readlines(queryfile);
    } else{
        input_files = {queryfile};
    }
    for(string file : input_files) check_readable(file);

    // Interpret output file
    string outfile = opts["out-file"].as<string>();
    bool gzip_output = opts["gzip-output"].as<bool>();
    vector<string> output_files;
    if(multi_file){
        output_files = readlines(outfile);
    } else{
        output_files = {outfile};
    }
    for(string file : output_files) check_writable(file);

    vector<string> variants = get_available_variants();

    throwing_ifstream in(indexfile, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(std::find(variants.begin(), variants.end(), variant) == variants.end()){
        cerr << "Error loading index from file: unrecognized variant specified in the file" << endl;
        return 1;
    }

    write_log("Loading the index variant " + variant, LogLevel::MAJOR);
    int64_t number_of_queries = 0;

    if (variant == "plain-matrix"){
        plain_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    if (variant == "rrr-matrix"){
        rrr_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    if (variant == "mef-matrix"){
        mef_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    if (variant == "plain-split"){
        plain_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    if (variant == "rrr-split"){
        rrr_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    if (variant == "mef-split"){
        mef_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    if (variant == "plain-concat"){
        plain_concat_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    if (variant == "mef-concat"){
        mef_concat_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    if (variant == "plain-subsetwt"){
        plain_sswt_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    if (variant == "rrr-subsetwt"){
        rrr_sswt_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }

    int64_t total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);

    return 0;

}
