#include <iostream>
#include "throwing_streams.hh"
#include "globals.hh"
#include "variants.hh"
#include "SeqIO/SeqIO.hh"

using namespace std;
using namespace sbwt;

void sequential_access_benchmark(const plain_matrix_sbwt_t& sbwt){
    int64_t max_queries = 1e5;
    write_log("Accessing up to " + to_string(max_queries) + " k-mers sequentially", LogLevel::MAJOR);
    int64_t k = sbwt.get_k();
    char buf[k];
    int64_t buf0_sum = 0; // To prevent the compiler from optimizing the queries away
    int64_t micros_start = cur_time_micros();
    int64_t n_queries = 0;
    for(int64_t i = 0; i < min((int64_t)1e5, sbwt.number_of_subsets()); i++){
        sbwt.get_kmer(i, buf);
        buf0_sum += buf[0];
        n_queries++;
    }
    int64_t micros_end = cur_time_micros();
    cout << "Total queries: " << n_queries << endl;
    cout << "Sequential access us / kmer: " << (double) (micros_end - micros_start) / n_queries << endl;
    cout << "Checksum: " << buf0_sum << endl;
}

void sequential_access_with_select_support_benchmark(const plain_matrix_sbwt_t& sbwt){
    write_log("Building select support", LogLevel::MAJOR);
    SubsetMatrixSelectSupport ss = sbwt.get_subset_rank_structure().build_select_support();
    int64_t max_queries = 1e5;
    write_log("Accessing up to " + to_string(max_queries) + " k-mers sequentially with select support", LogLevel::MAJOR);
    int64_t k = sbwt.get_k();
    char buf[k];
    int64_t buf0_sum = 0; // To prevent the compiler from optimizing the queries away
    int64_t micros_start = cur_time_micros();
    int64_t n_queries = 0;
    for(int64_t i = 0; i < min((int64_t)1e5, sbwt.number_of_subsets()); i++){
        sbwt.get_kmer_fast(i, buf, ss);
        buf0_sum += buf[0];
        n_queries++;
    }
    int64_t micros_end = cur_time_micros();
    cout << "Total queries: " << n_queries << endl;
    cout << "Sequential access us / kmer with select support: " << (double) (micros_end - micros_start) / n_queries << endl;
    cout << "Checksum: " << buf0_sum << endl;
}

void search_benchmark(const plain_matrix_sbwt_t& sbwt, const string& queryfile){
    int64_t k = sbwt.get_k();
    int64_t max_queries = 1e5;
    write_log("Searching up to " + to_string(max_queries) + " k-mers individually", LogLevel::MAJOR);
    seq_io::Reader<seq_io::Buffered_ifstream<std::ifstream>> reader(queryfile);

    int64_t micros_start = cur_time_micros();
    int64_t total_colex_rank = 0; // To prevent the compiler from optimizing the queries away
    int64_t n_queries = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        bool done = false;
        for(int64_t i = 0; i < len - k + 1; i++){
            total_colex_rank += sbwt.search(reader.read_buf + i);
            n_queries++;
            if(n_queries == max_queries){
                done = true;
                break;
            }
        }
        if(done) break;
    }
    int64_t micros_end = cur_time_micros();

    cout << "Total queries: " << n_queries << endl;
    cout << "Individual search us / kmer: " << (double) (micros_end - micros_start) / n_queries << endl;
    cout << "Checksum: " << total_colex_rank << endl;
}

void streaming_search_benchmark(const plain_matrix_sbwt_t& sbwt, const string& queryfile){
    int64_t k = sbwt.get_k();
    int64_t max_queries = 1e5;
    write_log("Searching up to " + to_string(max_queries) + " k-mers with streaming search", LogLevel::MAJOR);
    seq_io::Reader<seq_io::Buffered_ifstream<std::ifstream>> reader(queryfile);

    int64_t micros_start = cur_time_micros();
    int64_t total_colex_rank = 0; // To prevent the compiler from optimizing the queries away
    int64_t n_queries = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        for(int64_t x : sbwt.streaming_search(reader.read_buf, len)){
            n_queries++;
            total_colex_rank += x * (n_queries <= max_queries); // Don't add to checksum the extra queries
        }
        if(n_queries >= max_queries) break;
    }
    int64_t micros_end = cur_time_micros();

    cout << "Total queries: " << n_queries << endl;
    cout << "Individual search us / kmer: " << (double) (micros_end - micros_start) / n_queries << endl;
    cout << "Checksum: " << total_colex_rank << endl;

}

int main(int argc, char** argv){
    if(argc != 3){
        cerr << "Usage: " << argv[0] << " index.sbwt queries.fna" << endl;
        cerr << "Currently only supports the plain matrix variant and uncompressed queries" << endl;
        return 1;
    }

    string indexfile = argv[1];
    string queryfile = argv[2];

    throwing_ifstream in(indexfile, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if (variant != "plain-matrix"){
        cerr << "Error: only plain-matrix variant is supported currently" << endl;
        return 1;
    }

    write_log("Loading the index", LogLevel::MAJOR);
    plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);

    sequential_access_benchmark(sbwt);
    sequential_access_with_select_support_benchmark(sbwt);
    search_benchmark(sbwt, queryfile);
    streaming_search_benchmark(sbwt, queryfile);

}