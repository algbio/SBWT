#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "NodeBOSSInMemoryConstructor.hh"
#include "throwing_streams.hh"
#include "suffix_group_optimization.hh"
#include "kmc_construct.hh"
#include "libwheeler/BOSS.hh"
#include "globals.hh"
#include "Kmer.hh"
#include <map>

using namespace std;

namespace sbwt{

// Assumes that a root node always exists
template <typename subset_rank_t>
class SBWT{

public:

    struct BuildConfig{
        vector<string> input_files;
        int k = 30;
        bool add_reverse_complements = false;
        bool build_streaming_support = true;
        int n_threads = 1;
        int min_abundance = 1;
        int max_abundance = 1e9;
        int ram_gigas = 2;
        string temp_dir = ".";
    };

    subset_rank_t subset_rank;
    sdsl::bit_vector suffix_group_starts; // Used for streaming queries

    // C-array
    vector<int64_t> C;

    int64_t n_nodes;
    int64_t k;

    SBWT() : n_nodes(0), k(0) {}
    SBWT(const sdsl::bit_vector& A_bits, 
             const sdsl::bit_vector& C_bits, 
             const sdsl::bit_vector& G_bits, 
             const sdsl::bit_vector& T_bits, 
             const sdsl::bit_vector& streaming_support, // Streaming support may be empty
             int64_t k);
    SBWT(const BuildConfig& config);

    int64_t search(const string& kmer) const; // Search for std::string
    int64_t search(const char* S) const; // Search for C-string

    // Query for all k-mers in the input
    vector<int64_t> streaming_search(const string& input) const;
    vector<int64_t> streaming_search(const char* input, int64_t len) const;
    bool has_streaming_query_support() const {return suffix_group_starts.size() > 0;}
    
    int64_t serialize(ostream& out) const; // Returns the number of bytes written
    int64_t serialize(const string& filename) const; // Returns the number of bytes written
    void load(istream& in);
    void load(const string& filename);

};


template <typename subset_rank_t>
SBWT<subset_rank_t>::SBWT(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits, const sdsl::bit_vector& streaming_support, int64_t k){
    subset_rank = subset_rank_t(A_bits, C_bits, G_bits, T_bits);

    this->n_nodes = A_bits.size();
    this->k = k;
    this->suffix_group_starts = streaming_support;

    // Get the C-array
    C.clear(); C.resize(4);
    C[0] = 1; // There is one incoming ghost-dollar to the root node
    C[1] = C[0] + subset_rank.rank(n_nodes, 'A');
    C[2] = C[1] + subset_rank.rank(n_nodes, 'C');
    C[3] = C[2] + subset_rank.rank(n_nodes, 'G');

}

template <typename subset_rank_t>
SBWT<subset_rank_t>::SBWT(const BuildConfig& config){
    string old_temp_dir = get_temp_file_manager().get_dir();
    get_temp_file_manager().set_dir(config.temp_dir);

    NodeBOSSKMCConstructor<SBWT<subset_rank_t>> builder;
    builder.build(config.input_files, *this, config.k, config.n_threads, config.ram_gigas, config.build_streaming_support, config.min_abundance, config.max_abundance);

    get_temp_file_manager().set_dir(old_temp_dir); // Return the old temporary directory
}

template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::search(const string& kmer) const{
    assert(kmer.size() == k);
    return search(kmer.c_str());
}

template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::search(const char* kmer) const{
    int64_t node_left = 0;
    int64_t node_right = n_nodes-1;
    for(int64_t i = 0; i < k; i++){

        char char_idx = 0;
        if(toupper(kmer[i]) == 'A') char_idx = 0;
        else if(toupper(kmer[i]) == 'C') char_idx = 1;
        else if(toupper(kmer[i]) == 'G') char_idx = 2;
        else if(toupper(kmer[i]) == 'T') char_idx = 3;
        else return -1; // Invalid character

        node_left = C[char_idx] + subset_rank.rank(node_left, kmer[i]);
        node_right = C[char_idx] + subset_rank.rank(node_right+1, kmer[i]) - 1;

        if(node_left > node_right) return -1; // Not found
    }
    if(node_left != node_right){
        cerr << "Bug: node_left != node_right" << endl;
        exit(1);
    }
    return node_left;
}

// Utility function: Serialization for a std::vector
// Returns number of bytes written
template<typename T>
int64_t serialize_std_vector(const vector<T>& v, ostream& os){
    // Write C-array
    int64_t n_bytes = sizeof(T) * v.size();
    os.write((char*)&n_bytes, sizeof(n_bytes));
    os.write((char*)v.data(), n_bytes);
    return sizeof(n_bytes) + n_bytes;
}

template<typename T>
vector<T> load_std_vector(istream& is){
    int64_t n_bytes = 0;
    is.read((char*)&n_bytes, sizeof(n_bytes));
    assert(n_bytes % sizeof(T) == 0);
    vector<T> v(n_bytes / sizeof(T));
    is.read((char*)(v.data()), n_bytes);
    return v;
}


template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::serialize(ostream& os) const{
    int64_t written = 0;
    written += subset_rank.serialize(os);
    written += suffix_group_starts.serialize(os);

    written += serialize_std_vector(C, os);

    // Write number of nodes
    os.write((char*)&n_nodes, sizeof(n_nodes));
    written += sizeof(n_nodes);

    // Write k
    os.write((char*)&k, sizeof(k));
    written += sizeof(k);

    return written;
}

template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::serialize(const string& filename) const{
    throwing_ofstream out(filename, ios::binary);
    return serialize(out.stream);
}


template <typename subset_rank_t>
void SBWT<subset_rank_t>::load(istream& is){
    subset_rank.load(is);
    suffix_group_starts.load(is);
    C = load_std_vector<int64_t>(is);
    is.read((char*)&n_nodes, sizeof(n_nodes));
    is.read((char*)&k, sizeof(k));
}

template <typename subset_rank_t>
void SBWT<subset_rank_t>::load(const string& filename){
    throwing_ifstream in(filename, ios::binary);
    load(in.stream);
}

template <typename subset_rank_t>
vector<int64_t> SBWT<subset_rank_t>::streaming_search(const char* input, int64_t len) const{
    if(suffix_group_starts.size() == 0)
        throw std::runtime_error("Error: streaming search support not built");
    
    vector<int64_t> ans;
    if(len < k) return ans;

    ans.push_back(search(input)); // Search the first k-mer
    for(int64_t i = 1; i < len - k + 1; i++){
        if(ans.back() == -1){
            // Need to search from scratch
            ans.push_back(search(input + i));
        } else{
            // Got to the start of the suffix group and do one search iteration
            int64_t colex = ans.back();
            while(suffix_group_starts[colex] == 0) colex--; // can not go negative because the first column is always marked

            char char_idx = -1;
            char c = toupper(input[i+k-1]);
            if(c == 'A') char_idx = 0;
            else if(c == 'C') char_idx = 1;
            else if(c == 'G') char_idx = 2;
            else if(c == 'T') char_idx = 3;
        
            if(char_idx == -1) ans.push_back(-1); // Not found
            else{
                int64_t node_left = colex;
                int64_t node_right = colex;
                node_left = C[char_idx] + subset_rank.rank(node_left, c);
                node_right = C[char_idx] + subset_rank.rank(node_right+1, c) - 1;
                if(node_left == node_right) ans.push_back(node_left);
                else ans.push_back(-1);
                // Todo: could save one subset rank query if we have fast access to the SBWT columns
            }
        }
    }
    return ans;
}

template <typename subset_rank_t>
vector<int64_t> SBWT<subset_rank_t>::streaming_search(const string& input) const{
    return streaming_search(input.c_str(), input.size());
}

} // namespace sbwt