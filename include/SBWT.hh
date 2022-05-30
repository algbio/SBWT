#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "NodeBOSSInMemoryConstructor.hh"
#include "throwing_streams.hh"
#include "suffix_group_optimization.hh"
#include "kmc_construct.hh"
#include "globals.hh"
#include "Kmer.hh"
#include <map>

/*

This file contains a class that implements the SBWT index described in the paper:

Alanko, J. N., Puglisi, S. J., & Vuohtoniemi, J. (2022). Succinct k-mer Set 
Representations Using Subset Rank Queries on the Spectral Burrows-Wheeler 
Transform (SBWT). bioRxiv.
*/

using namespace std;

namespace sbwt{

// Assumes that a root node always exists
template <typename subset_rank_t>
class SBWT{

private:


    // Is the index for the reverse (lex-sorted k-mers) or forward (colex-sorted k-mers)?
    // The paper describes the colex version, but if we construct the index via KMC, then
    // we get the lex version, because KMC sorts in lex-order. Then the search will go backward
    // instead of forward.
    bool colex;

    subset_rank_t subset_rank; // The subset rank query implementation
    sdsl::bit_vector suffix_group_starts; // Marks the first column of every suffix group (see paper)
    vector<int64_t> C; // The array of cumulative character counts
    int64_t n_nodes; // Number of nodes (= columns) in the data structure
    int64_t n_kmers; // Number of k-mers indexed in the data structure
    int64_t k; // The k-mer k

public:

    struct BuildConfig{
        vector<string> input_files; /** List of paths to input filenames. */
        int k = 30; /** Length of the k-mers. */
        bool add_reverse_complements = false; /** Whether we should also add the reverse complemented k-mers to the index. */
        bool build_streaming_support = true; /** Whether we should build the streaming query support. */
        int n_threads = 1; /** Number of parallel threads in construction. */
        int min_abundance = 1; /** k-mers occurring fewer than this many times are discarded. */
        int max_abundance = 1e9; /** k-mers occurring more than this many times are discarded */
        int ram_gigas = 2; /** RAM budget in gigabytes. Not strictly enforced. */
        string temp_dir = "."; /** Path to the directory for the temporary files. */
    };

    /**
     * @brief Construct an empty SBWT.
     * 
     */
    SBWT() : colex(true), n_nodes(0), n_kmers(0), k(0){}

    /**
     * @brief Construct SBWT from precomputed data.
     * 
     * @param A_bits Row of character A in the plain matrix SBWT.
     * @param C_bits Row of character C in the plain matrix SBWT.
     * @param G_bits Row of character G in the plain matrix SBWT.
     * @param T_bits Row of character T in the plain matrix SBWT.
     * @param streaming_support The streaming support bit vector. Can be empty.
     * @param k Length of the k-mers.
     * @param number_of_kmers Number of k-mers in the data structure.
     * @param colex Whether the index is colex- or lex-sorted.
     */
    SBWT(const sdsl::bit_vector& A_bits, 
         const sdsl::bit_vector& C_bits, 
         const sdsl::bit_vector& G_bits, 
         const sdsl::bit_vector& T_bits, 
         const sdsl::bit_vector& streaming_support, // Streaming support may be empty
         int64_t k, 
         int64_t number_of_kmers,
         bool colex);

    /**
     * @brief Construct SBWT using the KMC-based construction algorithm.
     * 
     * @param config construction paramters.
     */
    SBWT(const BuildConfig& config); 

    // Accessors

    /**
     * @brief Return whether this index uses colex- or lex-sorted k-mers.
     * 
     * @return true If colex-sorted.
     * @return false If lex-sorted.
     */
    bool is_colex() const {return colex;}

    /**
     * @brief Get a const reference to the internal subset rank structure.
     */
    const subset_rank_t& get_subset_rank_structure() const {return subset_rank;}

    /**
     * @brief Get a const reference to the internal streaming support bit vector.
     */
    const sdsl::bit_vector& get_streaming_support() const {return suffix_group_starts;}

    /**
     * @brief Get a const reference to the cumulative character count array.
     */
    const vector<int64_t>& get_C_array() const {return C;}

    /**
     * @brief Get the number of subsets in the SBWT (= number of columns in the plain matrix representation).
     */
    int64_t number_of_subsets() const {return n_nodes;}

    /**
     * @brief Get the number of k-mers indexed in the data structure.
     */
    int64_t number_of_kmers() const {return n_kmers;}

    /**
     * @brief Get the length of the k-mers.
     */
    int64_t get_k() const {return k;}

    /**
     * @brief Search for a k-mer as an std::string.
     * 
     * @param kmer The k-mer to search for. If the length of this is longer than k, then only the first k-mer is searched.
     * @return The rank of the k-mer in the data structure, or -1 if the k-mer is not in the index.
     * @see streaming_search()
     */
    int64_t search(const string& kmer) const;

    /**
     * @brief Search for a k-mer as C-string.
     * 
     * @param kmer The k-mer to search for.
     * @return The rank of the k-mer in the data structure, or -1 if the k-mer is not in the index.
     * @see streaming_search()
     */
    int64_t search(const char* kmer) const;

    /**
     * @brief Query all k-mers of the input std::string. Requires that the streaming support had been built.
     * 
     * @throws std::runtime_error If the streaming support has not been built.
     * @param input The input string 
     * @return vector<int64_t> The ranks of the k-mers of the input in the data structure, with -1 for those that are not found in the index. These result will be the same as if search() was called for each k-mer of the input from left to right in order.
     * @see search()
     */
    vector<int64_t> streaming_search(const string& input) const;

    /**
     * @brief Query all k-mers of the input C-string. Requires that the streaming support had been built.
     * 
     * @throws std::runtime_error If the streaming support has not been built.
     * @param input The input string 
     * @param len Length of the input string
     * @return vector<int64_t> The ranks of the k-mers of the input in the data structure, with -1 for those that are not found in the index. These result will be the same as if search() was called for each k-mer of the input from left to right in order.
     * @see search()
     */
    vector<int64_t> streaming_search(const char* input, int64_t len) const;

    /**
     * @brief Whether streaming support is built for the data structure.
     * 
     * @return true If streaming support has been built.
     * @return false If streaming support has not been built.
     */
    bool has_streaming_query_support() const {return suffix_group_starts.size() > 0;}
    
    /**
     * @brief Write the data structure into the given output stream.
     * 
     * @param out The output stream.
     * @return int64_t Number of bytes written.
     * @see load()
     */
    int64_t serialize(ostream& out) const; // Returns the number of bytes written

    /**
     * @brief Write the data structure into the given file.
     * 
     * @param filename The output file.
     * @return int64_t Number of bytes written.
     * @see load()
     */
    int64_t serialize(const string& filename) const; // Returns the number of bytes written

    /**
     * @brief Load the serialized data structure from an input stream.
     * 
     * @param in The input stream.
     * @see serialize()
     */
    void load(istream& in);

    /**
     * @brief Load the serialized data structure from an input file.
     * 
     * @param filename The input file.
     * @see serialize()
     */
    void load(const string& filename);

};


template <typename subset_rank_t>
SBWT<subset_rank_t>::SBWT(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits, const sdsl::bit_vector& streaming_support, int64_t k, int64_t n_kmers, bool colex){
    subset_rank = subset_rank_t(A_bits, C_bits, G_bits, T_bits);

    this->n_nodes = A_bits.size();
    this->k = k;
    this->suffix_group_starts = streaming_support;
    this->colex = colex;
    this->n_kmers = n_kmers;

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
        char c = colex ? kmer[i] : kmer[k-1-i];
        char char_idx = 0;
        if(toupper(c) == 'A') char_idx = 0;
        else if(toupper(c) == 'C') char_idx = 1;
        else if(toupper(c) == 'G') char_idx = 2;
        else if(toupper(c) == 'T') char_idx = 3;
        else return -1; // Invalid character

        node_left = C[char_idx] + subset_rank.rank(node_left, c);
        node_right = C[char_idx] + subset_rank.rank(node_right+1, c) - 1;

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

    // Write number of k-mers
    os.write((char*)&n_kmers, sizeof(n_kmers));
    written += sizeof(n_kmers);

    // Write k
    os.write((char*)&k, sizeof(k));
    written += sizeof(k);

    // Write colex flag
    char flag = colex;
    os.write(&flag, sizeof(flag));
    written += sizeof(flag);

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
    is.read((char*)&n_kmers, sizeof(n_kmers));
    is.read((char*)&k, sizeof(k));

    char colex_flag;
    is.read(&colex_flag, sizeof(colex_flag));
    colex = colex_flag;
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

    // Search the first k-mer
    const char* first_kmer_start = colex ? input : input + len - k;
    ans.push_back(search(first_kmer_start)); 

    for(int64_t i = 1; i < len - k + 1; i++){
        if(ans.back() == -1){
            // Need to search from scratch
            ans.push_back(search(first_kmer_start + (colex ? i : -i)));
        } else{
            // Got to the start of the suffix group and do one search iteration
            int64_t column = ans.back();
            while(suffix_group_starts[column] == 0) column--; // can not go negative because the first column is always marked

            char c = toupper(input[colex ? i+k-1 : len-k-i]);
            char char_idx = -1;
            if(c == 'A') char_idx = 0;
            else if(c == 'C') char_idx = 1;
            else if(c == 'G') char_idx = 2;
            else if(c == 'T') char_idx = 3;
        
            if(char_idx == -1) ans.push_back(-1); // Not found
            else{
                int64_t node_left = column;
                int64_t node_right = column;
                node_left = C[char_idx] + subset_rank.rank(node_left, c);
                node_right = C[char_idx] + subset_rank.rank(node_right+1, c) - 1;
                if(node_left == node_right) ans.push_back(node_left);
                else ans.push_back(-1);
                // Todo: could save one subset rank query if we have fast access to the SBWT columns
            }
        }
    }
    if(!colex) std::reverse(ans.begin(), ans.end());
    return ans;
}

template <typename subset_rank_t>
vector<int64_t> SBWT<subset_rank_t>::streaming_search(const string& input) const{
    return streaming_search(input.c_str(), input.size());
}

} // namespace sbwt