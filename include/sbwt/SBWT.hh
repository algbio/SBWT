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

const std::string SBWT_VERSION = "v0.1"; // Update this after breaking changes. This is serialized with the index and checked when loading.

// Assumes that a root node always exists
template <typename subset_rank_t>
class SBWT{

private:

    subset_rank_t subset_rank; // The subset rank query implementation
    sdsl::bit_vector suffix_group_starts; // Marks the first column of every suffix group (see paper)
    vector<int64_t> C; // The array of cumulative character counts

    vector<pair<int64_t,int64_t> > kmer_prefix_precalc; // SBWT intervals for all p-mers with p = precalc_k.
    int64_t precalc_k = 0;

    int64_t n_nodes; // Number of nodes (= columns) in the data structure
    int64_t n_kmers; // Number of k-mers indexed in the data structure
    int64_t k; // The k-mer k

    int64_t get_char_idx(char c) const{
        switch(c){
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return -1;
        }
    }

public:

    struct BuildConfig{
        vector<string> input_files; /**< List of paths to input filenames. */
        int k = 30; /**< Length of the k-mers. */
        bool build_streaming_support = true; /**< Whether we should build the streaming query support. */
        int n_threads = 1; /**< Number of parallel threads in construction. */
        int min_abundance = 1; /**< k-mers occurring fewer than this many times are discarded. */
        int max_abundance = 1e9; /**< k-mers occurring more than this many times are discarded */
        int ram_gigas = 2; /**< RAM budget in gigabytes. Not strictly enforced. */
        int precalc_k = 0; /**< We will precalculate and store the SBWT intervals of all DNA-strings of this length */
        string temp_dir = "."; /**< Path to the directory for the temporary files. */
    };

    /**
     * @brief Construct an empty SBWT.
     * 
     */
    SBWT() : n_nodes(0), n_kmers(0), k(0){}

    /**
     * @brief Construct SBWT from plain matrix SBWT bit vectors.
     * 
     * @param A_bits Row of character A in the plain matrix SBWT.
     * @param C_bits Row of character C in the plain matrix SBWT.
     * @param G_bits Row of character G in the plain matrix SBWT.
     * @param T_bits Row of character T in the plain matrix SBWT.
     * @param streaming_support The streaming support bit vector. Can be empty.
     * @param k Length of the k-mers.
     * @param number_of_kmers Number of k-mers in the data structure.
     */
    SBWT(const sdsl::bit_vector& A_bits, 
         const sdsl::bit_vector& C_bits, 
         const sdsl::bit_vector& G_bits, 
         const sdsl::bit_vector& T_bits, 
         const sdsl::bit_vector& streaming_support, // Streaming support may be empty
         int64_t k, 
         int64_t number_of_kmers,
         int64_t precalc_k);

    /**
     * @brief Construct SBWT using the KMC-based construction algorithm.
     * 
     * @param config construction paramters.
     */
    SBWT(const BuildConfig& config); 

    // Accessors

    /**
     * @brief Get a const reference to the internal subset rank structure.
     */
    const subset_rank_t& get_subset_rank_structure() const {return subset_rank;}

    /**
     * @brief Get a const reference to the internal streaming support bit vector.
     */
    const sdsl::bit_vector& get_streaming_support() const {return suffix_group_starts;}

    /**
     * @brief Compute and return a bit vector that marks which nodes do not correspond to a full k-mer.
     */
    sdsl::bit_vector compute_dummy_node_marks() const;

    /**
     * @brief Get a const reference to the cumulative character count array.
     */
    const vector<int64_t>& get_C_array() const {return C;}

    /**
     * @brief Get a const reference to the k-mer prefix precalc
     */
    const vector<pair<int64_t, int64_t>>& get_precalc() const {return kmer_prefix_precalc;}

    /**
     * @brief Get the precalc k-mer prefix length
     */
    int64_t get_precalc_k() const {return precalc_k;}


    /**
     * @brief Precalculate all SBWT intervals of all strings of length prefix_length. These will be used in search.
     */
    void do_kmer_prefix_precalc(int64_t prefix_length);

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
     * @brief Follow an edge in the de Bruijn graph. If called on a dummy node, follows an edge in the dummy node tree.
     * 
     * @param node The node to move from.
     * @param c The character to follow.
     * @return The node id of the node at the end of the edge from `node` labeled with `c`, or -1 if does not exist.
     */
    int64_t forward(int64_t node, char c) const;

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
     * @brief Run SBWT search iterations for characters in S starting from interval I.
     * 
     * @param S The string to search for. Can be of any length.
     * @param I The SBWT interval to start from.
     * @return The updated interval, or {-1,-1} if it was not found.
     */
    std::pair<int64_t,int64_t> update_sbwt_interval(const string& S, std::pair<int64_t,int64_t> I) const;

    /**
     * @brief Run SBWT search iterations for characters in S starting from interval I.
     * 
     * @param S The string to search for.
     * @param S_length The length of S.
     * @param I The SBWT interval to start from.
     * @return The updated interval, or {-1,-1} if it was not found.
     */
    std::pair<int64_t,int64_t> update_sbwt_interval(const char* S, int64_t S_length, std::pair<int64_t,int64_t> I) const;

    /**
     * @brief Query all k-mers of the input std::string. Requires that the streaming support had been built.
     * 
     * @throws std::runtime_error If the streaming support has not been built.
     * @param input The input string 
     * @return vector<int64_t> The ranks of the k-mers of the input in the data structure, with -1 for those that are not found in the index. The result will be the same as if search() was called for each k-mer of the input from left to right in order.
     * @see search()
     */
    vector<int64_t> streaming_search(const string& input) const;

    /**
     * @brief Query all k-mers of the input C-string. Requires that the streaming support had been built.
     * 
     * @throws std::runtime_error If the streaming support has not been built.
     * @param input The input string 
     * @param len Length of the input string
     * @return vector<int64_t> The ranks of the k-mers of the input in the data structure, with -1 for those that are not found in the index. The result will be the same as if search() was called for each k-mer of the input from left to right in order.
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
SBWT<subset_rank_t>::SBWT(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits, const sdsl::bit_vector& streaming_support, int64_t k, int64_t n_kmers, int64_t precalc_k){
    subset_rank = subset_rank_t(A_bits, C_bits, G_bits, T_bits);

    this->n_nodes = A_bits.size();
    this->k = k;
    this->suffix_group_starts = streaming_support;
    this->n_kmers = n_kmers;

    // Get the C-array
    C.clear(); C.resize(4);
    C[0] = 1; // There is one incoming ghost-dollar to the root node
    C[1] = C[0] + subset_rank.rank(n_nodes, 'A');
    C[2] = C[1] + subset_rank.rank(n_nodes, 'C');
    C[3] = C[2] + subset_rank.rank(n_nodes, 'G');

    do_kmer_prefix_precalc(precalc_k);

}

template <typename subset_rank_t>
SBWT<subset_rank_t>::SBWT(const BuildConfig& config){
    string old_temp_dir = get_temp_file_manager().get_dir();
    get_temp_file_manager().set_dir(config.temp_dir);

    NodeBOSSKMCConstructor<SBWT<subset_rank_t>> builder;
    builder.build(config.input_files, *this, config.k, config.n_threads, config.ram_gigas, config.build_streaming_support, config.min_abundance, config.max_abundance, config.precalc_k);

    get_temp_file_manager().set_dir(old_temp_dir); // Return the old temporary directory

    // Precalc will be done in the other constructor which is called by the builder
}

template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::forward(int64_t node, char c) const{
    if(!has_streaming_query_support())
        throw std::runtime_error("Error: Streaming support required for SBWT::forward");

    // Go to start of the suffix group.
    while(!suffix_group_starts[node]) node--; // Guaranteed to terminate because the first node is always marked

    int64_t r1 = subset_rank.rank(node, c);
    int64_t r2 = subset_rank.rank(node+1, c);
    if(r1 == r2) return -1; // No edge found. TODO: could save one rank query if we had direct access to the SBWT sets

    return C[get_char_idx(c)] + r1;
}

template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::search(const string& kmer) const{
    assert(kmer.size() == k);
    return search(kmer.c_str());
}

template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::search(const char* kmer) const{
    pair<int64_t, int64_t> I;
    
    if(precalc_k > 0){ // Precalc is available
        // Find the index of the k-mer prefix in the precalc table (see do_kmer_prefix_precalc).
        // Todo: do this in a more bit-parallel way
        uint64_t precalc_idx = 0;
        for(int64_t i = 0; i < precalc_k; i++){
            int64_t char_idx = DNA_to_char_idx(kmer[precalc_k-1-i]);
            if(char_idx == -1) return -1; // non-ACGT character
            precalc_idx = (precalc_idx << 2) | char_idx; // Add the character
        }

        // Continue search from precalculated interval
        I = update_sbwt_interval(kmer + precalc_k, k - precalc_k, kmer_prefix_precalc[precalc_idx]);
    }

    else // No precalc
        I = update_sbwt_interval(kmer, k, {0, n_nodes-1});

    if(I.first != I.second){
        cerr << "Bug: k-mer search did not give a singleton interval: " << I.first << " " << I.second << endl;
        exit(1);
    }
    return I.first;
}

template<typename subset_rank_t>
std::pair<int64_t,int64_t> SBWT<subset_rank_t>::update_sbwt_interval(const string& S, pair<int64_t,int64_t> I) const{
    return update_sbwt_interval(S.c_str(), S.size(), I);
}

template<typename subset_rank_t>
std::pair<int64_t,int64_t> SBWT<subset_rank_t>::update_sbwt_interval(const char* S, int64_t S_length, pair<int64_t,int64_t> I) const{
    if(I.first == -1) return I;
    for(int64_t i = 0; i < S_length; i++){
        char c = toupper(S[i]);
        int64_t char_idx = get_char_idx(S[i]);
        if(char_idx == -1) return {-1,-1}; // Invalid character

        I.first = C[char_idx] + subset_rank.rank(I.first, c);
        I.second = C[char_idx] + subset_rank.rank(I.second+1, c) - 1;

        if(I.first > I.second) return {-1,-1}; // Not found
    }

    return {I.first, I.second};
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

    written += serialize_string(SBWT_VERSION, os);

    written += subset_rank.serialize(os);
    written += suffix_group_starts.serialize(os);

    written += serialize_std_vector(C, os);

    // Write precalc
    written += serialize_std_vector(kmer_prefix_precalc, os);
    os.write((char*)&precalc_k, sizeof(precalc_k));
    written += sizeof(precalc_k);

    // Write number of nodes
    os.write((char*)&n_nodes, sizeof(n_nodes));
    written += sizeof(n_nodes);

    // Write number of k-mers
    os.write((char*)&n_kmers, sizeof(n_kmers));
    written += sizeof(n_kmers);

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
    string version = load_string(is);
    if(version != SBWT_VERSION){
        throw std::runtime_error("Error: Corrupt index file, or the index was constructed with an incompatible version of SBWT.");
    }

    subset_rank.load(is);
    suffix_group_starts.load(is);
    C = load_std_vector<int64_t>(is);
    kmer_prefix_precalc = load_std_vector<pair<int64_t, int64_t>>(is);
    is.read((char*)&precalc_k, sizeof(precalc_k));
    is.read((char*)&n_nodes, sizeof(n_nodes));
    is.read((char*)&n_kmers, sizeof(n_kmers));
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

    // Search the first k-mer
    const char* first_kmer_start = input;
    ans.push_back(search(first_kmer_start)); 

    for(int64_t i = 1; i < len - k + 1; i++){
        if(ans.back() == -1){
            // Need to search from scratch
            ans.push_back(search(first_kmer_start + i));
        } else{
            // Got to the start of the suffix group and do one search iteration
            int64_t column = ans.back();
            while(suffix_group_starts[column] == 0) column--; // can not go negative because the first column is always marked

            char c = toupper(input[i+k-1]);
            int64_t char_idx = get_char_idx(c);
        
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
    return ans;
}

template <typename subset_rank_t>
vector<int64_t> SBWT<subset_rank_t>::streaming_search(const string& input) const{
    return streaming_search(input.c_str(), input.size());
}

template <typename subset_rank_t>
sdsl::bit_vector SBWT<subset_rank_t>::compute_dummy_node_marks() const{
    int64_t count = 0;
    vector<pair<int64_t,int64_t>> dfs_stack; // pairs (node, depth)
    dfs_stack.push_back({0, 0}); // Root node
    // dfs to depth k-1
    // the dummy part is a tree so no visited-list is required

    string ACGT = "ACGT";
    sdsl::bit_vector marks(n_nodes, 0);
    int64_t v,d; // node,depth
    while(!dfs_stack.empty()){
        tie(v,d) = dfs_stack.back();
        dfs_stack.pop_back();
        if(d < k){
            count++;
            marks[v] = 1;
        }
        if(d < k-1){ // Push children
            for(char c : ACGT){
                int64_t u = forward(v, c);
                if(u != -1) dfs_stack.push_back({u,d+1});
            }
        }
    }
    return marks;
}

template <typename subset_rank_t>
void SBWT<subset_rank_t>::do_kmer_prefix_precalc(int64_t prefix_length){
    if(prefix_length == 0) return;
    if(prefix_length > 20){
        throw std::runtime_error("Error: Can't precalc longer than 20-mers (would take over 4^20 = 2^40 bytes");
    }

    if(prefix_length > k)
        throw std::runtime_error("Error: Precalc length is longer than k (" + to_string(prefix_length) + " > " + to_string(k) + ")");
    
    uint64_t n_kmers_to_precalc = 1 << (2*prefix_length); // Four to the power prefix_length

    // Initialize member variables
    kmer_prefix_precalc.resize(n_kmers_to_precalc);
    this->precalc_k = prefix_length;

    uint64_t data = 0; // K-mer packed 2 bits per nucleotide
    string prefix(prefix_length, '\0');

    while(n_kmers_to_precalc--){
        for(int64_t i = 0; i < prefix_length; i++){
            char c = char_idx_to_DNA((data >> (2*i)) & 0x3); // Decode the i-th character
            prefix[i] = c;
        }

        kmer_prefix_precalc[data] = update_sbwt_interval(prefix, {0, n_nodes-1});
        data++;
    }

}

} // namespace sbwt