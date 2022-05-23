#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "NodeBOSSInMemoryConstructor.hh"
#include "throwing_streams.hh"
#include "suffix_group_optimization.hh"
#include "libwheeler/BOSS.hh"
#include "globals.hh"
#include "Kmer.hh"
#include <map>

using namespace std;

// Assumes that a root node always exists
template <typename subset_rank_t>
class NodeBOSS{

    public:

    subset_rank_t subset_rank;
    sdsl::bit_vector suffix_group_starts; // Used for streaming queries

    // C-array
    vector<int64_t> C;

    int64_t n_nodes;
    int64_t k;

    NodeBOSS() : n_nodes(0), k(0) {}
    void build_from_strings(const vector<string>& input, int64_t k, bool streaming_support); // This sorts all k-mers in memory and thus takes a lot of memory. Not optimized at all.
    void build_from_WheelerBOSS(const BOSS<sdsl::bit_vector>& boss, bool streaming_support);
    void build_from_bit_matrix(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits, int64_t k, bool streaming_support);
    void build_streaming_query_support(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits); // The NodeBOSS must be already built before calling this
    int64_t search(const string& kmer) const; // Search for std::string
    int64_t search(const char* S) const; // Search for C-string

    // Query for all k-mers in the input
    vector<int64_t> streaming_search(const string& input) const;
    vector<int64_t> streaming_search(const char* input, int64_t len) const;
    bool has_streaming_query_support() const {return suffix_group_starts.size() > 0;}

    // Return the label on the incoming edges to the given node.
    // If the node is the root node, returns a dollar.
    char incoming_label(int64_t node) const; 
    
    int64_t serialize(ostream& out) const; // Returns the number of bytes written
    int64_t serialize(const string& filename) const; // Returns the number of bytes written
    void load(istream& in);
    void load(const string& filename);

};


template <typename subset_rank_t>
void NodeBOSS<subset_rank_t>::build_from_bit_matrix(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits, int64_t k, bool streaming_support){
    subset_rank = subset_rank_t(A_bits, C_bits, G_bits, T_bits);
    n_nodes = A_bits.size();

    // Get the C-array
    C.clear(); C.resize(4);
    C[0] = 1; // There is one incoming ghost-dollar to the root node
    C[1] = C[0] + subset_rank.rank(n_nodes, 'A');
    C[2] = C[1] + subset_rank.rank(n_nodes, 'C');
    C[3] = C[2] + subset_rank.rank(n_nodes, 'G');

    this->k = k;

    if(streaming_support) build_streaming_query_support(A_bits, C_bits, G_bits, T_bits);
}

template <typename subset_rank_t>
void NodeBOSS<subset_rank_t>::build_from_WheelerBOSS(const BOSS<sdsl::bit_vector>& boss, bool streaming_support){

    n_nodes = boss.number_of_nodes();

    // Takes the colex-smallest edge to each node.

    // Build the bit vectors. Later build the bitvector_t versions
    sdsl::bit_vector BWT_A_plain(boss.number_of_nodes(), 0);
    sdsl::bit_vector BWT_C_plain(boss.number_of_nodes(), 0);
    sdsl::bit_vector BWT_G_plain(boss.number_of_nodes(), 0);
    sdsl::bit_vector BWT_T_plain(boss.number_of_nodes(), 0);

    int64_t node_id = -1; // Wheeler rank
    int64_t edge_id = -1; // Wheeler rank
    int64_t edge_count_to_current_node = 0;
    Progress_printer pp(boss.indegs_size(), 10);

    for(int64_t i = 0; i < boss.indegs_size(); i++){
        if(boss.indegs_at(i) == 1){
            node_id++;
            edge_count_to_current_node = 0;
        } else{
            // Process edge
            edge_id++;
            edge_count_to_current_node++;
            char label = boss.incoming_character(node_id);

            if(edge_count_to_current_node == 1){
                int64_t source_node = boss.edge_source(edge_id);
                switch(label){
                    case 'A': BWT_A_plain[source_node] = 1; break;
                    case 'C': BWT_C_plain[source_node] = 1; break;
                    case 'G': BWT_G_plain[source_node] = 1; break;
                    case 'T': BWT_T_plain[source_node] = 1; break;
                    default: cerr << "Invalid character: " << label << "\n"; exit(1);
                }                
            }
        }
        pp.job_done();
    }

    sdsl::bit_vector sg = mark_suffix_groups(BWT_A_plain, BWT_C_plain, BWT_G_plain, BWT_T_plain, boss.get_k());
    push_bits_left(BWT_A_plain, BWT_C_plain, BWT_G_plain, BWT_T_plain, sg);

    build_from_bit_matrix(BWT_A_plain, BWT_C_plain, BWT_G_plain, BWT_T_plain, boss.get_k(), streaming_support);

}

template <typename subset_rank_t>
int64_t NodeBOSS<subset_rank_t>::search(const string& kmer) const{
    assert(kmer.size() == k);
    return search(kmer.c_str());
}

template <typename subset_rank_t>
int64_t NodeBOSS<subset_rank_t>::search(const char* kmer) const{
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
int64_t NodeBOSS<subset_rank_t>::serialize(ostream& os) const{
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
int64_t NodeBOSS<subset_rank_t>::serialize(const string& filename) const{
    throwing_ofstream out(filename, ios::binary);
    return serialize(out.stream);
}


template <typename subset_rank_t>
void NodeBOSS<subset_rank_t>::load(istream& is){
    subset_rank.load(is);
    suffix_group_starts.load(is);
    C = load_std_vector<int64_t>(is);
    is.read((char*)&n_nodes, sizeof(n_nodes));
    is.read((char*)&k, sizeof(k));
}

template <typename subset_rank_t>
void NodeBOSS<subset_rank_t>::load(const string& filename){
    throwing_ifstream in(filename, ios::binary);
    load(in.stream);
}

template <typename subset_rank_t>
char NodeBOSS<subset_rank_t>::incoming_label(int64_t node) const{
    if(node < C[0]) return '$';
    else if(node < C[1]) return 'A';
    else if(node < C[2]) return 'C';
    else if(node < C[3]) return 'G';
    else return 'T';
}


template <typename subset_rank_t>
void NodeBOSS<subset_rank_t>::build_from_strings(const vector<string>& input, int64_t k, bool streaming_support){
    NodeBOSSInMemoryConstructor<NodeBOSS<subset_rank_t>> builder;
    builder.build(input, *this, k, streaming_support);
}

template <typename subset_rank_t>
void NodeBOSS<subset_rank_t>::build_streaming_query_support(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits){
    suffix_group_starts = mark_suffix_groups(A_bits, C_bits, G_bits, T_bits, k);
}

template <typename subset_rank_t>
vector<int64_t> NodeBOSS<subset_rank_t>::streaming_search(const char* input, int64_t len) const{
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
vector<int64_t> NodeBOSS<subset_rank_t>::streaming_search(const string& input) const{
    return streaming_search(input.c_str(), input.size());
}
