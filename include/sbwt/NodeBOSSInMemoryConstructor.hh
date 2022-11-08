#pragma once

#include <string>
#include <algorithm>
#include <unordered_set>
#include "Kmer.hh"
#include "sdsl/bit_vectors.hpp"

namespace sbwt{

/*
This file implements an in-memory construction algorithm for NodeBOSS. It is not intended for
production use, but rather as a reference to debug other construction algorithm.
*/

template <typename nodeboss_t>
class NodeBOSSInMemoryConstructor{

    typedef Kmer<MAX_KMER_LENGTH> kmer_t;
    

    public:

    struct Node{
        kmer_t kmer;
        char edge_flags;

        Node() : kmer(), edge_flags(0) {}
        Node(kmer_t kmer) : kmer(kmer), edge_flags(0) {}

        void set(char c){
            if(c == 'A') edge_flags |= 1 << 0;
            else if(c == 'C') edge_flags |= 1 << 1;
            else if(c == 'G') edge_flags |= 1 << 2;
            else if(c == 'T') edge_flags |= 1 << 3;
        }

        bool has(char c){
            if(c == 'A') return edge_flags & (1 << 0);
            else if(c == 'C') return edge_flags & (1 << 1);
            else if(c == 'G') return edge_flags & (1 << 2);
            else if(c == 'T') return edge_flags & (1 << 3);
            return false;
        }

        bool operator==(const Node &other) const{
            return this->kmer == other.kmer && this->edge_flags == other.edge_flags;
        }

        bool operator!=(const Node &other) const{
            return !(*this == other);
        }

        bool operator<(const Node &other) const{
            if(this->kmer < other.kmer) return true;
            if(this->kmer == other.kmer && this->edge_flags < other.edge_flags) return true;
            return false;
        }

    };

    int64_t get_char_ptr(const vector<kmer_t>& kmers, char c){
        for(int64_t i = 0; i < kmers.size(); i++){
            if(kmers[i].last() == c) return i;
        }
        return kmers.size(); // Not found -> return one-past-the-end
    }

    // Appends the prefixes of x to nodes
    void add_prefixes(kmer_t z, vector<Node>& nodes){
        kmer_t prefix = z.copy();
        while(prefix.get_k() > 0){
            char edge_char = prefix.last();
            prefix.dropright();
            Node node(prefix);
            node.set(edge_char);
            nodes.push_back(node);
        }
    }

    // Input must be sorted
    void merge_equal_nodes(vector<Node>& nodes){
        string ACGT = "ACGT";
        int64_t j = 0; // Next index in the output array (which is the same as the input array as we are doing this in-place).
        for(int64_t i = 0; i < nodes.size(); i++){
            if(i > 0 && nodes[i].kmer == nodes[i-1].kmer){
                // Copy edges to the latest node in the output
                for(char c : ACGT) if(nodes[i].has(c)) nodes[j-1].set(c);
            } else{
                nodes[j++] = nodes[i];
            }
        }
        nodes.resize(j);
    }

    // kmers must be unique and be colex-sorted before calling.
    // Returns node list in colex order.
    vector<Node> get_nodes(const vector<kmer_t>& kmers){
        string ACGT = "ACGT";
        map<char, int64_t> char_ptrs;
        for(char c : ACGT) char_ptrs[c] = get_char_ptr(kmers, c);

        vector<Node> nodes;
        Node empty;
        nodes.push_back(empty); // Always have a root node.

        for(int64_t i = 0; i < kmers.size(); i++){
            bool suffix_group_start = false;
            if(i == 0 || kmers[i].copy().dropleft() != kmers[i-1].copy().dropleft())
                suffix_group_start = true;

            kmer_t x = kmers[i];
            Node x_node(x);
            if(suffix_group_start){
                // Try all edges
                for(char c : ACGT){
                    if(char_ptrs[c] == kmers.size()) continue; // This character is done
                    // Let y = x[1..k-1] c
                    kmer_t y = kmers[i].copy().dropleft().appendright(c);

                    kmer_t z = kmers[char_ptrs[c]];
                    while(y > z){
                        // z has no incoming edge -> need dummy nodes
                        add_prefixes(z, nodes);
                        char_ptrs[c]++;
                        if(char_ptrs[c] == kmers.size()) break;
                        z = kmers[char_ptrs[c]];
                    }
                    if(y == z){ // Found. Record the edge label to x_node
                        char_ptrs[c]++;
                        x_node.set(c);
                    } else if(y < z){
                        // This edge does not exist (no problem)
                    }
                        
                }
            }
            nodes.push_back(x_node);
        }

        // Process remaining nodes
        for(char c : ACGT){
            while(char_ptrs[c] < kmers.size() && kmers[char_ptrs[c]].last() == c){
                add_prefixes(kmers[char_ptrs[c]++], nodes);
            }
        }

        std::sort(nodes.begin(), nodes.end());

        merge_equal_nodes(nodes);

        return nodes;
 
    }

    bool is_valid_kmer(const string& S){
        for(char c : S) if(c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
        return true;
    }

    vector<kmer_t> get_distinct_kmers(const vector<string>& input, int64_t k){
        write_log("Hashing distinct k-mers", LogLevel::MAJOR);
        unordered_set<kmer_t> kmer_ht; // k-mer hash table
        for(const string& S : input){
            for(int64_t i = 0; i < (int64_t)S.size()-k+1; i++){
                if(is_valid_kmer(S.substr(i,k)))
                    kmer_ht.insert(kmer_t(S.substr(i,k)));
            }
        }
        vector<kmer_t> kmers(kmer_ht.begin(), kmer_ht.end());
        return kmers;
    }

    sdsl::bit_vector build_streaming_support(vector<Node>& nodes, int64_t k){
        sdsl::bit_vector bv(nodes.size(), 0);
        bv[0] = 1;
        for(int64_t i = 1; i < nodes.size(); i++){
            kmer_t A = nodes[i-1].kmer;
            kmer_t B = nodes[i].kmer;
            if(A.get_k() == k) A.dropleft();
            if(B.get_k() == k) B.dropleft();
            bv[i] = (A != B);
        }
        return bv;
    }

    // Construct the given nodeboss from the given input strings
    void build(const vector<string>& input, nodeboss_t& nodeboss, int64_t k, bool streaming_support){

        vector<kmer_t> kmers = get_distinct_kmers(input, k);
        std::sort(kmers.begin(), kmers.end());

        write_log("Sorting nodes", LogLevel::MAJOR);
        vector<Node> nodes = get_nodes(kmers);

        write_log("Building SBWT", LogLevel::MAJOR);

        sdsl::bit_vector A_bits(nodes.size(), 0);
        sdsl::bit_vector C_bits(nodes.size(), 0);
        sdsl::bit_vector G_bits(nodes.size(), 0);
        sdsl::bit_vector T_bits(nodes.size(), 0);
        for(int64_t i = 0; i < nodes.size(); i++){
            if(nodes[i].has('A')) A_bits[i] = 1;
            if(nodes[i].has('C')) C_bits[i] = 1;
            if(nodes[i].has('G')) G_bits[i] = 1;
            if(nodes[i].has('T')) T_bits[i] = 1;
        }

        sdsl::bit_vector ssupport;
        if(streaming_support) ssupport = build_streaming_support(nodes, k);
        nodeboss_t constructed(A_bits, C_bits, G_bits, T_bits, ssupport, k, kmers.size(), 0); // No precalc
        nodeboss = constructed;
    }
};

template<typename nodeboss_t>
void build_nodeboss_in_memory(const vector<string>& input, nodeboss_t& nodeboss, int64_t k, bool streaming_support){
    NodeBOSSInMemoryConstructor<nodeboss_t> builder;
    builder.build(input, nodeboss, k, streaming_support);
}

}