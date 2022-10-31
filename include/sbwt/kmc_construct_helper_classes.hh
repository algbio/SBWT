#pragma once

#include "Kmer.hh"
#include <sdsl/bit_vectors.hpp>
#include "buffered_streams.hh"
#include "EM_sort/EM_sort.hh"
#include <set>
#include <unordered_map>
#include <stdexcept>

class CKMCFile; // Defined in KMC
class CKmerAPI; // Defined in KMC

namespace sbwt{

namespace KMC_construction_helper_classes{

typedef long long LL;
typedef Kmer<MAX_KMER_LENGTH> kmer_t;

struct Node{
    kmer_t kmer;
    char edge_flags;

    Node();
    Node(kmer_t kmer);

    static inline int64_t size_in_bytes(){
        return kmer_t::size_in_bytes() + sizeof(char); // char is the edge flags
    }

    void set(char c);
    bool has(char c) const;
    bool operator==(const Node &other) const;
    bool operator!=(const Node &other) const;
    bool operator<(const Node &other) const;
    string to_string() const;
    void serialize(char* buf);
    void load(const char* buf);

};

class Argv{ // Class for turning a vector<string> into char**
private:

    // Forbid copying the class because it wont work right
    Argv(Argv const& other);
    Argv& operator=(Argv const& other);

public:

    char** array = NULL;
    int64_t size = 0;

    Argv(vector<string> v);

    ~Argv();

};

// Also gives reverse complements if asked
class Kmer_stream_from_KMC_DB{

private:

    CKMCFile* kmer_database;
    CKmerAPI* kmer_object;

    uint32_t _kmer_length;
    uint32_t _mode;
    uint32_t _counter_size;
    uint32_t _lut_prefix_length;
    uint32_t _signature_len;
    uint32_t _min_count;
    uint64_t _max_count;
    uint64_t _total_kmers;

    bool add_revcomps;
    std::string str;
    std::string str_revcomp;
    bool revcomp_next = false;

    char get_rc(char c);

    void reverse_complement(string& S);

public:

    Kmer_stream_from_KMC_DB(string KMC_db_path, bool add_revcomps);

    bool done();
    
    Kmer<MAX_KMER_LENGTH> next();

    ~Kmer_stream_from_KMC_DB();
};

// This stream will always start with an empty k-mer with an empty edge label set
class Disk_Instream{

private:

    // Have pointer members -> no copying
    Disk_Instream(Disk_Instream const&) = delete;
    Disk_Instream& operator=(Disk_Instream const&) = delete;

    bool all_read = false;
    Buffered_ifstream<> in;
    char* in_buffer;

    Node top; // Default-initialized to an empty k-mer and an empty edge set

    void update_top();

public:

    Disk_Instream(string filename);
    bool stream_done() const;
    Node stream_next();
    Node peek_next();
    ~Disk_Instream();

};

// A single sorted stream out of two sorted streams
class Node_stream_merger{

    Disk_Instream& A;
    Disk_Instream& B;

public:

    Node_stream_merger(Disk_Instream& A, Disk_Instream& B);

    bool stream_done();

    Node stream_next();

};

} // End of namespace KMC_construction_helper_classes
} // End of namepace sbwt