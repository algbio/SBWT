#pragma once

#include "Kmer.hh"
#include <sdsl/bit_vectors.hpp>
#include "SeqIO/SeqIO.hh"
#include "SeqIO/buffered_streams.hh"
#include "EM_sort/EM_sort.hh"
#include "kmc_construct_helper_classes.hh"
#include "run_kmc.hh"
#include <set>
#include <unordered_map>
#include <stdexcept>

namespace sbwt{

template <typename nodeboss_t>
class NodeBOSSKMCConstructor{

typedef KMC_construction_helper_classes::Disk_Instream Disk_Instream;
typedef KMC_construction_helper_classes::Argv Argv;
typedef KMC_construction_helper_classes::Kmer_stream_from_KMC_DB Kmer_stream_from_KMC_DB;
typedef KMC_construction_helper_classes::kmer_t kmer_t;
typedef KMC_construction_helper_classes::Node Node;
typedef KMC_construction_helper_classes::Node_stream_merger Node_stream_merger;
typedef KMC_construction_helper_classes::SimpleSortedKmerDB SimpleSortedKmerDB;

public:

    // Appends the prefixes of x to nodes
    void add_prefixes(kmer_t z, seq_io::Buffered_ofstream<>& out, char* buf){
        kmer_t prefix = z.copy();
        while(prefix.get_k() > 0){
            char edge_char = prefix.last();
            prefix.dropright();
            Node node(prefix);
            node.set(edge_char);
            node.serialize(buf);
            out.write(buf, Node::size_in_bytes());
        }
    }

    // The result is written to the given sdsl bit vectors. The bits vectors will be resized to fit all the bits.
    void build_bit_vectors_from_sorted_streams(const string& nodefile, const string& dummyfile,
            sdsl::bit_vector& A_bits_sdsl, sdsl::bit_vector& C_bits_sdsl, sdsl::bit_vector& G_bits_sdsl, sdsl::bit_vector& T_bits_sdsl, sdsl::bit_vector& suffix_group_starts_sdsl, int64_t k){
        vector<bool> A_bits, C_bits, G_bits, T_bits, suffix_group_starts;

        // These streams are such that the always start with an empty k-mer and an empty edge set.
        // This will always add the empty string to the graph even if the graph is cyclic. This is
        // intentional to ensure that the root node in the SBWT graph always exists - otherwise the
        // search would need a special case for cyclic graphs.
        Disk_Instream nodes_in(nodefile);
        Disk_Instream dummies_in(dummyfile);

        Node_stream_merger merger(nodes_in, dummies_in);

        Node prev_node;
        bool first = true;
        while(!merger.stream_done()){
            Node x = merger.stream_next();
            if(first || x.kmer != prev_node.kmer){
                // New column
                A_bits.push_back(0);
                C_bits.push_back(0);
                G_bits.push_back(0);
                T_bits.push_back(0);

                // Figure out if this is a suffix group start
                bool is_start = false;
                is_start |= first;
                kmer_t a = prev_node.kmer.copy();
                kmer_t b = x.kmer.copy();
                if(a.get_k() == k) a.dropleft();
                if(b.get_k() == k) b.dropleft();
                is_start |= (a != b);
                suffix_group_starts.push_back(is_start);
                first = false;
            }
            if(x.has('A')) A_bits.back() = 1;
            if(x.has('C')) C_bits.back() = 1;
            if(x.has('G')) G_bits.back() = 1;
            if(x.has('T')) T_bits.back() = 1;
            prev_node = x;
            
        }

        // Copy to sdsl vectors
        A_bits_sdsl.resize(A_bits.size());
        C_bits_sdsl.resize(C_bits.size());
        G_bits_sdsl.resize(G_bits.size());
        T_bits_sdsl.resize(T_bits.size());
        suffix_group_starts_sdsl.resize(suffix_group_starts.size());
        for(int64_t i = 0; i < A_bits.size(); i++){
            A_bits_sdsl[i] = A_bits[i];
            C_bits_sdsl[i] = C_bits[i];
            G_bits_sdsl[i] = G_bits[i];
            T_bits_sdsl[i] = T_bits[i];
            suffix_group_starts_sdsl[i] = suffix_group_starts[i];
        }        
    }

    // This deletes the KMC database on disk after use
    void write_nodes_and_dummies(const string& KMC_db_path, const string& nodes_outfile, const string& dummies_outfile, int64_t n_kmers){
        char node_serialize_buffer[Node::size_in_bytes()];
        
        seq_io::Buffered_ofstream nodes_out(nodes_outfile, ios::binary);
        seq_io::Buffered_ofstream dummies_out(dummies_outfile, ios::binary);

        // These streams are assumed to give k-mers in colex order
        Kmer_stream_from_KMC_DB kmc_db(KMC_db_path, false); // No reverse complements
        write_log("Uncompressing KMC database to disk", LogLevel::MAJOR);
        string uncompressed_db_filename = get_temp_file_manager().create_filename();
        SimpleSortedKmerDB all_stream(kmc_db, uncompressed_db_filename);

        // Delete the KMC database files. The temp file manager can not do this because
        // KMC appends suffixes to the filename and the manager does not know about that.
        std::filesystem::remove(KMC_db_path + ".kmc_pre");
        std::filesystem::remove(KMC_db_path + ".kmc_suf");

        vector<SimpleSortedKmerDB*> char_streams(255); // A k-mer database stream for each character ACGT

        string ACGT = "ACGT";
        for(char c : ACGT){
            char_streams[c] = new SimpleSortedKmerDB(all_stream);
            char_streams[c]->seek_to_char_block(c);
        }

        map<char, Kmer<MAX_KMER_LENGTH>> cur_kmers; // cur_kmers[c] = the colex-smallest k-mer that ends in c that has not been seen yet
        set<char> all_processed; // K-mers ending in these characters have all been processed

        // Initialize first k-mer of each character block
        for(char c : ACGT){
            if(char_streams[c]->done()){
                // This character is not the last character of any k-mer in the data
                all_processed.insert(c);
                break;
            }
            cur_kmers[c] = char_streams[c]->next();
        }
        
        kmer_t prev_x;
        int64_t x_idx = 0;

        write_log("Streaming",LogLevel::MAJOR);
        Progress_printer pp2(n_kmers, 100);
        // Figure out out-edges and which nodes need dummy prefixes
        while(!all_stream.done()){
            kmer_t x = all_stream.next();

            bool suffix_group_start = false;
            if(x_idx == 0 || x.copy().dropleft() != prev_x.copy().dropleft())
                suffix_group_start = true;

            Node x_node(x);            
            if(suffix_group_start){
                for(char c : ACGT){
                    if(all_processed.count(c)) continue;
                    kmer_t y = x.copy().dropleft().appendright(c);
                    kmer_t z = cur_kmers[c];

                    while(y > z){
                        add_prefixes(z,dummies_out,node_serialize_buffer);
                        if(char_streams[c]->done()){
                            all_processed.insert(c);
                            break;
                        } else{
                            z = char_streams[c]->next();
                            cur_kmers[c] = z;
                        }
                    }
                    if(y == z){
                        x_node.set(c);

                        // Advance pointer
                        if(char_streams[c]->done()){
                            all_processed.insert(c);
                            break;
                        } else{
                            cur_kmers[c] = char_streams[c]->next();
                        }
                    }
                }
            }
            x_node.serialize(node_serialize_buffer);
            nodes_out.write(node_serialize_buffer, Node::size_in_bytes());
            x_idx++;
            prev_x = x;
            pp2.job_done();
        }

        // Process the remaining k-mers
        for(char c : ACGT){
            if(all_processed.count(c)) continue;
            while(cur_kmers[c].last() == c){
                add_prefixes(cur_kmers[c], dummies_out, node_serialize_buffer);
                if(char_streams[c]->done()) break;
                else cur_kmers[c] = char_streams[c]->next();
            }
        }

        // Clean up
        for(char c : ACGT) delete char_streams[c];
        get_temp_file_manager().delete_file(uncompressed_db_filename);
    }

    // Construct the given nodeboss from the given input strings
    void build(const vector<string>& input_files, nodeboss_t& nodeboss, int64_t k, int64_t n_threads, int64_t ram_gigas, bool streaming_support, int64_t min_abundance, int64_t max_abundance, int64_t precalc_k){

        string KMC_db_path; int64_t n_kmers;
        std::tie(KMC_db_path, n_kmers) = run_kmc(input_files, k, n_threads, ram_gigas, min_abundance, max_abundance);

        string nodes_outfile = get_temp_file_manager().create_filename();
        string dummies_outfile = get_temp_file_manager().create_filename();

        write_log("Writing nodes and dummies to disk", LogLevel::MAJOR);
        write_nodes_and_dummies(KMC_db_path, nodes_outfile, dummies_outfile, n_kmers);

        write_log("Sorting dummies on disk", LogLevel::MAJOR);
        string dummies_sortedfile = get_temp_file_manager().create_filename();
        EM_sort_constant_binary(dummies_outfile, dummies_sortedfile, 
            [](const char* A, const char* B){
                Node Ax; Ax.load(A);
                Node Bx; Bx.load(B);
                return Ax < Bx;
            }, ram_gigas * ((int64_t)1 <<30), Node::size_in_bytes(), n_threads);
        
        write_log("Merging sorted streams", LogLevel::MAJOR);
        sdsl::bit_vector A_bits, C_bits, G_bits, T_bits, suffix_group_starts;
        build_bit_vectors_from_sorted_streams(nodes_outfile, dummies_sortedfile, A_bits, C_bits, G_bits, T_bits, suffix_group_starts, k);
        
        write_log("Building SBWT structure", LogLevel::MAJOR);
        if(streaming_support){
            nodeboss = nodeboss_t(A_bits, C_bits, G_bits, T_bits, suffix_group_starts, k, n_kmers, precalc_k);
        } else{
            sdsl::bit_vector empty;
            nodeboss = nodeboss_t(A_bits, C_bits, G_bits, T_bits, empty, k, n_kmers, precalc_k);
        }
            
    }
};

}