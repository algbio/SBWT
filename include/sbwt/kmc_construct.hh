#pragma once

#include "Kmer.hh"
#include <sdsl/bit_vectors.hpp>
#include "KMC/kmc_api/kmc_file.h"
#include "KMC/include/kmc_runner.h"
#include "KMC_code.hh"
#include "SeqIO.hh"
#include "buffered_streams.hh"
#include "EM_sort/EM_sort.hh"
#include "kmc_construct_helper_classes.hh"
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
typedef KMC_construction_helper_classes::LL LL;
typedef KMC_construction_helper_classes::Node Node;
typedef KMC_construction_helper_classes::Node_stream_merger Node_stream_merger;

public:

    // Appends the prefixes of x to nodes
    void add_prefixes(kmer_t z, Buffered_ofstream<>& out, char* buf){
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
            sdsl::bit_vector& A_bits_sdsl, sdsl::bit_vector& C_bits_sdsl, sdsl::bit_vector& G_bits_sdsl, sdsl::bit_vector& T_bits_sdsl, sdsl::bit_vector& suffix_group_starts_sdsl, LL k){
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
        for(LL i = 0; i < A_bits.size(); i++){
            A_bits_sdsl[i] = A_bits[i];
            C_bits_sdsl[i] = C_bits[i];
            G_bits_sdsl[i] = G_bits[i];
            T_bits_sdsl[i] = T_bits[i];
            suffix_group_starts_sdsl[i] = suffix_group_starts[i];
        }        
    }

    void sort_kmc_db(const string& input_db_file, const string& output_db_file){
        vector<string> args = {"kmc_tools", "transform", input_db_file, "sort", output_db_file};
        Argv argv(args);

        CParametersParser params_parser(argv.size, argv.array);
        params_parser.Parse();
        if (params_parser.validate_input_dbs())
        {
            params_parser.SetThreads();
            CApplication<KMER_WORDS> app(params_parser);
            app.Process();
        }
    }

    // Returns the KMC database prefix and the number of distinct k-mers that had abundance within the given bounds
    pair<string, int64_t> run_kmc(const vector<string>& input_files, LL k, LL n_threads, LL ram_gigas, int64_t min_abundance, int64_t max_abundance){

        write_log("Running KMC counter", LogLevel::MAJOR);

        string KMC_db_file_prefix = get_temp_file_manager().create_filename("kmers");

        KMC::Stage1Params stage1Params;

        string f = input_files[0]; // First input file
        SeqIO::FileFormat format = SeqIO::figure_out_file_format(f);

        for(string f2 : input_files){
            SeqIO::FileFormat format2 = SeqIO::figure_out_file_format(f2);
            if(format.format != format2.format || format.gzipped != format2.gzipped){
                throw std::runtime_error("Error: all input files must have the same format");
            }
        }

        stage1Params.SetInputFiles(input_files)
            .SetKmerLen(k)
            .SetNThreads(n_threads)
            .SetMaxRamGB(ram_gigas)
            .SetInputFileType(format.format == SeqIO::FASTA ? KMC::InputFileType::MULTILINE_FASTA : KMC::InputFileType::FASTQ)
            .SetCanonicalKmers(false)
            .SetTmpPath(get_temp_file_manager().get_dir());

        KMC::Runner kmc;

        auto stage1Results = kmc.RunStage1(stage1Params);

        uint32_t ramForStage2 = ram_gigas;
        KMC::Stage2Params stage2Params;
        stage2Params.SetNThreads(n_threads)
            .SetMaxRamGB(ramForStage2)
            .SetCutoffMin(min_abundance)
            .SetCutoffMax(max_abundance)
            .SetOutputFileName(KMC_db_file_prefix)
            .SetStrictMemoryMode(true);

        auto stage2Results = kmc.RunStage2(stage2Params);

        int64_t n_kmers = stage2Results.nUniqueKmers - stage2Results.nBelowCutoffMin - stage2Results.nAboveCutoffMax;

        write_log("Sorting KMC database", LogLevel::MAJOR);

        try{
            sort_kmc_db(KMC_db_file_prefix, KMC_db_file_prefix + "-sorted");
        } catch(KMCAlreadySortedException& e){
            // Just copy the database
            std::filesystem::copy(KMC_db_file_prefix + ".kmc_pre", KMC_db_file_prefix + "-sorted.kmc_pre");
            std::filesystem::copy(KMC_db_file_prefix + ".kmc_suf", KMC_db_file_prefix + "-sorted.kmc_suf");
        }

        // Delete the unsorted KMC database files. The temp file manager can not do this because
        // KMC appends suffixes to the filename and the manager does not know about that.
        std::filesystem::remove(KMC_db_file_prefix + ".kmc_pre");
        std::filesystem::remove(KMC_db_file_prefix + ".kmc_suf");

        // Clean up the KMC global singleton config state because it seems that it's left
        // in a partial state sometimes, which messes up our code if we call KMC again later.
        CConfig::GetInstance().input_desc.clear();
        CConfig::GetInstance().headers.clear();
        CConfig::GetInstance().simple_output_desc.clear();
        CConfig::GetInstance().transform_output_desc.clear();

        return {KMC_db_file_prefix + "-sorted", n_kmers};
    }

    void write_nodes_and_dummies(const string& KMC_db_path, const string& nodes_outfile, const string& dummies_outfile, LL n_kmers){
        char node_serialize_buffer[Node::size_in_bytes()];
        
        Buffered_ofstream nodes_out(nodes_outfile, ios::binary);
        Buffered_ofstream dummies_out(dummies_outfile, ios::binary);

        // These streams are assumed to give k-mers in colex order
        Kmer_stream_from_KMC_DB all_stream(KMC_db_path, false); // No reverse complements
        map<char, Kmer_stream_from_KMC_DB*> char_streams; // A KMC databse stream for each character

        string ACGT = "ACGT";
        for(char c : ACGT)
            char_streams[c] = new Kmer_stream_from_KMC_DB(KMC_db_path, false); // Kmer stream not be moved -> heap-allocate with new. Remember to delete.

        map<char, Kmer<MAX_KMER_LENGTH>> cur_kmers; // cur_kmers[c] = the colex-smallest k-mer that ends in c that has not been seen yet
        set<char> all_processed; // K-mers ending in these characters have all been processed

        // Rewind character streams to their starting positions
        for(char c : ACGT){
            write_log(string("Rewinding ") + c,LogLevel::MAJOR);
            Progress_printer pp1(n_kmers, 100);
            while(true){
                if(char_streams[c]->done()){
                    // This character is not the last character of any k-mer in the data
                    all_processed.insert(c);
                    break;
                }
                Kmer<MAX_KMER_LENGTH> x = char_streams[c]->next();
                if(x.last() == c){
                    cur_kmers[c] = x;
                    break;
                }
                pp1.job_done();
            }
        }
        
        kmer_t prev_x;
        LL x_idx = 0;

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
    }

    // Construct the given nodeboss from the given input strings
    void build(const vector<string>& input_files, nodeboss_t& nodeboss, LL k, LL n_threads, LL ram_gigas, bool streaming_support, int64_t min_abundance, int64_t max_abundance, LL precalc_k){

        string KMC_db_path; int64_t n_kmers;
        std::tie(KMC_db_path, n_kmers) = run_kmc(input_files, k, n_threads, ram_gigas, min_abundance, max_abundance);

        string nodes_outfile = get_temp_file_manager().create_filename();
        string dummies_outfile = get_temp_file_manager().create_filename();

        write_log("Writing nodes and dummies to disk", LogLevel::MAJOR);
        write_nodes_and_dummies(KMC_db_path, nodes_outfile, dummies_outfile, n_kmers);

        // Delete the KMC database files. The temp file manager can not do this because
        // KMC appends suffixes to the filename and the manager does not know about that.
        std::filesystem::remove(KMC_db_path + ".kmc_pre");
        std::filesystem::remove(KMC_db_path + ".kmc_suf");

        write_log("Sorting dummies on disk", LogLevel::MAJOR);
        string dummies_sortedfile = get_temp_file_manager().create_filename();
        EM_sort_constant_binary(dummies_outfile, dummies_sortedfile, 
            [](const char* A, const char* B){
                Node Ax; Ax.load(A);
                Node Bx; Bx.load(B);
                return Ax < Bx;
            }, ram_gigas * ((LL)1 <<30), Node::size_in_bytes(), n_threads);
        
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