#pragma once

#include "Kmer.hh"
#include <sdsl/bit_vectors.hpp>
#include "KMC/kmc_api/kmc_file.h"
#include "KMC/include/kmc_runner.h"
#include "KMC_code.hh"
#include "input_reading.hh"
#include "buffered_streams.hh"
#include "EM_sort/EM_sort.hh"
#include "kmc_construct_helper_classes.hh"
#include <set>
#include <unordered_map>
#include <stdexcept>

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
    void add_prefixes(kmer_t z, Buffered_ofstream& out, char* buf){
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
            sdsl::bit_vector& A_bits_sdsl, sdsl::bit_vector& C_bits_sdsl, sdsl::bit_vector& G_bits_sdsl, sdsl::bit_vector& T_bits_sdsl){
        vector<bool> A_bits, C_bits, G_bits, T_bits;

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
            }
            if(x.has('A')) A_bits.back() = 1;
            if(x.has('C')) C_bits.back() = 1;
            if(x.has('G')) G_bits.back() = 1;
            if(x.has('T')) T_bits.back() = 1;
            prev_node = x;
            first = false;
        }

        // Copy to sdsl vectors
        A_bits_sdsl.resize(A_bits.size());
        C_bits_sdsl.resize(C_bits.size());
        G_bits_sdsl.resize(G_bits.size());
        T_bits_sdsl.resize(T_bits.size());
        for(LL i = 0; i < A_bits.size(); i++){
            A_bits_sdsl[i] = A_bits[i];
            C_bits_sdsl[i] = C_bits[i];
            G_bits_sdsl[i] = G_bits[i];
            T_bits_sdsl[i] = T_bits[i];
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

    // Returns the KMC database prefix
    string run_kmc(const string& infile, LL k, LL n_threads, LL ram_gigas, int64_t min_abundance){

        write_log("Running KMC counter", LogLevel::MAJOR);

        string KMC_db_file_prefix = get_temp_file_manager().create_filename("kmers");

        std::vector<std::string> inputFiles {infile};
        KMC::Stage1Params stage1Params;

        string file_format = figure_out_file_format(infile);
        if(file_format != "fasta" && file_format != "fastq"){
            throw std::runtime_error("File format not supported: " + file_format);
        }

        stage1Params.SetInputFiles(inputFiles)
            .SetKmerLen(k)
            .SetNThreads(n_threads)
            .SetMaxRamGB(ram_gigas)
            .SetInputFileType(file_format == "fasta" ? KMC::InputFileType::MULTILINE_FASTA : KMC::InputFileType::FASTQ)
            .SetCanonicalKmers(false);

        KMC::Runner kmc;

        auto stage1Results = kmc.RunStage1(stage1Params);

        uint32_t ramForStage2 = ram_gigas;
        KMC::Stage2Params stage2Params;
        stage2Params.SetNThreads(n_threads)
            .SetMaxRamGB(ramForStage2)
            .SetCutoffMin(min_abundance)
            .SetOutputFileName(KMC_db_file_prefix).
            SetStrictMemoryMode(true);

        kmc.RunStage2(stage2Params);

        write_log("Sorting KMC database", LogLevel::MAJOR);

        try{
            sort_kmc_db(KMC_db_file_prefix, KMC_db_file_prefix + "-sorted");
        } catch(KMCAlreadySortedException& e){
            // Just copy the database
            std::filesystem::copy(KMC_db_file_prefix + ".kmc_pre", KMC_db_file_prefix + "-sorted.kmc_pre");
            std::filesystem::copy(KMC_db_file_prefix + ".kmc_suf", KMC_db_file_prefix + "-sorted.kmc_suf");

            // Clean up the KMC global singleton config state because if an exception is thrown
            // the state is left in a bad state
            CConfig::GetInstance().input_desc.clear();
            CConfig::GetInstance().headers.clear();
            CConfig::GetInstance().simple_output_desc.clear();
            CConfig::GetInstance().transform_output_desc.clear();
        }

        // Delete the unsorted KMC database files. The temp file manager can not do this because
        // KMC appends suffixes to the filename and the manager does not know about that.
        std::filesystem::remove(KMC_db_file_prefix + ".kmc_pre");
        std::filesystem::remove(KMC_db_file_prefix + ".kmc_suf");

        return KMC_db_file_prefix + "-sorted";
    }

    void write_nodes_and_dummies(const string& KMC_db_path, const string& nodes_outfile, const string& dummies_outfile){
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
            }
        }
        
        kmer_t prev_x;
        LL x_idx = 0;

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
    void build(const string& infile, nodeboss_t& nodeboss, LL k, LL n_threads, LL ram_gigas, bool streaming_support, int64_t min_abundance){

        string KMC_db_path = run_kmc(infile, k, n_threads, ram_gigas, min_abundance);

        string nodes_outfile = get_temp_file_manager().create_filename();
        string dummies_outfile = get_temp_file_manager().create_filename();

        write_log("Writing nodes and dummies to disk", LogLevel::MAJOR);
        write_nodes_and_dummies(KMC_db_path, nodes_outfile, dummies_outfile);

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
        sdsl::bit_vector A_bits, C_bits, G_bits, T_bits;
        build_bit_vectors_from_sorted_streams(nodes_outfile, dummies_sortedfile, A_bits, C_bits, G_bits, T_bits);
        
        write_log("Building SBWT structure", LogLevel::MAJOR);
        nodeboss.build_from_bit_matrix(A_bits, C_bits, G_bits, T_bits, k, streaming_support);
    }
};