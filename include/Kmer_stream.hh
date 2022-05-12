#pragma once

// Stream of k-mers from a parquet file

#include "arrow/io/file.h"
#include "parquet/stream_reader.h"

#include <iostream>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <fstream>
#include <queue>
#include <stdexcept>

using namespace std;


class Kmer_stream{

    private:

    vector<string> files;
    int64_t k;

    // Stream state
    int64_t current_file_idx;
    std::shared_ptr<arrow::io::ReadableFile> current_infile;
    parquet::StreamReader current_stream;
    string line;
    // End of stream state

    public:

    Kmer_stream(string directory, int64_t k) : current_file_idx(0), k(k){
        const std::filesystem::path path{directory};
        for (auto const& dir_entry : std::filesystem::directory_iterator{path}){
            string filename = dir_entry.path().filename().string();
            if(filename.size() >= 4 && filename.substr(0,4) == "part"){
                files.push_back(dir_entry.path().string());
            }
        }

        std::sort(files.begin(), files.end());

        if(files.size() == 0){
            cerr << "No files with filename starting with 'path' found in the directory" << endl;
            exit(1);
        }

        PARQUET_ASSIGN_OR_THROW(
            current_infile,
            arrow::io::ReadableFile::Open(files[0])
        );

        parquet::StreamReader os{parquet::ParquetFileReader::Open(current_infile)};
        current_stream = std::move(os);
    }

    // Assumes that S is a comma-separated string of tokens, and returns the token with index token_idx
    string get_token_from_line(const string& S, int64_t token_idx){
        vector<int64_t> commas = {-1}; // Start sentinel
        for(int64_t i = 0; i < S.size(); i++){
            if(S[i] == ',') commas.push_back(i);
        }
        commas.push_back(S.size()); // End sentinel
        return S.substr(commas[token_idx] + 1, commas[token_idx + 1] - commas[token_idx] - 1);
    }

    // Puts the next k-mer and its outgoing labels to the given strings. Returns false if there were no more k-mers to read.
    bool next(string& kmer, string& outlabels){
        if(current_file_idx == files.size()) return false;

        while(true){
            while(current_stream.eof()){
                current_file_idx++;
                if(current_file_idx == files.size()) return false; // Done
                PARQUET_ASSIGN_OR_THROW(
                    current_infile,
                    arrow::io::ReadableFile::Open(files[current_file_idx])
                );

                parquet::StreamReader os{parquet::ParquetFileReader::Open(current_infile)};
                current_stream = std::move(os);
            }

            int64_t packed_kmer;
            int32_t packed_edges;
            int32_t length;
            current_stream >> packed_kmer >> packed_edges >> length >> parquet::EndRow;
            kmer = unpack_kmer(packed_kmer, length);
            outlabels = unpack_edges(packed_edges);
            return true;
        }
    }

    // Next k-mer that ends with c. Returns false if not found
    bool next_that_ends_with(string& kmer, string& outlabels, char c){
        while(true){
            bool found = next(kmer, outlabels);
            if(!found) return false;
            if(kmer.back() == c) return true;
        }
    }


    string unpack_edges(int32_t packed_edges){
        string edges;
        if(packed_edges & 1) edges += 'A';
        if(packed_edges & 2) edges += 'C';
        if(packed_edges & 4) edges += 'G';
        if(packed_edges & 8) edges += 'T';
        if(packed_edges & 16) edges += '$';
        return edges;
    }

    // Length can be smaller than k. In that case we left-pad with dollars to length k
    string unpack_kmer(int64_t packed_kmer, int32_t length){
        uint64_t x = *(reinterpret_cast<uint64_t*>(&packed_kmer));
        x >>= (k-length)*2;
        
        string kmer;
        for(int64_t i = 0; i < length; i++){
            if(((x >> (2*i)) & 0x03) == 0x00)
                kmer += 'A';
            if(((x >> (2*i)) & 0x03) == 0x01)
                kmer += 'C';
            if(((x >> (2*i)) & 0x03) == 0x02)
                kmer += 'G';
            if(((x >> (2*i)) & 0x03) == 0x03)
                kmer += 'T';
        }
        return left_pad(kmer, k);
    }

    string left_pad(const string& S, int64_t to_length){
        assert(S.size() <= to_length);
        return string(to_length - S.size(), '$') + S;
    }

};

