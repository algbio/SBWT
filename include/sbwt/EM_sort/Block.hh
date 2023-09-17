#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <utility>
#include <tuple>
#include <cstring>
#include <cstdio>
#include <functional>
#include "bit_level_stuff.hh"
#include "SeqIO/buffered_streams.hh"

namespace sbwt{

using namespace std;

class Generic_Block{
    public:
    virtual void sort(const std::function<bool(const char* x, const char* y)>& cmp) = 0;
    virtual void write_to_file(string filename) = 0;
    virtual ~Generic_Block(){}
};

// Takes a pointer to a buffer. Read to the buffer the whole record including the 8-byte size value.
// The buffer might be realloc'd. Returns a bool whether read was successful.
// buffer_len must be greater than 0 (otherwise buffer doubling does not work)
bool read_variable_binary_record(seq_io::Buffered_ifstream<>& input, char** buffer, int64_t* buffer_len);

// Block of records of variable length. The first 8 bytes of a record are a big-endian integer L 
// that tells size of the record. Then follow L-8 bytes which is the "payload" of the record.
class Variable_binary_block : public Generic_Block{
    public:

    char* data; // A buffer. Not necessarily all full of data. Contains concateated records
    int64_t data_len = 0;
    int64_t next_start = 0;
    
    vector<int64_t> starts;

    Variable_binary_block(){
        data = (char*) malloc(1 * sizeof(char));
        data_len = 1;
    }

    ~Variable_binary_block(){
        free(data);
    }

    void double_space(){
        data = (char*)realloc(data, 2 * data_len * sizeof(char));
        data_len *= 2;
    }

    virtual void sort(const std::function<bool(const char* x, const char* y)>& cmp){
        auto cmp_wrap = [&](int64_t x, int64_t y){
            return cmp(data+x,data+y);
        };
        std::sort(starts.begin(), starts.end(), cmp_wrap);
    }

    void add_record(const char* record){
        int64_t space_left = data_len - next_start;
        int64_t rec_len = parse_big_endian_LL(record);
        while(space_left < rec_len){
            double_space();
            space_left = data_len - next_start;
        }

        memcpy(data+next_start, record, rec_len);
        starts.push_back(next_start);
        next_start += rec_len;
    }

    virtual void write_to_file(string filename){
        seq_io::Buffered_ofstream<> out(filename, ios::binary);
        for(int64_t i = 0; i < starts.size(); i++){
            int64_t length = parse_big_endian_LL(data + starts[i]);
            out.write(data + starts[i], length);
        }
    }

    int64_t estimate_size_in_bytes(){
        return next_start + starts.size() * sizeof(int64_t);
    }
};

// RETURN VALUE MUST BE FREED BY CALLER
Variable_binary_block* get_next_variable_binary_block(seq_io::Buffered_ifstream<>& input, int64_t B);

// Records of exactly n bytes each
class Constant_binary_block : public Generic_Block{
public:

    public:

    char* data; // A buffer. Not necessarily all full of data. Contains concateated records
    int64_t data_len = 0;
    int64_t next_start = 0;
    int64_t record_size;
    vector<int64_t> starts; // We sort this vector. The data array is not touched in sorting.

    Constant_binary_block(int64_t record_size) : record_size(record_size) {
        data = (char*) malloc(1 * sizeof(char));
        data_len = 1;
    }

    ~Constant_binary_block(){
        free(data);
    }

    void double_space(){
        data = (char*)realloc(data, 2 * data_len * sizeof(char));
        data_len *= 2;
    }

    virtual void sort(const std::function<bool(const char* x, const char* y)>& cmp){
        auto cmp_wrap = [&](int64_t x, int64_t y){
            return cmp(data+x,data+y);
        };
        std::sort(starts.begin(), starts.end(), cmp_wrap);
    }

    void add_record(const char* record){
        int64_t space_left = data_len - next_start;
        while(space_left < record_size){
            double_space();
            space_left = data_len - next_start;
        }

        memcpy(data+next_start, record, record_size);
        starts.push_back(next_start);
        next_start += record_size;
    }

    virtual void write_to_file(string filename){
        seq_io::Buffered_ofstream<> out(filename, ios::binary);
        for(int64_t i = 0; i < starts.size(); i++){
            out.write(data + starts[i], record_size);
        }
    }

    int64_t estimate_size_in_bytes(){
        return next_start + starts.size() * sizeof(int64_t);
    }

};

// Reads up to B bytes into a new block
// THE RETURN VALUE MUST BE FREED BY THE CALLER
Constant_binary_block* get_next_constant_binary_block(seq_io::Buffered_ifstream<>& input, int64_t B, int64_t record_size);

}