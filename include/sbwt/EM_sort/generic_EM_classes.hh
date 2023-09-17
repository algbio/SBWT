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
#include <cassert>
#include <functional>
#include "globals.hh"
#include "Block.hh"
#include "ParallelBoundedQueue.hh"
#include "SeqIO/buffered_streams.hh"

namespace sbwt{

class Generic_Block_Producer{
public:
    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, int64_t block_size) = 0;
    virtual ~Generic_Block_Producer(){}
};

class Generic_Block_Consumer{
public:
    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, const std::function<bool(const char* x, const char* y)>& cmp) = 0;
    virtual vector<string> get_outfilenames() = 0;
    virtual ~Generic_Block_Consumer(){}
};

class Block_Consumer : public Generic_Block_Consumer{ // Works for any type of block
public:

    vector<string> filenames;
    int64_t thread_id;

    Block_Consumer(int64_t thread_id) : thread_id(thread_id) {}

    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, const std::function<bool(const char* x, const char* y)>& cmp){
        while(true){
            Generic_Block* block  = Q.pop();
            if(block == nullptr){
                Q.push(nullptr, 0);
                break; // No more work available
            }
            write_log("Thread " + to_string(thread_id) + ": Starting to sort a block.", LogLevel::MINOR);
            block->sort(cmp);
            write_log("Thread " + to_string(thread_id) + ": Finished sorting a block.", LogLevel::MINOR);
            filenames.push_back(get_temp_file_manager().create_filename());
            block->write_to_file(filenames.back());
            delete block; // Initially allocated by a producer
        }
    }

    virtual vector<string> get_outfilenames(){
        return filenames;
    }
};

//
// Constant record classes below
//

class Constant_Block_Producer : public Generic_Block_Producer{
    public:

    seq_io::Buffered_ifstream<> in;
    int64_t record_size;

    Constant_Block_Producer(string infile, int64_t record_size) : in(infile, ios::binary), record_size(record_size) {}

    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, int64_t block_size){
        while(true){
            Constant_binary_block* block = get_next_constant_binary_block(in, block_size,record_size); // Freed by a consumer
            if(block->starts.size() == 0){
                Q.push(nullptr, 0);
                break; // No more work available
            }
            Q.push((Generic_Block*)block, block_size);
        }
    }
};

class Constant_Record_Reader{
public:

    vector<seq_io::Buffered_ifstream<>> inputs;
    int64_t record_size;

    Constant_Record_Reader(int64_t record_size) : record_size(record_size){}

    void open_files(vector<string> filenames){
        inputs.clear();
        inputs.resize(filenames.size());
        for(int64_t i = 0; i < filenames.size(); i++){
            inputs[i].open(filenames[i], ios::binary);
        }
    }

    void close_files(){
        for(seq_io::Buffered_ifstream<>& in : inputs) in.close();
    }

    int64_t get_num_files(){
        return inputs.size();
    }

    bool read_record(int64_t input_index, char** buffer, int64_t* buffer_size){
        if(*buffer_size < record_size){
            *buffer = (char*)realloc(*buffer, record_size);
            *buffer_size = record_size;
        }
        return inputs[input_index].read(*buffer, record_size);
    }

};

class Constant_Record_Writer{
public:
    seq_io::Buffered_ofstream<> out;
    int64_t record_size;

    Constant_Record_Writer(int64_t record_size) : record_size(record_size){}

    void open_file(string filename){
        out.open(filename, ios::binary);
    }

    void close_file(){
        out.close();
    }

    void write(char* record){
        out.write(record, record_size);
    }

};

//
// Variable binary record classes below
//

class Variable_Block_Producer : public Generic_Block_Producer{
    public:

    seq_io::Buffered_ifstream<> in;

    Variable_Block_Producer(string infile) : in(infile, ios::binary) {}

    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, int64_t block_size){
        while(true){
            Variable_binary_block* block = get_next_variable_binary_block(in,block_size); // Freed by a consumer
            if(block->starts.size() == 0){
                Q.push(nullptr, 0);
                break; // No more work available
            }
            Q.push((Generic_Block*)block, block_size);
        }
    }
};

class Variable_Record_Reader{
public:

    vector<seq_io::Buffered_ifstream<>> inputs;

    Variable_Record_Reader() {}

    void open_files(vector<string> filenames){
        inputs.clear();
        inputs.resize(filenames.size());
        for(int64_t i = 0; i < filenames.size(); i++){
            inputs[i].open(filenames[i], ios::binary);
        }
    }

    void close_files(){
        for(seq_io::Buffered_ifstream<>& in : inputs) in.close();
    }

    int64_t get_num_files(){
        return inputs.size();
    }

    bool read_record(int64_t input_index, char** buffer, int64_t* buffer_size){
        return read_variable_binary_record(inputs[input_index], buffer, buffer_size);
    }

};

class Variable_Record_Writer{
public:
    seq_io::Buffered_ofstream<> out;

    Variable_Record_Writer() {}

    void open_file(string filename){
        out.open(filename, ios::binary);
    }

    void close_file(){
        out.close();
    }

    void write(char* record){
        int64_t rec_len = parse_big_endian_LL(record);
        out.write(record, rec_len);
    }

};

}