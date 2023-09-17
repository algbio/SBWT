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
#include <gtest/gtest.h>
#include "globals.hh"
#include "EM_sort/EM_sort.hh"
#include "setup_tests.hh"

using namespace sbwt;

string generate_variable_binary_testcase(int64_t max_record_len_bytes, int64_t n_records){
    string outfile = get_temp_file_manager().create_filename();
    seq_io::Buffered_ofstream out(outfile, ios::binary);
    for(int64_t i = 0; i < n_records; i++){
        int64_t record_len = max((int64_t)8, rand() % (max_record_len_bytes + 1));
        //cout << "Add record of length " << record_len << endl;
        write_big_endian_LL(out, record_len);
        for(int64_t b = 0; b < record_len - 8; b++){
            char byte = static_cast<char>(rand() % 256);
            out.write(&byte, 1);
        }
    }
    return outfile;
}

string generate_constant_binary_testcase(int64_t record_len, int64_t n_records){
    string outfile = get_temp_file_manager().create_filename();
    seq_io::Buffered_ofstream out(outfile, ios::binary);
    for(int64_t i = 0; i < n_records; i++){
        for(int64_t b = 0; b < record_len; b++){
            char byte = static_cast<char>(rand() % 256);
            out.write(&byte, 1);
        }
    }
    return outfile;
}

string record_to_string(const char* rec){
    stringstream ss;
    int64_t len = parse_big_endian_LL(rec);
    ss << len << ": ";
    for(int64_t i = 0; i < len-8; i++){
        int value = static_cast<int>(*reinterpret_cast<const unsigned char*>(rec+8+i));
        ss << value << " ";
    }
    return ss.str();
}

string binary_sort_stdlib(string infile, const std::function<bool(const char*, const char*)> & cmp){
    
    // Using the Block class to help reading in the records
    seq_io::Buffered_ifstream in(infile);
    Variable_binary_block* block = get_next_variable_binary_block(in,(int64_t)1e16);
    auto cmp_wrap = [&](int64_t x, int64_t y){
        return cmp(block->data+x,block->data+y);
    };
    std::sort(block->starts.begin(), block->starts.end(), cmp_wrap);

    // Verify that the sort worked
    for(int64_t i = 1; i < block->starts.size(); i++){
        // record i must not be strictly less than record i-1
        EXPECT_FALSE(cmp(block->data + block->starts[i], block->data + block->starts[i-1]));
    }

    string outfile = get_temp_file_manager().create_filename();
    seq_io::Buffered_ofstream out(outfile, ios::binary);
    for(int64_t i = 0; i < block->starts.size(); i++){
        int64_t length = parse_big_endian_LL(block->data + block->starts[i]);
        out.write(block->data + block->starts[i], length);
    }
    delete block;
    return outfile;
}

string constant_binary_sort_stdlib(string infile, int64_t rec_len, const std::function<bool(const char*, const char*)> & cmp){
    
    // Using the Block class to help reading in the records
    seq_io::Buffered_ifstream in(infile);
    Constant_binary_block* block = get_next_constant_binary_block(in,(int64_t)1e16, rec_len);
    auto cmp_wrap = [&](int64_t x, int64_t y){
        return cmp(block->data+x,block->data+y);
    };
    std::sort(block->starts.begin(), block->starts.end(), cmp_wrap);

    // Verify that the sort worked
    for(int64_t i = 1; i < block->starts.size(); i++){
        // record i must not be strictly less than record i-1
        EXPECT_FALSE(cmp(block->data + block->starts[i], block->data + block->starts[i-1]));
    }

    string outfile = get_temp_file_manager().create_filename();
    seq_io::Buffered_ofstream out(outfile, ios::binary);
    for(int64_t i = 0; i < block->starts.size(); i++){
        out.write(block->data + block->starts[i], rec_len);
    }
    delete block;
    return outfile;
}

void test_variable_binary_sort(string infile, const std::function<bool(const char*, const char*)> & cmp){
    
    string fileA = binary_sort_stdlib(infile, cmp);
    string fileB = get_temp_file_manager().create_filename();
    int64_t ram = rand() % 10000 + 1;
    logger << "Sorting variable binary records with " << ram << " RAM" << endl;
    EM_sort_variable_length_records(infile, fileB, cmp, ram, 3);
    ASSERT_TRUE(files_are_equal(fileA, fileB));

    get_temp_file_manager().delete_file(fileA);
    get_temp_file_manager().delete_file(fileB);
    
}

void test_constant_binary_sort(string infile, int64_t record_len, const std::function<bool(const char*, const char*)> & cmp){
    
    string fileA = constant_binary_sort_stdlib(infile, record_len, cmp);
    string fileB = get_temp_file_manager().create_filename();
    int64_t ram = rand() % 10000 + 1;
    logger << "Sorting constant-binary records with " << ram << " RAM" << endl;
    EM_sort_constant_binary(infile, fileB, cmp, ram, record_len, 3);
    ASSERT_TRUE(files_are_equal(fileA, fileB));

    get_temp_file_manager().delete_file(fileA);
    get_temp_file_manager().delete_file(fileB);
    
}

TEST(TEST_EM_SORT, variable_binary_sort){
    for(int64_t max_record_len_bytes = 8; max_record_len_bytes <= 1e6; max_record_len_bytes *= 2){
        // Max record length must be at least 8 because of the long long at the start that tells the length
        for(int64_t n_records = 1; n_records <= 1e6 && max_record_len_bytes*n_records <= 1e6; n_records *= 2){
            logger << max_record_len_bytes << " " << n_records << endl;
            string infile = generate_variable_binary_testcase(max_record_len_bytes, n_records);
            test_variable_binary_sort(infile, memcmp_variable_binary_records);
        }
    }
}

TEST(TEST_EM_SORT, constant_binary_sort){
    
    for(int64_t record_len = 1; record_len <= 1e6; record_len *= 2){
        for(int64_t n_records = 1; n_records <= 1e6 && record_len*n_records <= 1e6; n_records *= 2){
            logger << record_len << " " << n_records << endl;

            auto cmp = [&](const char* x, const char* y){
                return memcmp(x,y,record_len) < 0;
            };
            
            string infile = generate_constant_binary_testcase(record_len, n_records);
            test_constant_binary_sort(infile, record_len, cmp);
        }
    }
}