#pragma once

#include <gtest/gtest.h>
#include "globals.hh"
#include "setup_tests.hh"
#include "SeqIO.hh"

using namespace sbwt;

TEST(TEST_BUFFERED_IO, getline){
    string filename = get_temp_file_manager().create_filename();
    throwing_ofstream out(filename, ios::binary);
    out << "\n111\n\nA\n2222222"; 
    char zero = 0;
    out.write(&zero, 1); // Null in the middle of a string
    out << "------\n\n333333\n4444"; // Does not end in newline
    out.close();

    // Check that behaviour is identical std::ifstream

    vector<string> correct_lines;
    ifstream in(filename, ios::binary);
    string line;
    while(getline(in,line)) correct_lines.push_back(line);

    vector<string> our_lines;
    Buffered_ifstream our_in(filename, ios::binary);
    while(our_in.getline(line)) our_lines.push_back(line);

    ASSERT_EQ(correct_lines.size(), our_lines.size());
    for(int64_t i = 0; i < correct_lines.size(); i++){
        cout << "'" << our_lines[i] << "' '" << correct_lines[i] << "'" << endl;
        ASSERT_EQ(our_lines[i], correct_lines[i]);
    }
}

TEST(TEST_BUFFERED_IO, write_and_read){
    string filename = get_temp_file_manager().create_filename();

    int64_t n = 1e6 * 5; // Code below assumes that this is a multiple of ten
    vector<char> data(n);
    for(int64_t i = 0; i < n; i++) data[i] = rand() % 256;

    {
        Buffered_ofstream<> out(filename, ios::binary);
        
        out.set_buffer_capacity(123);
        
        for(int64_t i = 0; i < n; i += 10)
            out.write(data.data() + i, 10);
    } // End of scope should flush and close the stream

    vector<char> data2(n);
    
    Buffered_ifstream<> in(filename, ios::binary);
    in.set_buffer_capacity(114);

    for(int64_t i = 0; i < n; i += 10)
        in.read(data2.data() + i, 10);

    ASSERT_EQ(data, data2);
} 