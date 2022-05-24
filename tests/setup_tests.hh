#pragma once

#include <iostream>
#include <gtest/gtest.h>
#include "../globals.hh"
#include "version.h"
#include "throwing_streams.hh"

class TestLogger{
    public:
    bool verbose = false;

    // This is to make std::endl compile and work with the logger
    TestLogger& operator<<(std::ostream& (*f)(std::ostream&)){
        if(verbose) f(std::cerr);
        return *this;
    }
};

template <typename T>
TestLogger& operator<<(TestLogger& L, const T& t){
    if(L.verbose) cerr << t;
    return L;
}

TestLogger logger; // Pipe things you want to print into this object with the '<<' operator

set<string> get_all_kmers(const vector<string>& input, int64_t k){
    set<string> kmers;
    for(string x : input)
        for(int64_t i = 0; i < x.size() - k + 1; i++)
            kmers.insert(x.substr(i,k));
    return kmers;
}

const std::string generate_random_kmer(int64_t k) {
    std::string s;
    for (int64_t i = 0; i < k; i++) {
        const int r = std::rand() % 4;
        switch (r) {
            case (0): s += 'A'; break;
            case (1): s += 'C'; break;
            case (2): s += 'G'; break;
            case (3): s += 'T'; break;
            default: break;
        }
    }
    return s;
}

void write_seqs_to_fasta_file(const vector<string>& v, const string& filename){
    throwing_ofstream out(filename);
    for(string S : v) out.stream << ">\n" << S << "\n";
}


void enable_test_logging(){logger.verbose = true; }
void disable_test_logging(){logger.verbose = false; }

void setup_tests(int argc, char** argv){
    
    if(system("mkdir -p temp") != 0){
        cerr << "Error creating directory ./temp" << endl;
        exit(1);
    }

    if(system("mkdir -p test_data") != 0){
        cerr << "Error creating directory ./test_data" << endl;
        exit(1);
    }

    if(system("mkdir -p test_out") != 0){
        cerr << "Error creating directory ./test_out" << endl;
        exit(1);
    }
    
    bool verbose = false;
    for(int64_t i = 1; i < argc; i++)
        if(argv[i] == string("--verbose") || argv[i] == string("-v")) verbose = true;

    get_temp_file_manager().set_dir("temp");

    verbose ? enable_test_logging() : disable_test_logging(); // test logger
    verbose ? set_log_level(LogLevel::DEBUG) : set_log_level(LogLevel::OFF); // main logger

    ::testing::InitGoogleTest(&argc, argv);

    srand(247829347);

}