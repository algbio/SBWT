#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include <set>
#include <unordered_map>
#include <map>
#include <gtest/gtest.h>
#include "stdlib_printing.hh"
#include "globals.hh"
#include "setup_tests.hh"
#include "test_kmer.hh"
#include "test_small.hh"
#include "test_large.hh"
#include "test_misc.hh"
#include "test_EM_sort.hh"
#include "test_CLI.hh"
#include <cassert>

int main(int argc, char **argv) {
    try{
        setup_tests(argc, argv);
        return RUN_ALL_TESTS();
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    } catch(const std::exception& e){
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}