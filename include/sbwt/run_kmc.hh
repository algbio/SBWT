#pragma once

#include <string>
#include <vector>
#include <utility>

// A wrapper to KMC construction and sorting

namespace sbwt{

using namespace std;

// Returns the KMC database prefix and the number of distinct k-mers that had abundance within the given bounds
pair<string, int64_t> run_kmc(const vector<string>& input_files, int64_t k, int64_t n_threads, int64_t ram_gigas, int64_t min_abundance, int64_t max_abundance);

// Sort a KMC database
void sort_kmc_db(const string& input_db_file, const string& output_db_file);

}