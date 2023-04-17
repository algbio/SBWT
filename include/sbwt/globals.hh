#pragma once

#include <string>
#include <cmath> 
#include <cassert>
#include "TempFileManager.hh"

namespace sbwt{

using namespace std;

#ifndef MAX_KMER_LENGTH
#define MAX_KMER_LENGTH 32
#endif

// Table mapping ascii values of characters to their reverse complements,
// lower-case to lower case, upper-case to upper-case. Non-ACGT characters
// are mapped to themselves.
static constexpr unsigned char rc_table[256] =
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 66, 71, 68, 69, 70, 67, 72, 73,
74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 65, 85, 86, 87, 88, 89, 90, 91,
92, 93, 94, 95, 96, 116, 98, 103, 100, 101, 102, 99, 104, 105, 106, 107,
108, 109, 110, 111, 112, 113, 114, 115, 97, 117, 118, 119, 120, 121, 122,
123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137,
138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,
153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182,
183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212,
213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227,
228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242,
243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};

// Table mapping ACGT to 0123
static constexpr int8_t from_ACGT_to_0123_lookup_table[256] = 
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

// Table mapping 0123 to ACGT
static constexpr char from_0123_to_ACGT_lookup_table[4] = {'A','C','G','T'};

// ACGT -> 0123
constexpr int64_t DNA_to_char_idx(char c){
    return from_ACGT_to_0123_lookup_table[(uint8_t)c];
}

// 0123 -> ACGT
constexpr char char_idx_to_DNA(int64_t i){
    assert(i >= 0 && i < 4);   
    return from_0123_to_ACGT_lookup_table[i];
}

// ACGT -> TGCA
constexpr char get_rc(char c){
    return rc_table[(unsigned char)c];
}


enum LogLevel {OFF = 0, MAJOR = 1, MINOR = 2, DEBUG = 3};

vector<string> readlines(string filename);
long long cur_time_millis();
long long cur_time_micros();
double seconds_since_program_start();
string getTimeString();
void set_log_level(LogLevel level);
LogLevel get_log_level();
void write_log(string message, LogLevel level);
void check_true(bool condition, string error_message);
string get_rc(const string& s); // Reverse complement
Temp_File_Manager& get_temp_file_manager();

void check_readable(string filename);
void check_writable(string filename);

Temp_File_Manager& get_temp_file_manager();

int64_t serialize_string(const string& S, ostream& out); // Returns the number of bytes written
string load_string(istream& in); // Loads string serialized by serialize_string

class Progress_printer{

    public:

    int64_t n_jobs;
    int64_t processed;
    int64_t total_prints;
    int64_t next_print;
    bool first_print;

    Progress_printer(int64_t n_jobs, int64_t total_prints) : n_jobs(n_jobs), processed(0), total_prints(total_prints), next_print(0), first_print(true) {}

    void job_done(){
        if(sbwt::get_log_level() >= LogLevel::MINOR){
            if(next_print == processed){
                //string erase(current_string.size() + 1, '\b'); // Backspace characters. +1 For the endline
                if(!first_print) cerr << '\r' << flush; // Erase current line
                first_print = false;
                
                int64_t progress_percent = round(100 * ((double)processed / n_jobs));
                cerr << to_string(progress_percent) + "%" << flush;

                next_print += n_jobs / total_prints;
            }
            processed++;
            if(processed == n_jobs) cerr << "\r100%" << endl; // No more prints coming
        }
    }

};


class Argv{ // Class for turning a vector<string> into char**
private:

    // Forbid copying the class because it wont work right
    Argv(Argv const& other);
    Argv& operator=(Argv const& other);

public:

    char** array = NULL;
    int64_t size = 0;

    Argv(vector<string> v){
        array = (char**)malloc(sizeof(char*) * v.size());
        // Copy contents of v into array
        for(int64_t i = 0; i < v.size(); i++){
            char* s = (char*)malloc(sizeof(char) * (v[i].size() + 1)); // +1: space for '\0' at the end
            for(int64_t j = 0; j < v[i].size(); j++){
                s[j] = v[i][j]; // Can't use strcpy because s.c_str() is const
            }
            s[v[i].size()] = '\0';
            array[i] = s;
        }
        size = v.size();
    }

    ~Argv(){
        for(int64_t i = 0; i < size; i++) free(array[i]);
        free(array);
    }

};

}