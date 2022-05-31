#pragma once

#include <string>
#include <cmath> 
#include "TempFileManager.hh"

namespace sbwt{

using namespace std;

#ifndef MAX_KMER_LENGTH
#define MAX_KMER_LENGTH 32
#endif

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
char get_rc(char c); // Reverse complement
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
        if(sbwt::get_log_level() > LogLevel::MINOR){
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