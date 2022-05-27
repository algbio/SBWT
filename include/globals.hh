#pragma once

#include <string>
#include <cmath> 
#include "TempFileManager.hh"

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

    Progress_printer(int64_t n_jobs, int64_t total_prints) : n_jobs(n_jobs), processed(0), total_prints(total_prints), next_print(0) {}

    void job_done(){
        if(next_print == processed){
            int64_t progress_percent = round(100 * ((double)processed / n_jobs));
            write_log("Progress: " + to_string(progress_percent) + "%", MINOR);
            next_print += n_jobs / total_prints;
        }
        processed++;
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

