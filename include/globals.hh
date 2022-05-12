#pragma once

#include <string>
#include <cmath> 
using namespace std;

enum LogLevel {OFF = 0, MAJOR = 1, MINOR = 2, DEBUG = 3};

long long cur_time_millis();
long long cur_time_micros();
double seconds_since_program_start();
string getTimeString();
void set_log_level(LogLevel level);
LogLevel get_log_level();
void write_log(string message, LogLevel level);
void check_true(bool condition, string error_message);

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

