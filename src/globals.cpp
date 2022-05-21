#include "globals.hh"
#include <algorithm>
#include <chrono>
#include <string>
#include <mutex>
#include <iostream>
#include <iomanip>

using namespace std::chrono;

Temp_File_Manager& get_temp_file_manager(){
    static Temp_File_Manager temp_file_manager; // Singleton
    return temp_file_manager;
}

long long cur_time_millis(){
	return (std::chrono::duration_cast< milliseconds >(high_resolution_clock::now().time_since_epoch())).count();
}

long long cur_time_micros(){
	return (std::chrono::duration_cast< microseconds >(high_resolution_clock::now().time_since_epoch())).count();
}

long long int program_start_millis = cur_time_millis();
long long int program_start_micros = cur_time_micros();

double seconds_since_program_start(){
	return (cur_time_micros() - program_start_micros) / 1000.0;
}

string getTimeString(){
    std::time_t result = std::time(NULL);
    string time = std::asctime(std::localtime(&result));
    return time.substr(0,time.size() - 1); // Trim the trailing newline
}

static LogLevel loglevel = MAJOR;
void set_log_level(LogLevel level){
    loglevel = level;
}
LogLevel get_log_level(){
    return loglevel;
}

static std::mutex write_log_mutex;
void write_log(string message, LogLevel level){
    if(level <= loglevel){
        std::lock_guard<std::mutex> lock(write_log_mutex);
        std::streamsize default_precision = std::cout.precision();

        std::cerr << 
        std::setprecision(4) << std::fixed <<
        seconds_since_program_start() <<
        std::setprecision(default_precision) << 
        " " << getTimeString() << " " << message << std::endl;
    }
}

void check_true(bool condition, string error_message){
    if(!condition){
        throw std::runtime_error(error_message);
    }
}

