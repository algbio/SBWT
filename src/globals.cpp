#include "globals.hh"
#include "throwing_streams.hh"
#include "SeqIO/SeqIO.hh"
#include "SeqIO/buffered_streams.hh"
#include <algorithm>
#include <chrono>
#include <string>
#include <mutex>
#include <iostream>
#include <iomanip>
#include "zstr/zstr.hpp"

using namespace std::chrono;
using namespace sbwt;

string sbwt::get_rc(const string& S){
    string T = S;
    std::reverse(T.begin(), T.end());
    for(char& c : T) c = get_rc(c);
    return T;
}

vector<string> sbwt::readlines(string filename){
    vector<string> lines;
    string line;
    throwing_ifstream in(filename);
    while(getline(in.stream,line)){
        lines.push_back(line);
    }
    return lines;
}


Temp_File_Manager& sbwt::get_temp_file_manager(){
    static Temp_File_Manager temp_file_manager; // Singleton
    return temp_file_manager;
}

void sbwt::check_readable(string filename){
    throwing_ifstream F(filename); // Throws on failure
}

// Also clears the file
void sbwt::check_writable(string filename){
    throwing_ofstream F(filename, std::ofstream::out | std::ofstream::app); // Throws on failure
}

// Returns the number of bytes written
int64_t sbwt::serialize_string(const string& S, ostream& out){
    int64_t size = S.size();
    out.write((char*)&size, sizeof(size));
    out.write(S.data(), size);
    return sizeof(size) + size;
}

string sbwt::load_string(istream& in){
    int64_t size;
    in.read((char*)&size, sizeof(size));
    string S(size, '\0');
    in.read((char*)&S[0], size); // The C++ standard guarantees that std::string is stored contiguously in memory
    return S;
}

long long sbwt::cur_time_millis(){
    return (std::chrono::duration_cast< milliseconds >(high_resolution_clock::now().time_since_epoch())).count();
}

long long sbwt::cur_time_micros(){
    return (std::chrono::duration_cast< microseconds >(high_resolution_clock::now().time_since_epoch())).count();
}

static long long int program_start_millis = cur_time_millis();
static long long int program_start_micros = cur_time_micros();

double sbwt::seconds_since_program_start(){
    return (cur_time_micros() - program_start_micros) / 1000.0;
}

string sbwt::getTimeString(){
    std::time_t result = std::time(NULL);
    string time = std::asctime(std::localtime(&result));
    return time.substr(0,time.size() - 1); // Trim the trailing newline
}

static LogLevel loglevel = MAJOR;
void sbwt::set_log_level(LogLevel level){
    loglevel = level;
}
LogLevel sbwt::get_log_level(){
    return loglevel;
}

static std::mutex write_log_mutex;
void sbwt::write_log(string message, LogLevel level){
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

void sbwt::check_true(bool condition, string error_message){
    if(!condition){
        throw std::runtime_error(error_message);
    }
}
