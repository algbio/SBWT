#include "globals.hh"
#include "throwing_streams.hh"
#include "buffered_streams.hh"
#include "SeqIO.hh"
#include <algorithm>
#include <chrono>
#include <string>
#include <mutex>
#include <iostream>
#include <iomanip>
#include "zstr/zstr.hpp"

using namespace std::chrono;
using namespace sbwt;

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


char sbwt::get_rc(char c){
    return rc_table[(unsigned char)c];
}

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

int64_t DNA_to_char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

char char_idx_to_DNA(int64_t i){
    assert(i >= 0 && i < 4);   
    static string ACGT = "ACGT";
    return ACGT[i];
}