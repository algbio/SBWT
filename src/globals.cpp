#include "globals.hh"
#include "throwing_streams.hh"
#include "buffered_streams.hh"
#include "input_reading.hh"
#include "sequence_writers.hh"
#include <algorithm>
#include <chrono>
#include <string>
#include <mutex>
#include <iostream>
#include <iomanip>

using namespace std::chrono;

char get_rc(char c){
    if(toupper(c) == 'A') return 'T';
    else if(toupper(c) == 'C') return 'G';
    else if(toupper(c) == 'G') return 'C';
    else if(toupper(c) == 'T') return 'A';
    else return c; // Don't touch characters like 'N'
}

string get_rc(const string& S){
    string T = S;
    std::reverse(T.begin(), T.end());
    for(char& c : T) c = get_rc(c);
    return T;
}

vector<string> readlines(string filename){
    vector<string> lines;
	string line;
	throwing_ifstream in(filename);
	while(getline(in.stream,line)){
		lines.push_back(line);
	}
    return lines;
}


Temp_File_Manager& get_temp_file_manager(){
    static Temp_File_Manager temp_file_manager; // Singleton
    return temp_file_manager;
}

void check_readable(string filename){
    throwing_ifstream F(filename); // Throws on failure
}

// Also clears the file
void check_writable(string filename){
    throwing_ofstream F(filename, std::ofstream::out | std::ofstream::app); // Throws on failure
}

// Returns the number of bytes written
int64_t serialize_string(const string& S, ostream& out){
    int64_t size = S.size();
    out.write((char*)&size, sizeof(size));
    out.write(S.data(), size);
    return sizeof(size) + size;
}

string load_string(istream& in){
    int64_t size;
    in.read((char*)&size, sizeof(size));
    string S(size, '\0');
    in.read((char*)&S[0], size); // The C++ standard guarantees that std::string is stored contiguously in memory
    return S;
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

vector<string> create_reverse_complement_files(const vector<string>& files){
    vector<string> newfiles;
    for(string f : files){
        Sequence_Reader_Buffered sr(f);
        int64_t mode = sr.get_mode();

        string f_rev = get_temp_file_manager().create_filename("", mode == FASTA_MODE ? ".fna" : ".fastq");
        newfiles.push_back(f_rev);
        Sequence_Writer_Buffered out(f_rev);
        
        while(true) { 
            LL len = sr.get_next_read_to_buffer();
            if(len == 0) break;

            // Reverse complement
            std::reverse(sr.read_buf, sr.read_buf + len);
            for(LL i = 0; i < len; i++) sr.read_buf[i] = get_rc(sr.read_buf[i]);

            out.write_sequence(sr.read_buf, len);
        }
    }
    return newfiles;
}
