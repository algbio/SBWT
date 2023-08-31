#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "commands.hh"
#include "globals.hh"

using namespace std;

static vector<string> commands = {"build", "build-variant", "search"};

void print_help(int argc, char** argv){
    (void) argc; // Unused parameter
    cerr << "Available commands: " << endl;
    for(string S : commands) cerr << "   " << argv[0] << " " << S << endl;
    cerr << "Running a command without arguments prints the usage instructions for the command." << endl;
}

int main(int argc, char** argv){

    #ifndef __BMI2__
    cerr << "WARNING: This program was compiled for a CPU without support for the BMI2 instruction set. The performance of the Elias-Fano variants will be very bad." << endl;
    #endif

    sbwt::write_log("Maximum k-mer length is set to " + to_string(MAX_KMER_LENGTH), sbwt::LogLevel::MAJOR);

    if(argc == 1){
        print_help(argc, argv);
        return 0;
    }

    string command = argv[1];
    if(command == "--help" || command == "-h"){
        print_help(argc, argv);
        return 0;
    }

    // Drop the first element of argv
    for(int64_t i = 1; i < argc; i++) argv[i-1] = argv[i];
    argc--;

    try{
        if(command == "build") return build_main(argc, argv);
        else if(command == "search") return search_main(argc, argv);
        else if(command == "build-variant") return build_from_plain_main(argc, argv);
        else{
            throw std::runtime_error("Invalid command: " + command);
            return 1;
        }
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    } catch(const std::exception& e){
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

}
