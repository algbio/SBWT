#include "input_reading.hh"

string figure_out_file_format(string filename){
    for(int64_t i = (int64_t)filename.size()-1; i >= 0; i--){
        if(filename[i] == '.'){
            string end = filename.substr(i);
            
            if(end == ".fasta") return "fasta";
            if(end == ".fna") return "fasta";
            if(end == ".ffn") return "fasta";
            if(end == ".faa") return "fasta";
            if(end == ".frn") return "fasta";
            if(end == ".fa") return "fasta";

            if(end == ".fastq") return "fastq";
            if(end == ".fq") return "fastq";

            if(end == ".gz") return "gzip";

            throw(runtime_error("Unknown file format: " + filename));
        }
    }
    throw(runtime_error("Unknown file format: " + filename));
    return "unknown";
}