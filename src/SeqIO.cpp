#include "SeqIO.hh"
#include <algorithm>

namespace sbwt{
namespace SeqIO{

using namespace sbwt;
using namespace sbwt::SeqIO;

const vector<string> fasta_suffixes = {".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"};
const vector<string> fastq_suffixes = {".fastq", ".fq"};

void reverse_complement_c_string(char* S, int64_t len){
    std::reverse(S, S + len);
    for(int64_t i = 0; i < len; i++)
        S[i] = sbwt::get_rc(S[i]);
}

FileFormat figure_out_file_format(string filename){
    Format fasta_or_fastq;
    bool gzipped = false;
    string extension;

    string gzip_suffix = "";
    if(filename.size() >= 3 && filename.substr(filename.size()-3) == ".gz"){
        filename = filename.substr(0, filename.size()-3); // Drop .gz
        gzipped = true;
    }

    for(int64_t i = (int64_t)filename.size()-1; i >= 0; i--){
        if(filename[i] == '.'){
            string ending = filename.substr(i);
            
            if(std::find(fasta_suffixes.begin(), fasta_suffixes.end(), ending) != fasta_suffixes.end()){
                fasta_or_fastq = FASTA;
                extension = ending + (gzipped ? ".gz" : "");
                return {fasta_or_fastq, gzipped, extension};
            }
            
            if(std::find(fastq_suffixes.begin(), fastq_suffixes.end(), ending) != fastq_suffixes.end()){
                fasta_or_fastq = FASTQ;
                extension = ending + (gzipped ? ".gz" : "");
                return {fasta_or_fastq, gzipped, extension};
            }

            throw(runtime_error("Unknown file format: " + filename + (gzipped ? ".gz" : "")));
        }
    }
    throw(runtime_error("Unknown file format: " + filename + (gzipped ? ".gz" : "")));
}

int64_t count_sequences(const string& filename){
    int64_t count = 0;
    if(figure_out_file_format(filename).gzipped){
        sbwt::SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> reader(filename);
        while(reader.get_next_read_to_buffer()) count++;
    } else{
        sbwt::SeqIO::Reader<Buffered_ifstream<std::ifstream>> reader(filename);
        while(reader.get_next_read_to_buffer()) count++;
    }
    return count;
}

} // namespace SeqIO
} // namespace sbwt