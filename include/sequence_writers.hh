#pragma once

#include "input_reading.hh" // FASTA_MODE and FASTQ mode come from here. Todo: refactor

class Sequence_Writer_Buffered{

    string fasta_header = ">\n";
    string fastq_header = "@\n";
    string newline = "\n";
    string plus = "+";

    public:

    Buffered_ofstream out;
    LL mode;

    // Writes either fasta or fastq based on the file extension
    Sequence_Writer_Buffered(string filename) : out(filename) {
        string format = figure_out_file_format(filename);
        if(format == "fasta") mode = FASTA_MODE;
        else if(format == "fastq") mode = FASTQ_MODE;
        else throw(runtime_error("Unknown file format: " + filename));
    }

    void write_sequence(const char* seq, LL len){
        if(mode == FASTA_MODE){
            // FASTA format
            out.write(fasta_header.c_str(), 2);
            out.write(seq, len);
            out.write(newline.c_str(), 1);
        } else{
            // FASTQ
            out.write(fastq_header.c_str(), 2);
            out.write(seq, len);
            out.write(newline.c_str(), 1);
            out.write(plus.c_str(), 1);
            out.write(newline.c_str(), 1);
            out.write(seq, len); // Use the read again for the quality values
            out.write(newline.c_str(), 1);
        }
    }

    // Flush the stream. The stream is also automatically flushed when the object is destroyed.
    void flush(){
        out.flush();
    }
};