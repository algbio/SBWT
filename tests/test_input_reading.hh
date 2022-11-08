#pragma once

#include <gtest/gtest.h>
#include "globals.hh"
#include "setup_tests.hh"
#include "SeqIO.hh"

using namespace sbwt;

void check_sequence_reader_output(const vector<string>& seqs, int64_t mode, string fastafile){
    SeqIO::Reader sr(fastafile, mode);
    int64_t n_seqs_read = 0;
    for(string seq : seqs){
        ASSERT_EQ(sr.get_next_read(), seq);
        n_seqs_read++;
    }
    ASSERT_EQ(n_seqs_read, seqs.size());
    ASSERT_EQ(sr.get_next_read().size(), 0);
}

void check_buffered_sequence_reader_output(const vector<string>& seqs, int64_t mode, string filename){
    SeqIO::Reader sr(filename, mode);
    for(string seq : seqs){
        string next = sr.get_next_read();
        ASSERT_EQ(next, seq);
    }
    ASSERT_TRUE(sr.get_next_read() == ""); // Done
}

TEST(INPUT_PARSING, fasta_basic){
    vector<string> seqs = {"AAGTGCTGTANAYA","ACGTURYKMSWBDHVN-"};
    string fasta = ">\n" + seqs[0] + "\n>\n" + seqs[1] + "\n";
    logger << fasta << endl << seqs << endl;
    string filename = string_to_temp_file(fasta);
    check_sequence_reader_output(seqs, SeqIO::FASTA, filename);
    check_buffered_sequence_reader_output(seqs, SeqIO::FASTA, filename);
}

TEST(INPUT_PARSING, fasta_reverse_complements){
    vector<string> seqs = {"AAGTGCTGTANAYA","ACGTURYKMSWBDHVN-"};
    string fasta = ">\n" + seqs[0] + "\n>\n" + seqs[1] + "\n";
    logger << fasta << endl << seqs << endl;
    string filename = string_to_temp_file(fasta, ".fna");
    SeqIO::Reader reader(filename);
    reader.enable_reverse_complements();
    vector<string> seqs_read;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        ASSERT_EQ(len, strlen(reader.read_buf));
        seqs_read.push_back(string(reader.read_buf));
    }

    ASSERT_EQ(seqs_read[0], "AAGTGCTGTANAYA");
    ASSERT_EQ(seqs_read[1], "TYTNTACAGCACTT");
    ASSERT_EQ(seqs_read[2], "ACGTURYKMSWBDHVN-");
    ASSERT_EQ(seqs_read[3], "-NVHDBWSMKYRUACGT");

}

TEST(INPUT_PARSING, fasta_multiple_lines){
    vector<string> seqs = {"AAGTGCTGTANAYA","ACGTURYKMSWBDHVN-"};
    string fasta;

    // Write 3 chars per line
    for(string seq : seqs){
        fasta += ">\n";
        for(int64_t i = 0; i < (int64_t)seq.size(); i += 3){
            fasta += seq.substr(i,3) + "\n";
        }
    }
    logger << fasta << endl << seqs << endl;
    string filename = string_to_temp_file(fasta);
    check_sequence_reader_output(seqs, SeqIO::FASTA, filename);
    check_buffered_sequence_reader_output(seqs, SeqIO::FASTA, filename);
}

TEST(INPUT_PARSING, fasta_upper_case){
    vector<string> seqs = {"AagTGCtGTaNAYA","AcGTURYKmSWbDHVn-"};
    string fasta;

    for(string seq : seqs) fasta += ">\n" + seq + "\n";
    
    logger << fasta << endl << seqs << endl;
    string filename = string_to_temp_file(fasta);

    for(string& seq : seqs) for(char& c : seq) c = toupper(c); // Upper case for validation
    
    check_sequence_reader_output(seqs, SeqIO::FASTA, filename);
    check_buffered_sequence_reader_output(seqs, SeqIO::FASTA, filename);
}

TEST(INPUT_PARSING, fasta_super_long_line){
    vector<string> seqs;
    seqs.push_back(string(3e6, 'A'));
    seqs.push_back(string(4e5, 'G'));

    string fasta;
    for(string seq : seqs) fasta += ">\n" + seq + "\n";
    string filename = string_to_temp_file(fasta);
    check_sequence_reader_output(seqs, SeqIO::FASTA, filename);
    check_buffered_sequence_reader_output(seqs, SeqIO::FASTA, filename);
}


TEST(INPUT_PARSING, fasta_headers_legacy){
    // Using the legacy SeqIO::Unbuffered_Reader because SeqIO::Reader does not parse SeqIO::Unbuffered_Read_stream
    vector<string> seqs;
    seqs.push_back(string(3e6, 'A'));
    seqs.push_back(string(4e5, 'G'));

    vector<string> headers;
    headers.push_back(string(1e5, 'h'));
    headers.push_back(string(1e6, 'H'));

    string fasta;
    for(int64_t i = 0; i < seqs.size(); i++) fasta += ">" + headers[i] + "\n" + seqs[i] + "\n";
    string filename = string_to_temp_file(fasta);
    SeqIO::Unbuffered_Reader sr(filename, SeqIO::FASTA);
    SeqIO::Unbuffered_Read_stream rs = sr.get_next_query_stream();
    ASSERT_EQ(rs.header, headers[0]);
    rs.get_all();
    rs = sr.get_next_query_stream();
    ASSERT_EQ(rs.header, headers[1]);
    rs.get_all();
}

TEST(INPUT_PARSING, fasta_headers){
    // Using the legacy SeqIO::Unbuffered_Reader because SeqIO::Reader does not parse SeqIO::Unbuffered_Read_stream
    vector<string> seqs;
    seqs.push_back(string(3e6, 'A'));
    seqs.push_back(string(4e5, 'G'));
    seqs.push_back(string(512, 'T')); // Power of two special case?

    vector<string> headers;
    headers.push_back(string(1e5, 'h'));
    headers.push_back(string(1e6, 'H'));
    headers.push_back(string(512, 't')); // Power of two special case?

    string fasta;
    for(int64_t i = 0; i < seqs.size(); i++) fasta += ">" + headers[i] + "\n" + seqs[i] + "\n";
    string filename = string_to_temp_file(fasta, ".fna");
    SeqIO::Reader sr(filename);

    string header;

    sr.get_next_read_to_buffer();
    header = sr.header_buf;
    ASSERT_EQ(headers[0], header);
    
    sr.get_next_read_to_buffer();
    header = sr.header_buf;
    ASSERT_EQ(headers[1], header);

    sr.get_next_read_to_buffer();
    header = sr.header_buf;
    ASSERT_EQ(headers[2], header);
}


TEST(INPUT_PARSING, fastq_basic){
    vector<string> seqs =  {"AAGTGCTGTANAYAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA","ACGTURYKMSWBDHVN-"};
    vector<string> quals = {"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", "IIIIIIIIIIIIIIIII"};
    string fastq = "@\n" + seqs[0] + "\n+\n" + quals[0] + "\n" +
                   "@\n" + seqs[1] + "\n+\n" + quals[1] + "\n";
    logger << fastq << endl << seqs << " " << quals << endl;
    string filename = string_to_temp_file(fastq);
    check_buffered_sequence_reader_output(seqs, SeqIO::FASTQ, filename);
}

TEST(INPUT_PARSING, fastq_reverse_complements){
    vector<string> seqs =  {"AAGTGCTGTANAYA","ACGTURYKMSWBDHVN-"};
    vector<string> quals = {"IIIIIIIIIIIIII","IIIIIIIIIIIIIIIII"};
    string fastq = "@\n" + seqs[0] + "\n+\n" + quals[0] + "\n"
                 + "@\n" + seqs[1] + "\n+\n" + quals[1] + "\n";
    logger << fastq << endl << seqs << endl;
    string filename = string_to_temp_file(fastq, ".fq");
    SeqIO::Reader reader(filename);
    reader.enable_reverse_complements();
    vector<string> seqs_read;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        ASSERT_EQ(len, strlen(reader.read_buf));
        seqs_read.push_back(string(reader.read_buf));
    }

    ASSERT_EQ(seqs_read[0], "AAGTGCTGTANAYA");
    ASSERT_EQ(seqs_read[1], "TYTNTACAGCACTT");
    ASSERT_EQ(seqs_read[2], "ACGTURYKMSWBDHVN-");
    ASSERT_EQ(seqs_read[3], "-NVHDBWSMKYRUACGT");

}


TEST(INPUT_PARSING, fastq_upper_case){
    vector<string> seqs =  {"AAGtGcygTANAynAAaAAAAAAAAAAAAAAAAAAAAAAAAAAaaaAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA","ACGTURYKMSWBDHVN-"};
    vector<string> quals = {"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", "IIIIIIIIIIIIIIIII"};
    string fastq = "@\n" + seqs[0] + "\n+\n" + quals[0] + "\n" +
                   "@\n" + seqs[1] + "\n+\n" + quals[1] + "\n";
    logger << fastq << endl << seqs << " " << quals << endl;
    string filename = string_to_temp_file(fastq);
    for(string& seq : seqs) for(char& c : seq) c = toupper(c); // Upper case for validation
    check_buffered_sequence_reader_output(seqs, SeqIO::FASTQ, filename);
}


TEST(INPUT_PARSING, fastq_super_long_line){
    vector<string> seqs;
    seqs.push_back(string(1e6, 'A'));
    seqs.push_back(string(1e5, 'G'));
    vector<string> quals;
    quals.push_back(string(1e6, 'I'));
    quals.push_back(string(1e5, 'I'));

    string fastq = "@\n" + seqs[0] + "\n+\n" + quals[0] + "\n" +
                   "@\n" + seqs[1] + "\n+\n" + quals[1] + "\n";
    string filename = string_to_temp_file(fastq);
    check_buffered_sequence_reader_output(seqs, SeqIO::FASTQ, filename);
}


// Headers are not stored anymore
TEST(INPUT_PARSING, fastq_headers_legacy){
    // Using the legacy SeqIO::Unbuffered_Reader because SeqIO::Reader does not parse SeqIO::Unbuffered_Read_stream
    vector<string> seqs;
    seqs.push_back(string(1e6, 'A'));
    seqs.push_back(string(1e5, 'G'));

    vector<string> quals;
    quals.push_back(string(1e6, 'I'));
    quals.push_back(string(1e5, 'I'));

    vector<string> headers;
    headers.push_back(string(1e5, 'h'));
    headers.push_back(string(1e6, 'H'));

    string fastq = "@" + headers[0] + "\n" + seqs[0] + "\n+\n" + quals[0] + "\n" +
                   "@" + headers[1] + "\n" + seqs[1] + "\n+\n" + quals[1] + "\n";
    string filename = string_to_temp_file(fastq);
    SeqIO::Unbuffered_Reader sr(filename, SeqIO::FASTQ);
    SeqIO::Unbuffered_Read_stream rs = sr.get_next_query_stream();
    ASSERT_EQ(rs.header, headers[0]);
    rs.get_all();
    rs = sr.get_next_query_stream();
    ASSERT_EQ(rs.header, headers[1]);
    rs.get_all();
}

TEST(INPUT_PARSING, fastq_headers){
    vector<string> seqs;
    seqs.push_back(string(3e6, 'A'));
    seqs.push_back(string(4e5, 'G'));
    seqs.push_back(string(512, 'T')); // Power of two special case?

    vector<string> quals;
    quals.push_back(string(1e6, 'I'));
    quals.push_back(string(1e5, 'I'));
    quals.push_back(string(512, 'I')); // Power of two special case?

    vector<string> headers;
    headers.push_back(string(1e5, 'h'));
    headers.push_back(string(1e6, 'H'));
    headers.push_back(string(512, 't')); // Power of two special case?

    string fastq = "@" + headers[0] + "\n" + seqs[0] + "\n+\n" + quals[0] + "\n" +
                   "@" + headers[1] + "\n" + seqs[1] + "\n+\n" + quals[1] + "\n" +
                   "@" + headers[2] + "\n" + seqs[2] + "\n+\n" + quals[2] + "\n";
    string filename = string_to_temp_file(fastq, ".fq");

    SeqIO::Reader sr(filename);
    string header;

    sr.get_next_read_to_buffer();
    header = sr.header_buf;
    cout << header.size() << endl;
    ASSERT_EQ(headers[0], header);
    
    sr.get_next_read_to_buffer();
    header = sr.header_buf;
    ASSERT_EQ(headers[1], header);

    sr.get_next_read_to_buffer();
    header = sr.header_buf;
    ASSERT_EQ(headers[2], header);
}

TEST(INPUT_PARSING, fastq_things_after_plus){
    vector<string> seqs =  {"AAGTGCTGTANAYAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA","ACGTURYKMSWBDHVN-"};
    vector<string> quals = {"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", "IIIIIIIIIIIIIIIII"};
    string fastq = "@\n" + seqs[0] + "\n+SOMETHING\n" + quals[0] + "\n" +
                   "@\n" + seqs[1] + "\n+SOMETHING2\n" + quals[1] + "\n";
    logger << fastq << endl << seqs << " " << quals << endl;
    string filename = string_to_temp_file(fastq);
    check_buffered_sequence_reader_output(seqs, SeqIO::FASTQ, filename);
}

TEST(INPUT_PARSING, fasta_empty_sequence){
    string filename = string_to_temp_file(">\nAAA\n>\n>\nAAA", ".fna");
    SeqIO::Reader sr(filename);
    try{
        while(true) { 
            string read = sr.get_next_read();
            if(read.size() == 0) break;
        }
        ASSERT_TRUE(false); // Should not come here
    } catch(std::runtime_error& e){
        // This is what was supposed to happen
        logger << "Error thrown as expected" << endl;
        return;
    }
}

TEST(INPUT_PARSING, fastq_empty_sequence){
    string filename = string_to_temp_file("@\nAAA\n+\nIII\n@\n\n+\n\n@\nAAA\n+\nIII\n", ".fq");
    SeqIO::Reader sr(filename);
    try{
        while(true) { 
            string read = sr.get_next_read();
            if(read.size() == 0) break;
        }
        ASSERT_TRUE(false); // Should not come here
    } catch(std::runtime_error& e){
        // This is what was supposed to happen
        logger << "Error thrown as expected" << endl;
        return;
    }
}

TEST(INPUT_PARSING, first_char_sanity_check_fastq){
    string filename = string_to_temp_file("ATGCTAGCTGACTGATCGTACA", ".fq");
    try{
        SeqIO::Reader sr(filename);
        ASSERT_TRUE(false); // Should not come here
    } catch(std::runtime_error& e){
        // This is what was supposed to happen
        logger << "Error thrown as expected" << endl;
        return;
    }
}

TEST(INPUT_PARSING, first_char_sanity_check_fasta){
    string filename = string_to_temp_file("ATGCTAGCTGACTGATCGTACA", ".fna");
    try{
        SeqIO::Reader sr(filename);
        ASSERT_TRUE(false); // Should not come here
    } catch(std::runtime_error& e){
        // This is what was supposed to happen
        logger << "Error thrown as expected" << endl;
        return;
    }
}

TEST(INPUT_PARSING, multi_file){
    vector<string> seqs;
    vector<string> headers;
    for(int64_t i = 0; i < 20; i++){
        seqs.push_back(generate_random_kmer(30));
        headers.push_back(to_string(i));
    }

    vector<string> filenames;
    // Put the sequences into files, 3 sequences per file
    for(int64_t i = 0; i < seqs.size(); i += 3){
        string filename = get_temp_file_manager().create_filename("",".fna");
        filenames.push_back(filename);

        // Add up to 3 sequences
        vector<string> file_seqs;
        vector<string> file_headers;
        for(int64_t j = i; j < i + 3; j++){
            if(j >= seqs.size()) break;
            file_seqs.push_back(seqs[j]);
            file_headers.push_back(headers[j]);
        }

        write_seqs_to_fasta_file(file_seqs, file_headers, filename);
    }

    SeqIO::Multi_File_Reader<> reader(filenames);

    // Test that we can get the sequences and headers back
    int64_t seq_idx = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        logger << reader.header_buf << endl << reader.read_buf << endl;
        ASSERT_EQ(string(reader.header_buf), headers[seq_idx]);
        ASSERT_EQ(string(reader.read_buf), seqs[seq_idx]);
        seq_idx++;
    }

    // Test rewind
    reader.rewind_to_start();
    seq_idx = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        logger << reader.header_buf << endl << reader.read_buf << endl;
        ASSERT_EQ(string(reader.header_buf), headers[seq_idx]);
        ASSERT_EQ(string(reader.read_buf), seqs[seq_idx]);
        seq_idx++;
    }

    // Test reverse complements
    reader.enable_reverse_complements();
    reader.rewind_to_start();
    seq_idx = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        logger << reader.header_buf << endl << reader.read_buf << endl;
        
        string true_read = (seq_idx % 2 == 0) ? seqs[seq_idx/2] : get_rc(seqs[seq_idx/2]);
        string true_header = headers[seq_idx/2];
        ASSERT_EQ(string(reader.header_buf), true_header);
        ASSERT_EQ(string(reader.read_buf), true_read);
        seq_idx++;
    }



}