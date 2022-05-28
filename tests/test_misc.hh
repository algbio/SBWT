#include "setup_tests.hh"
#include "kmc_construct.hh"
#include "globals.hh"
#include <gtest/gtest.h>

using namespace sbwt;

TEST(MISC, test_rc){
    ASSERT_EQ(get_rc('A'), 'T');
    ASSERT_EQ(get_rc('C'), 'G');
    ASSERT_EQ(get_rc('G'), 'C');
    ASSERT_EQ(get_rc('T'), 'A');
    ASSERT_EQ(get_rc('a'), 't');
    ASSERT_EQ(get_rc('c'), 'g');
    ASSERT_EQ(get_rc('g'), 'c');
    ASSERT_EQ(get_rc('t'), 'a');
    ASSERT_EQ(get_rc('N'), 'N');
}

void create_rc_file_test(const string& file_extension){

    vector<string> seqs1 = {"ACAGT", "CGAG", "CGGACG"};
    vector<string> seqs2 = {"AGAT", "GAGA", "AAAAAA"};
    string f1 = get_temp_file_manager().create_filename("",file_extension);
    string f2 = get_temp_file_manager().create_filename("",file_extension);
    vector<string> oldfiles = {f1,f2};

    // Write to files
    SeqIO::Writer w1(f1);
    SeqIO::Writer w2(f2);
    for(string seq : seqs1) w1.write_sequence(seq.c_str(), seq.size());
    for(string seq : seqs2) w2.write_sequence(seq.c_str(), seq.size());
    w1.flush();
    w2.flush();

    // Create file list
    string filelist = get_temp_file_manager().create_filename("",".txt");
    throwing_ofstream filelist_out(filelist);
    filelist_out.stream << f1 << "\n" << f1 << "\n";
    filelist_out.close();

    // Reverse complement
    vector<string> newfiles = SeqIO::create_reverse_complement_files<
                SeqIO::Reader<Buffered_ifstream<std::ifstream>>,
                SeqIO::Writer<Buffered_ofstream<std::ofstream>>>(oldfiles);
    ASSERT_EQ(newfiles.size(), oldfiles.size());

    // Check
    for(LL i = 0; i < newfiles.size(); i++){
        SeqIO::Unbuffered_Reader sr1(oldfiles[i]);
        SeqIO::Unbuffered_Reader sr2(newfiles[i]);
        while(!sr1.done()){
            string s1 = sr1.get_next_query_stream().get_all();
            ASSERT_FALSE(sr2.done());
            string s2 = sr2.get_next_query_stream().get_all();
            logger << s1 << endl << get_rc(s2) << endl << "--" << endl;
            ASSERT_EQ(s1, get_rc(s2));
        }
        ASSERT_TRUE(sr2.done());
    }
 
}

TEST(MISC, create_rc_files){
    create_rc_file_test(".fna");
    create_rc_file_test(".fq");
}