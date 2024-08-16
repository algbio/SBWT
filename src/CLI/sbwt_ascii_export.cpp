#include "globals.hh"
#include "throwing_streams.hh"
#include "cxxopts.hpp"
#include "SBWT.hh"
#include "SubsetMatrixRank.hh"
#include "SeqIO/SeqIO.hh"
#include "variants.hh"
#include "commands.hh"
#include "throwing_streams.hh"

using namespace std;
using namespace sbwt;

int ascii_export_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Export the SBWT subset sequence in an ASCII format. Each set is written as a string of characters. A non-empty set is a string of DNA characters (ACGT), such that the last character of the set is written in lower case. Empty sets are represented as a single '$'. The representations of the sets are concatenated together. For example, the sequence {A,C}, {A,T}, {}, {C}, {A,G,T} is represented as AcAt$cAGt.");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,index-file", "Index input file.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string indexfile = opts["index-file"].as<string>();
    string outfile = opts["out-file"].as<string>();

    vector<string> variants = get_available_variants();

    throwing_ifstream in(indexfile, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(std::find(variants.begin(), variants.end(), variant) == variants.end()){
        cerr << "Error loading index from file: unrecognized variant specified in the file" << endl;
        return 1;
    }

    write_log("Loading and exporting the index variant " + variant, LogLevel::MAJOR);
    seq_io::Buffered_ofstream<> out(outfile);

    if (variant == "plain-matrix"){
        plain_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        sbwt.ascii_export_sets(out);
    }
    if (variant == "rrr-matrix"){
        rrr_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        sbwt.ascii_export_sets(out);
    }
    if (variant == "mef-matrix"){
        cerr << "Error: Index export does not work for mef-matrix because mef does not implement access to the sets" << endl;
        return 1;
    }
    if (variant == "plain-split"){
        plain_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        sbwt.ascii_export_sets(out);
    }
    if (variant == "rrr-split"){
        rrr_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        sbwt.ascii_export_sets(out);
    }
    if (variant == "mef-split"){
        cerr << "Error: Index export does not work for mef-split because mef does not implement access to the sets" << endl;
        return 1;
    }
    if (variant == "plain-concat"){
        plain_concat_sbwt_t sbwt;
        sbwt.load(in.stream);
        sbwt.ascii_export_sets(out);
    }
    if (variant == "mef-concat"){
        cerr << "Error: Index export does not work for mef-concat because mef does not implement access to the sets" << endl;
        return 1;
    }
    if (variant == "plain-subsetwt"){
        plain_sswt_sbwt_t sbwt;
        sbwt.load(in.stream);
        sbwt.ascii_export_sets(out);
    }
    if (variant == "rrr-subsetwt"){
        rrr_sswt_sbwt_t sbwt;
        sbwt.load(in.stream);
        sbwt.ascii_export_sets(out);
    }

    return 0;
}