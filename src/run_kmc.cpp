
#include "run_kmc.hh"
#include "globals.hh"
#include "KMC/kmc_api/kmc_file.h"
#include "KMC/include/kmc_runner.h"
#include "KMC_code.hh"
#include "kmc_construct_helper_classes.hh"
#include "SeqIO/SeqIO.hh"

namespace sbwt{

using namespace std;
using namespace kmc_tools;

// See header for description
pair<string, int64_t> run_kmc(const vector<string>& input_files, int64_t k, int64_t n_threads, int64_t ram_gigas, int64_t min_abundance, int64_t max_abundance){

    write_log("Running KMC counter", LogLevel::MAJOR);

    string KMC_db_file_prefix = get_temp_file_manager().create_filename("kmers");

    KMC::Stage1Params stage1Params;

    string f = input_files[0]; // First input file
    seq_io::FileFormat format = seq_io::figure_out_file_format(f);

    for(string f2 : input_files){
        seq_io::FileFormat format2 = seq_io::figure_out_file_format(f2);
        if(format.format != format2.format || format.gzipped != format2.gzipped){
            throw std::runtime_error("Error: all input files must have the same format");
        }
    }

    stage1Params.SetInputFiles(input_files)
        .SetKmerLen(k)
        .SetNThreads(n_threads)
        .SetMaxRamGB(ram_gigas)
        .SetInputFileType(format.format == seq_io::FASTA ? KMC::InputFileType::MULTILINE_FASTA : KMC::InputFileType::FASTQ)
        .SetCanonicalKmers(false)
        .SetTmpPath(get_temp_file_manager().get_dir());

    KMC::Runner kmc;

    auto stage1Results = kmc.RunStage1(stage1Params);

    uint32_t ramForStage2 = ram_gigas;
    KMC::Stage2Params stage2Params;
    stage2Params.SetNThreads(n_threads)
        .SetMaxRamGB(ramForStage2)
        .SetCutoffMin(min_abundance)
        .SetCutoffMax(max_abundance)
        .SetOutputFileName(KMC_db_file_prefix)
        .SetStrictMemoryMode(true);

    auto stage2Results = kmc.RunStage2(stage2Params);

    int64_t n_kmers = stage2Results.nUniqueKmers - stage2Results.nBelowCutoffMin - stage2Results.nAboveCutoffMax;

    write_log("Sorting KMC database", LogLevel::MAJOR);

    try{
        sort_kmc_db(KMC_db_file_prefix, KMC_db_file_prefix + "-sorted");
    } catch(KMCAlreadySortedException& e){
        // Just copy the database
        std::filesystem::copy(KMC_db_file_prefix + ".kmc_pre", KMC_db_file_prefix + "-sorted.kmc_pre");
        std::filesystem::copy(KMC_db_file_prefix + ".kmc_suf", KMC_db_file_prefix + "-sorted.kmc_suf");
    }

    // Delete the unsorted KMC database files. The temp file manager can not do this because
    // KMC appends suffixes to the filename and the manager does not know about that.
    std::filesystem::remove(KMC_db_file_prefix + ".kmc_pre");
    std::filesystem::remove(KMC_db_file_prefix + ".kmc_suf");

    // Clean up the KMC global singleton config state because it seems that it's left
    // in a partial state sometimes, which messes up our code if we call KMC again later.
    CConfig::GetInstance().input_desc.clear();
    CConfig::GetInstance().headers.clear();
    CConfig::GetInstance().simple_output_desc.clear();
    CConfig::GetInstance().transform_output_desc.clear();

    return {KMC_db_file_prefix + "-sorted", n_kmers};
}

void sort_kmc_db(const string& input_db_file, const string& output_db_file){
    vector<string> args = {"kmc_tools", "transform", input_db_file, "sort", output_db_file};
    Argv argv(args);

    CParametersParser params_parser(argv.size, argv.array);
    params_parser.Parse();
    if (params_parser.validate_input_dbs())
    {
        params_parser.SetThreads();
        CApplication<KMER_WORDS> app(params_parser);
        app.Process();
    }
}

} // Namespace