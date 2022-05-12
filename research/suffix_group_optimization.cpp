#include "libwheeler/BOSS.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "globals.hh"
#include "NodeBOSS.hh"
#include "SubsetWT.hh"
#include "SubsetSplitRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetConcatRank.hh"
#include "suffix_group_optimization.hh"
#include "combinatorics.hh"
#include <filesystem>
#include <vector>

//#include "MEF.hpp"

using namespace sdsl;

using namespace std;
typedef long long LL;

void check_readable(string filename){
    throwing_ifstream F(filename); // Throws on failure
}

// Also clears the file
void check_writable(string filename){
    throwing_ofstream F(filename, std::ofstream::out | std::ofstream::app); // Throws on failure
}

class EquivClass{

    public:

    int width; // How many columns are in this class
    set<char> labels; // Tha labels in this class

    bool operator<(const EquivClass& other) const{
        if(this->width != other.width) return this->width < other.width;
        return this->labels < other.labels;
    }

    string to_string() const{
        stringstream ss;
        ss << "(" << width << ", ";
        for(char c : labels) ss << c;
        ss << ")";
        return ss.str();
    }

};

class EquivClassConfig{
    public:

    vector<set<char> > V; // List of columns. A column is represented as the set of characters in that column.

    EquivClassConfig(const vector<set<char> >& V) : V(V) {}

    string to_string() const{
        stringstream ss;
        ss << "(";
        bool first = true;
        for(set<char> S : V){
            if(!first) ss << ", ";
            for(char c : S)
                ss << c;
            first = false;
        }
        ss << ")";
        return ss.str();
    }

};

class WorkingSolution{

    public:

    vector<EquivClass> classes; // Available equivalence classes
    vector<LL> counts; // Number of occurrences of each class
    vector<vector<EquivClassConfig> > available_configs; // Available configs for each class

    WorkingSolution(vector<EquivClass>& classes, vector<LL>& counts)
        : classes(classes), counts(counts) {

        assert(classes.size() == counts.size());
        available_configs.resize(classes.size());

        const set<char> emptyset;
        LL class_index = 0;
        for (const EquivClass& EC : classes) {
            if(EC.labels.size() == 0){
                // Empty set must be considered as a special case because generate_all_set_partitions returns
                // an empty list for an empty set, whereas we would have wanted a list with an empty set
                vector<set<char> > choices;
                while(choices.size() < EC.width) choices.push_back(emptyset);
                available_configs[class_index].push_back(EquivClassConfig(choices));
            } else{ // Nonempty set
                for (vector<set<char>> partition : generate_all_set_partitions(EC.labels)) {
                    if(partition.size() <= EC.width){
                        while(partition.size() < EC.width) // Fill empty columns
                            partition.push_back(emptyset);
                        available_configs[class_index].push_back(EquivClassConfig(partition));
                    }
                }
            }
            class_index++;
        }

        for(LL i = 0; i < classes.size(); i++){
            cout << classes[i].to_string() << " " << counts[i] << endl;
            for(EquivClassConfig ECG : available_configs[i])
                cout << ECG.to_string() << endl;
            cout << "--" << endl;
        }
    }

    string choices_to_string(const vector<LL>& choices){
        stringstream ss;
        for(LL i = 0; i < classes.size(); i++){
            ss << classes[i].to_string() << " " << counts[i] << " " 
            << available_configs[i][choices[i]].to_string() << endl;
        }
        return ss.str();
    }

    // Notes for myself:
    // Initial column entropy: 2.21751
    // Column entropy after pushing to left: 2.22634
    // Current entropy: 2.22634

    void optimize(){
        vector<LL> choices(classes.size());
        for(LL class_idx = 0; class_idx < classes.size(); class_idx++){
            choices[class_idx] = available_configs[class_idx].size()-1; // This is the one where everything is in the same column, I hope
        }
        while(true){
            double current_entropy = get_current_column_entropy(choices);
            cout << "Starting new round" << endl;
            cout << "Current entropy: " << current_entropy << endl;

            // Try to find a change that improves the entropy
            bool improvement_found = false;
            for(LL class_idx = 0; class_idx < classes.size(); class_idx++){
                cout << class_idx << endl;
                for(LL choice_idx = 0; choice_idx < available_configs[class_idx].size(); choice_idx++){;
                    vector<LL> choices_copy = choices;
                    choices_copy[class_idx] = choice_idx; // Make an edit
                    double new_entropy = get_current_column_entropy(choices_copy);
                    if(new_entropy < current_entropy - 1e-9){
                        // Better than epsilon = 1e-9 improvement required to deal with floating point rounding errors
                        improvement_found = true;
                        current_entropy = new_entropy;
                        choices = choices_copy; // Update current choices
                        cout << "Current entropy improved: " << current_entropy << endl;
                    }
                }
            }
            if(!improvement_found) break; // The algorithm has converged
        }
        cout << choices_to_string(choices) << endl;
        cout << "Final entropy: " << get_current_column_entropy(choices) << endl;
    }

    double get_current_column_entropy(const vector<LL>& choices){
        
        map<set<char>, LL> column_counts;
        LL total_count = 0;
        for(LL class_idx = 0; class_idx < classes.size(); class_idx++){
            LL config_idx = choices[class_idx];
            for(set<char> column : available_configs[class_idx][config_idx].V){
                column_counts[column] += counts[class_idx]; // Add one column for every copy of this class
                total_count += counts[class_idx]; // Add one column for every copy of this class
            }
        }

        double entropy = 0;
        for(const auto &[column, count] : column_counts){
            double p = count / (double) total_count;
            if(p != 0 && p != 1) entropy += - p * log2(p);
        }
        return entropy;

    }
};



int main(int argc, char** argv){

    set_log_level(LogLevel::MINOR);

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "This program prints a NodeBOSS representation of the index. Not production quality code.");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,input-matrixboss", "The input matrixBOSS", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("save-boundaries", "File to save suffix group boundaries (optional)", cxxopts::value<string>()->default_value(""))
        ("load-boundaries", "File to load group boundaries from (optional)", cxxopts::value<string>()->default_value(""))
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string outfile = opts["out-file"].as<string>();
    string infile = opts["input-matrixboss"].as<string>();
    string temp_dir = opts["temp-dir"].as<string>();

    string boundaries_outfile = opts["save-boundaries"].as<string>();
    string boundaries_infile = opts["load-boundaries"].as<string>();

    check_writable(outfile);
    check_readable(infile);
    if(boundaries_infile != "") check_readable(boundaries_infile);
    if(boundaries_outfile != "") check_writable(boundaries_outfile);
    std::filesystem::create_directory(temp_dir);

    write_log("Loading the matrixBOSS", LogLevel::MAJOR);
    NodeBOSS<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>> matrixboss;
    throwing_ifstream in(infile, ios::binary);
    matrixboss.load(in.stream);

    LL n_nodes = matrixboss.n_nodes;

    sdsl::bit_vector suffix_group_starts;
    if(boundaries_infile == ""){
        // Mark suffix group starts
        write_log("Computing suffix group boundaries", LogLevel::MAJOR);
        suffix_group_starts = mark_suffix_groups(matrixboss.subset_rank.A_bits,
                                                                matrixboss.subset_rank.C_bits,
                                                                matrixboss.subset_rank.G_bits,
                                                                matrixboss.subset_rank.T_bits,
                                                                matrixboss.C,
                                                                matrixboss.k);
        if(boundaries_outfile != ""){
            write_log("Saving suffix group boundaries to file", LogLevel::MAJOR);
            throwing_ofstream out(boundaries_outfile, ios::binary);
            suffix_group_starts.serialize(out.stream);
        }
    } else{
        throwing_ifstream in(boundaries_infile, ios::binary);
        write_log("Loading suffix group boundaries from file", LogLevel::MAJOR);
        suffix_group_starts.load(in.stream);
    }

    sdsl::bit_vector A_copy = matrixboss.subset_rank.A_bits;
    sdsl::bit_vector C_copy = matrixboss.subset_rank.C_bits;
    sdsl::bit_vector G_copy = matrixboss.subset_rank.G_bits;
    sdsl::bit_vector T_copy = matrixboss.subset_rank.T_bits;

    cout << "Initial column entropy: "
         << compute_column_entropy(A_copy, C_copy, G_copy, T_copy) << endl;

    write_log("Arranging bits", LogLevel::MAJOR);
    push_bits_left(A_copy, C_copy, G_copy, T_copy, suffix_group_starts);

    cout << "Column entropy after pushing to left: "
         << compute_column_entropy(A_copy, C_copy, G_copy, T_copy) << endl;

    map<EquivClass, LL> counts;

//    cout << A_copy << endl << C_copy << endl << G_copy << endl << T_copy << endl;
//    cout << suffix_group_starts << endl;

    for(LL i = 0; i < n_nodes; i++){
        if(suffix_group_starts[i] == 1){
            set<char> labels;
            if(A_copy[i]) labels.insert('A');
            if(C_copy[i]) labels.insert('C');
            if(G_copy[i]) labels.insert('G');
            if(T_copy[i]) labels.insert('T');
            LL width = 1;
            while(i+width <= n_nodes-1 && suffix_group_starts[i+width] == 0) width++;
            EquivClass EC;
            EC.labels = labels;
            EC.width = width;
            counts[EC]++;
        }
    }


    vector<EquivClass> classes_vec;
    vector<LL> counts_vec;

    for(const auto &[EC, count] : counts){
        classes_vec.push_back(EC);
        counts_vec.push_back(count);
    }

    WorkingSolution WS(classes_vec, counts_vec);
    WS.optimize();

}