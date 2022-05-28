#include "suffix_group_optimization.hh"
#include <iostream>

using namespace std;

namespace sbwt{

// Entropy of distribution P
double entropy(const vector<double>& P){
    double ans = 0;
    for(double p : P){
        if(p != 0 && p != 1){
            ans += p * log2(1.0 / p);
        }
    }
    return ans;
}

// Pushes the bits to the left end of the suffix group
void push_bits_left(sdsl::bit_vector& A_bits, 
                    sdsl::bit_vector& C_bits, 
                    sdsl::bit_vector& G_bits, 
                    sdsl::bit_vector& T_bits, 
                    const sdsl::bit_vector& suffix_group_marks){

    for(int64_t i = (int64_t)A_bits.size() - 1; i >= 1; i--){
        if(suffix_group_marks[i] == 0){

            // Push left
            A_bits[i-1] = A_bits[i-1] | A_bits[i];
            C_bits[i-1] = C_bits[i-1] | C_bits[i];
            G_bits[i-1] = G_bits[i-1] | G_bits[i];
            T_bits[i-1] = T_bits[i-1] | T_bits[i];

            // Clear
            A_bits[i] = 0;
            C_bits[i] = 0;
            G_bits[i] = 0;
            T_bits[i] = 0;
        }
    }
}

// Maximally spreads the bits inside a suffix group. Assumes the bits have already been pushed to the left
void spread_bits_after_push_left(sdsl::bit_vector& A_bits,
                                 sdsl::bit_vector& C_bits, 
                                 sdsl::bit_vector& G_bits, 
                                 sdsl::bit_vector& T_bits, 
                                 const sdsl::bit_vector& suffix_group_marks){

    vector<sdsl::bit_vector*> M = {&A_bits, &C_bits, &G_bits, &T_bits}; // Matrix
    for(int64_t i = 0; i < (int64_t)A_bits.size()-1; i++){
        if(suffix_group_marks[i+1] == 0){ // Column i and i+1 have the same suffix group
            // Keep topmost 1-bit where it is, move everything else to the right
            int64_t top = 0;
            while(top < 4 && (*M[top])[i] == 0) top++;
            // top is now the row of the topmost 1-bit in the column
            for(int64_t j = top+1; j < 4; j++){
                (*M[j])[i+1] = (*M[j])[i];
                (*M[j])[i] = 0;
            }
        }
    }
}

sdsl::bit_vector mark_suffix_groups(const sdsl::bit_vector& A_bits,
                                    const sdsl::bit_vector& C_bits, 
                                    const sdsl::bit_vector& G_bits, 
                                    const sdsl::bit_vector& T_bits,
                                    int64_t k){

    int64_t n_nodes = A_bits.size();
    vector<int64_t> C_array(4);

    vector<char> last; // last[i] = incoming character to node i
    last.push_back('$');

    C_array[0] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(A_bits[i]) last.push_back('A');

    C_array[1] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(C_bits[i]) last.push_back('C');

    C_array[2] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(G_bits[i]) last.push_back('G');
    
    C_array[3] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(T_bits[i]) last.push_back('T');

    if(last.size() != n_nodes){
        cerr << "BUG " << last.size() << " " << n_nodes << endl;
        exit(1);
    }

    // Mark suffix group starts
    sdsl::bit_vector suffix_group_starts(n_nodes);
    for(int64_t i = 0; i < n_nodes; i++) suffix_group_starts[i] = 0;

    for(int64_t round = 0; round < k-1; round++){
        for(int64_t i = 0; i < n_nodes; i++){
            if(i == 0 || last[i] != last[i-1])
                suffix_group_starts[i] = 1;
        }

        // Propagate the labels one step forward in the graph
        vector<char> propagated(n_nodes, '$');
        int64_t A_ptr = C_array[0];
        int64_t C_ptr = C_array[1];
        int64_t G_ptr = C_array[2];
        int64_t T_ptr = C_array[3];
        for(int64_t i = 0; i < n_nodes; i++){
            if(A_bits[i]) propagated[A_ptr++] = last[i];
            if(C_bits[i]) propagated[C_ptr++] = last[i];
            if(G_bits[i]) propagated[G_ptr++] = last[i];
            if(T_bits[i]) propagated[T_ptr++] = last[i];
        }
        last = propagated;
    }

    return suffix_group_starts;
}

double compute_column_entropy(const sdsl::bit_vector& A_bits,
                              const sdsl::bit_vector& C_bits, 
                              const sdsl::bit_vector& G_bits, 
                              const sdsl::bit_vector& T_bits){
    map<vector<bool>, int64_t> counts;
    for(int64_t i = 0; i < A_bits.size(); i++){
        vector<bool> column = {(bool)A_bits[i], (bool)C_bits[i], (bool)G_bits[i], (bool)T_bits[i]};
        counts[column]++;
    }
    vector<double> P;
    for(const auto &[key, value] : counts){
        P.push_back((double)value / A_bits.size());
    }
    return entropy(P);
}

} // Namespace sbwt