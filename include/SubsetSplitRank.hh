#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "globals.hh"
#include "SBWT.hh"
#include <map>
#include <sdsl/wavelet_trees.hpp>

namespace sbwt{

template<typename X_bitvector_t, typename X_bitvector_rank_t, typename Z_bitvector_t, typename Z_rank_support_t>
class SubsetSplitRank{

typedef sdsl::wt_blcd<sdsl::bit_vector, // Underlying bit vector type
                      sdsl::rank_support_v5<1>, // Rank support
                      sdsl::select_support_scan<1>, // Scan = no fast select support
                      sdsl::select_support_scan<0>> WT_type; // Scan = no fast select support

public:

    X_bitvector_t X; // Marks which columns have other than 1 outgoing edge (0 or >= 2)
    WT_type Y; // The outgoing labels in columns with exactly one outgoing edge
    Z_bitvector_t Z_A; // The row of 'A' in the Z matrix, in columns with != 1 one-bit
    Z_bitvector_t Z_C; // The row of 'C' in the Z matrix, in columns with != 1 one-bit
    Z_bitvector_t Z_G; // The row of 'G' in the Z matrix, in columns with != 1 one-bit
    Z_bitvector_t Z_T; // The row of 'T' in the Z matrix, in columns with != 1 one-bit

    X_bitvector_rank_t X_rs; 
    Z_rank_support_t Z_A_rs; 
    Z_rank_support_t Z_C_rs;
    Z_rank_support_t Z_G_rs;
    Z_rank_support_t Z_T_rs;

    // Returns the number of bytes written
    int64_t serialize(ostream& os) const{
        int64_t written = 0;
        written += X.serialize(os);
        written += Y.serialize(os);
        written += Z_A.serialize(os);
        written += Z_C.serialize(os);
        written += Z_G.serialize(os);
        written += Z_T.serialize(os);
        
        written += X_rs.serialize(os);
        written += Z_A_rs.serialize(os);
        written += Z_C_rs.serialize(os);
        written += Z_G_rs.serialize(os);
        written += Z_T_rs.serialize(os);
        return written;
    }

    void load(istream& is){
        X.load(is);
        Y.load(is);
        Z_A.load(is);
        Z_C.load(is);
        Z_G.load(is);
        Z_T.load(is);

        if(std::is_same<sdsl::rank_support_v5<>, X_bitvector_rank_t>::value){
            // Special handling needed for rank_support_v5 because of a design flaw in sdsl
            X_rs.load(is, &X);
        } else{
            X_rs.load(is);
            X_rs.set_vector(&X);
        }
        if(std::is_same<sdsl::rank_support_v5<>, Z_rank_support_t>::value){
            // Special handling needed for rank_support_v5 because of a design flaw in sdsl
            Z_A_rs.load(is, &Z_A);
            Z_C_rs.load(is, &Z_C);
            Z_G_rs.load(is, &Z_G);
            Z_T_rs.load(is, &Z_T);
        } else{
            Z_A_rs.load(is);
            Z_C_rs.load(is);
            Z_G_rs.load(is);
            Z_T_rs.load(is);
                        
            Z_A_rs.set_vector(&Z_A);
            Z_C_rs.set_vector(&Z_C);
            Z_G_rs.set_vector(&Z_G);
            Z_T_rs.set_vector(&Z_T);
        }
    }

    SubsetSplitRank() {}

    SubsetSplitRank(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits){
        int64_t n = A_bits.size();
        int64_t n_b = 0; // Number of branching nodes plus the nodes that do not have outedges
        int64_t n_u = 0; // Number of nodes with exactly one outgoing edge
        for(int64_t i = 0; i < n; i++){
            if(A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 1) n_u++;
            else n_b++;
        }

        //cout << "n_b n_u " << n_b << " " << n_u << endl;

        sdsl::bit_vector X_plain(n);
        sdsl::bit_vector Z_A_plain(n_b);
        sdsl::bit_vector Z_C_plain(n_b);
        sdsl::bit_vector Z_G_plain(n_b);
        sdsl::bit_vector Z_T_plain(n_b);
        std::string Y_str(n_u, '\0');

        for(int64_t i = 0, Y_str_idx = 0, Z_idx = 0; i < n; i++){
            if(A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 1){
                // One outgoing label
                X_plain[i] = 0;
                if(A_bits[i] == 1) Y_str[Y_str_idx] = 'A';
                if(C_bits[i] == 1) Y_str[Y_str_idx] = 'C';
                if(G_bits[i] == 1) Y_str[Y_str_idx] = 'G';
                if(T_bits[i] == 1) Y_str[Y_str_idx] = 'T';
                Y_str_idx++;
            } else{
                // 0 or > 1 outgoing labels
                X_plain[i] = 1;
                Z_A_plain[Z_idx] = A_bits[i];
                Z_C_plain[Z_idx] = C_bits[i];
                Z_G_plain[Z_idx] = G_bits[i];
                Z_T_plain[Z_idx] = T_bits[i];
                Z_idx++;
            }
        }

        X = X_bitvector_t(X_plain);
        sdsl::construct_im(Y, Y_str.c_str(), 1); // 1: file format is a sequence, not a serialized sdsl object
        Z_A = Z_bitvector_t(Z_A_plain);
        Z_C = Z_bitvector_t(Z_C_plain);
        Z_G = Z_bitvector_t(Z_G_plain);
        Z_T = Z_bitvector_t(Z_T_plain);

        sdsl::util::init_support(X_rs, &X);
        sdsl::util::init_support(Z_A_rs, &Z_A);
        sdsl::util::init_support(Z_C_rs, &Z_C);
        sdsl::util::init_support(Z_G_rs, &Z_G);
        sdsl::util::init_support(Z_T_rs, &Z_T);

    }

    SubsetSplitRank(const SubsetSplitRank& other){
        assert(&other != this); // What on earth are you trying to do?
        operator=(other);
    }

    SubsetSplitRank& operator=(const SubsetSplitRank& other){
        if(&other != this){
            this->X = other.X;
            this->Y = other.Y;
            this->Z_A = other.Z_A;
            this->Z_C = other.Z_C;
            this->Z_G = other.Z_G;
            this->Z_T = other.Z_T;

            this->X_rs = other.X_rs;
            this->Z_A_rs = other.Z_A_rs;
            this->Z_C_rs = other.Z_C_rs;
            this->Z_G_rs = other.Z_G_rs;
            this->Z_T_rs = other.Z_T_rs;

            this->X_rs.set_vector(&(this->X)); 
            this->Z_A_rs.set_vector(&(this->Z_A)); 
            this->Z_C_rs.set_vector(&(this->Z_C)); 
            this->Z_G_rs.set_vector(&(this->Z_G)); 
            this->Z_T_rs.set_vector(&(this->Z_T)); 

            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

    int64_t rank(int64_t pos, char c) const{
        int64_t rank1 = X_rs.rank(pos);
        int64_t rank0 = pos - rank1;
        int64_t Y_count = Y.rank(rank0, c);
        switch(c){
            case 'A': return Y_count + Z_A_rs.rank(rank1);
            case 'C': return Y_count + Z_C_rs.rank(rank1);
            case 'G': return Y_count + Z_G_rs.rank(rank1);
            case 'T': return Y_count + Z_T_rs.rank(rank1);
            default: cerr << "Error: Rank called with non-ACGT character" << endl; exit(1);
        }
    }

};

}