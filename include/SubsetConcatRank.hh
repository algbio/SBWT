#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "globals.hh"
#include <map>

namespace sbwt{

using namespace std;

template <typename L_bitvec_t, typename L_select0_t, typename concat_WT_t>
class SubsetConcatRank{

    concat_WT_t concat;
    L_bitvec_t L; // Marks the first element from each set in concat with a zero
    L_select0_t L_ss0;

    public:

    // Count of character c in subsets up to pos, not including pos
    int64_t rank(int64_t pos, char c) const{
        return concat.rank(L_ss0.select(pos+1), c);
    }

    SubsetConcatRank(){}

    SubsetConcatRank(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits){
        assert(A_bits.size() == C_bits.size() && C_bits.size() == G_bits.size() && G_bits.size() == T_bits.size());
    
        std::string concat_str;
        vector<bool> L_vec_bool;
        int64_t n = A_bits.size();
        for(int64_t i = 0; i < n; i++){
            if(A_bits[i] == 1) concat_str.push_back('A');
            if(C_bits[i] == 1) concat_str.push_back('C');
            if(G_bits[i] == 1) concat_str.push_back('G');
            if(T_bits[i] == 1) concat_str.push_back('T');

            if(A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 0) 
                concat_str.push_back('$');

            L_vec_bool.push_back(0);
            while(L_vec_bool.size() < concat_str.size()) L_vec_bool.push_back(1);
        }

        L_vec_bool.push_back(0); // End sentinel to avoid a special case in select.
        
        // Construct the final L bitvector from L_vec_bool
        sdsl::bit_vector L_sdsl(L_vec_bool.size());
        for(int64_t i = 0; i < L_vec_bool.size(); i++) L_sdsl[i] = L_vec_bool[i];
        L = L_bitvec_t(L_sdsl);

        // Init supports
        sdsl::util::init_support(this->L_ss0, &(this->L));
        sdsl::construct_im(concat, concat_str.c_str(), 1); // 1: file format is a sequence, not a serialized sdsl object
    }

    int64_t serialize(ostream& os) const{
        int64_t written = 0;
        written += concat.serialize(os);
        written += L.serialize(os);
        written += L_ss0.serialize(os);
        return written;
    }

    void load(istream& is){
        concat.load(is);
        L.load(is);
        L_ss0.load(is);
        L_ss0.set_vector(&L);
    }


    SubsetConcatRank(const SubsetConcatRank& other){
        assert(&other != this); // What on earth are you trying to do?
        operator=(other);
    }

    SubsetConcatRank& operator=(const SubsetConcatRank& other){
        if(&other != this){
            this->concat = other.concat;
            this->L = other.L;
            this->L_ss0 = other.L_ss0;

            this->L_ss0.set_vector(&(this->L));

            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

};

}