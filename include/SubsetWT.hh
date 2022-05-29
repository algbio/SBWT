#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "sdsl/wavelet_trees.hpp"
#include "globals.hh"
#include <map>
#include <cassert>

namespace sbwt{

template<typename WT_type>
class SubsetWT{

public:

    WT_type ACGT_wt;
    WT_type AC_wt;
    WT_type GT_wt;

    char to_char(bool left, bool right) const{
        if(!left && !right) return '0';
        if(!left && right) return '1';
        if(left && !right) return '2';
        if(left && right) return '3';
        return 0; // Make compiler happy
    }

    std::string bitvector_pair_to_string(sdsl::bit_vector& v1, sdsl::bit_vector& v2) const{
        assert(v1.size() == v2.size());
        string S(v1.size(), '\0');
        for(int64_t i = 0; i < v1.size(); i++){
            S[i] = to_char(v1[i], v2[i]);
        }
        return S;
    }

    SubsetWT(){}

    SubsetWT(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits){
        sdsl::bit_vector AC_bv;
        sdsl::bit_vector GT_bv;
        sdsl::bit_vector A_bv;
        sdsl::bit_vector C_bv;
        sdsl::bit_vector G_bv;
        sdsl::bit_vector T_bv;

        assert(A_bits.size() == C_bits.size() && C_bits.size() == G_bits.size() && G_bits.size() == T_bits.size());
        int64_t n = A_bits.size();
        AC_bv.resize(n);
        GT_bv.resize(n);
        int64_t AC_total = 0;
        int64_t GT_total = 0;
        for(int64_t i = 0; i < n; i++){
            AC_bv[i] = A_bits[i] || C_bits[i];
            GT_bv[i] = G_bits[i] || T_bits[i];
            AC_total += AC_bv[i];
            GT_total += GT_bv[i];
        }

        A_bv.resize(AC_total);
        C_bv.resize(AC_total);
        G_bv.resize(GT_total);
        T_bv.resize(GT_total);

        for(int64_t i = 0, j = 0; i < n; i++){
            if(AC_bv[i]){
                A_bv[j] = A_bits[i];
                C_bv[j] = C_bits[i];
                j++;
            }
        }
        
        for(int64_t i = 0, j = 0; i < n; i++){
            if(GT_bv[i]){
                G_bv[j] = G_bits[i];
                T_bv[j] = T_bits[i];
                j++;
            }
        }

        string ACGT_string = bitvector_pair_to_string(AC_bv, GT_bv);
        string AC_string = bitvector_pair_to_string(A_bv, C_bv);
        string GT_string = bitvector_pair_to_string(G_bv, T_bv);

        sdsl::construct_im(ACGT_wt, ACGT_string.c_str(), 1); // 1: file format is a sequence, not a serialized sdsl object
        sdsl::construct_im(AC_wt, AC_string.c_str(), 1); // 1: file format is a sequence, not a serialized sdsl object
        sdsl::construct_im(GT_wt, GT_string.c_str(), 1); // 1: file format is a sequence, not a serialized sdsl object

    }

    // Count of character c in subsets up to pos, not including pos
    int64_t rank(int64_t pos, char c) const{
        assert(c == 'A' || c == 'C' || c == 'G' || c == 'T');
        if(c == 'A'){
            int64_t x = ACGT_wt.rank(pos, to_char(1,0)) + ACGT_wt.rank(pos, to_char(1,1));
            return AC_wt.rank(x, to_char(1,0)) + AC_wt.rank(x, to_char(1,1));
        }
        if(c == 'C'){
            int64_t x = ACGT_wt.rank(pos, to_char(1,0)) + ACGT_wt.rank(pos, to_char(1,1));
            return AC_wt.rank(x, to_char(0,1)) + AC_wt.rank(x, to_char(1,1));
        }
        if(c == 'G'){
            int64_t x = ACGT_wt.rank(pos, to_char(0,1)) + ACGT_wt.rank(pos, to_char(1,1));
            return GT_wt.rank(x, to_char(1,0)) + GT_wt.rank(x, to_char(1,1));
        }
        if(c == 'T'){
            int64_t x = ACGT_wt.rank(pos, to_char(0,1)) + ACGT_wt.rank(pos, to_char(1,1));
            return GT_wt.rank(x, to_char(0,1)) + GT_wt.rank(x, to_char(1,1));
        }
        return 0;
    }

    int64_t serialize(ostream& os) const{
        int64_t written = 0;
        written += ACGT_wt.serialize(os);
        written += AC_wt.serialize(os);
        written += GT_wt.serialize(os);
        return written;
    }

    void load(istream& is){
        ACGT_wt.load(is);
        AC_wt.load(is);
        GT_wt.load(is);
    }


};

}