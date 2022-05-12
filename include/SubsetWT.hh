#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "libwheeler/BOSS.hh"
#include "globals.hh"
#include <map>
#include <cassert>

// Handles the cases p = 0 and p = 1 also.
double p_log_p(double p){
    if(p == 0 || p == 1) return 0;
    return p * log2(p);
}

template<typename bitvector_t>
double get_entropy(const bitvector_t& v){
    double ones = 0;
    for(int64_t i = 0; i < v.size(); i++) ones += v[i];

    double n = v.size();
    return n * (- p_log_p(ones/n) - p_log_p((n-ones)/n));
}

template<typename bitvector_t>
double get_pair_entropy(const bitvector_t& v1, const bitvector_t& v2){
    assert(v1.size() == v2.size());
    double c00 = 0;
    double c01 = 0;
    double c10 = 0;
    double c11 = 0;

    for(int64_t i = 0; i < v1.size(); i++){
        if(v1[i] == 0 && v2[i] == 0) c00++;
        if(v1[i] == 1 && v2[i] == 0) c10++;
        if(v1[i] == 0 && v2[i] == 1) c01++;
        if(v1[i] == 1 && v2[i] == 1) c11++;
    }

    double n = v1.size();

    double p00 = c00/n;
    double p01 = c01/n;
    double p10 = c10/n;
    double p11 = c11/n;
    
    return n * (- p_log_p(p00) - p_log_p(p01) - p_log_p(p10) - p_log_p(p11));
}

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

        cout << "Subset WT bitvector entropy per k-mer: " << (get_entropy(AC_bv) + get_entropy(GT_bv) + get_entropy(A_bv) + get_entropy(C_bv) + get_entropy(G_bv) + get_entropy(T_bv)) / n << endl;
        cout << "Subset WT bitvector pair-entropy per k-mer: " << (get_pair_entropy(AC_bv, GT_bv) + get_pair_entropy(A_bv, C_bv) + get_pair_entropy(G_bv, T_bv)) / n << endl;

        cout << "Constructing subset WT rank supports" << endl;

        string ACGT_string = bitvector_pair_to_string(AC_bv, GT_bv);
        string AC_string = bitvector_pair_to_string(A_bv, C_bv);
        string GT_string = bitvector_pair_to_string(G_bv, T_bv);

        sdsl::construct_im(ACGT_wt, ACGT_string.c_str(), 1); // 1: file format is a sequence, not a serialized sdsl object
        sdsl::construct_im(AC_wt, AC_string.c_str(), 1); // 1: file format is a sequence, not a serialized sdsl object
        sdsl::construct_im(GT_wt, GT_string.c_str(), 1); // 1: file format is a sequence, not a serialized sdsl object

        cout << "...Done" << endl;

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