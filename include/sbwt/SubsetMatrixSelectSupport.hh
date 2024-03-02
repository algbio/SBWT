#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include "globals.hh"
#include <map>

namespace sbwt{

using namespace std;

// A subset select support based on bitvector select support on indicator bitvectors for each character.
// This class does not own the bit vectors its pointing to. But it does own the select support data.
template<typename bitvector_t>
class SubsetMatrixSelectSupport{

    public:

    // Select supports
    typename bitvector_t::select_1_t A_bits_ss;
    typename bitvector_t::select_1_t C_bits_ss;
    typename bitvector_t::select_1_t G_bits_ss;
    typename bitvector_t::select_1_t T_bits_ss;

    int64_t select(int64_t pos, char c) const{
        if(c == 'A') return A_bits_ss.select(pos);
        if(c == 'C') return C_bits_ss.select(pos);
        if(c == 'G') return G_bits_ss.select(pos);
        if(c == 'T') return T_bits_ss.select(pos);
        return 0;
    }

    SubsetMatrixSelectSupport(){}

    SubsetMatrixSelectSupport(const sdsl::bit_vector* A_bits, const sdsl::bit_vector* C_bits, const sdsl::bit_vector* G_bits, const sdsl::bit_vector* T_bits){
        sdsl::util::init_support(this->A_bits_ss, A_bits);
        sdsl::util::init_support(this->C_bits_ss, C_bits);
        sdsl::util::init_support(this->G_bits_ss, G_bits);
        sdsl::util::init_support(this->T_bits_ss, T_bits);
    }

    SubsetMatrixSelectSupport(const SubsetMatrixSelectSupport& other){
        assert(&other != this); // What on earth are you trying to do?
        operator=(other);
    }

    SubsetMatrixSelectSupport& operator=(const SubsetMatrixSelectSupport& other){
        if(&other != this){

            this->A_bits_ss = other.A_bits_ss;
            this->C_bits_ss = other.C_bits_ss;
            this->G_bits_ss = other.G_bits_ss;
            this->T_bits_ss = other.T_bits_ss;

            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

int64_t serialize(ostream& os) const{
    int64_t written = 0;

    written += A_bits_ss.serialize(os);
    written += C_bits_ss.serialize(os);
    written += G_bits_ss.serialize(os);
    written += T_bits_ss.serialize(os);

    return written;
}

void load(istream& is, const sdsl::bit_vector* A_bits, const sdsl::bit_vector* C_bits, const sdsl::bit_vector* G_bits, const sdsl::bit_vector* T_bits){

    A_bits_ss.load(is);
    C_bits_ss.load(is);
    G_bits_ss.load(is);
    T_bits_ss.load(is);

    A_bits_ss.set_vector(A_bits);
    C_bits_ss.set_vector(C_bits);
    G_bits_ss.set_vector(G_bits);
    T_bits_ss.set_vector(T_bits);
}

};

}