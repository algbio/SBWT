#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "globals.hh"
#include <map>

namespace sbwt{

using namespace std;

template <typename bitvector_t, typename rank_support_t>
class SubsetMatrixRank{

    public:

    // Bit vectors
    bitvector_t A_bits;
    bitvector_t C_bits;
    bitvector_t G_bits;
    bitvector_t T_bits;

    // Rank supports
    rank_support_t A_bits_rs;
    rank_support_t C_bits_rs;
    rank_support_t G_bits_rs;
    rank_support_t T_bits_rs;

    // Count of character c in subsets up to pos, not including pos
    int64_t rank(int64_t pos, char c) const{
        if(c == 'A') return A_bits_rs.rank(pos);
        if(c == 'C') return C_bits_rs.rank(pos);
        if(c == 'G') return G_bits_rs.rank(pos);
        if(c == 'T') return T_bits_rs.rank(pos);
        return 0;
    }

    SubsetMatrixRank(){}

    SubsetMatrixRank(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits)
        : A_bits(A_bits), C_bits(C_bits), G_bits(G_bits), T_bits(T_bits){
        sdsl::util::init_support(this->A_bits_rs, &(this->A_bits));
        sdsl::util::init_support(this->C_bits_rs, &(this->C_bits));
        sdsl::util::init_support(this->G_bits_rs, &(this->G_bits));
        sdsl::util::init_support(this->T_bits_rs, &(this->T_bits));
    }

    SubsetMatrixRank(const SubsetMatrixRank& other){
        assert(&other != this); // What on earth are you trying to do?
        operator=(other);
    }

    SubsetMatrixRank& operator=(const SubsetMatrixRank& other){
        if(&other != this){
            this->A_bits = other.A_bits;
            this->C_bits = other.C_bits;
            this->G_bits = other.G_bits;
            this->T_bits = other.T_bits;

            this->A_bits_rs = other.A_bits_rs;
            this->C_bits_rs = other.C_bits_rs;
            this->G_bits_rs = other.G_bits_rs;
            this->T_bits_rs = other.T_bits_rs;

            this->A_bits_rs.set_vector(&(this->A_bits)); 
            this->C_bits_rs.set_vector(&(this->C_bits)); 
            this->G_bits_rs.set_vector(&(this->G_bits)); 
            this->T_bits_rs.set_vector(&(this->T_bits)); 

            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

int64_t serialize(ostream& os) const{
    int64_t written = 0;
    written += A_bits.serialize(os);
    written += C_bits.serialize(os);
    written += G_bits.serialize(os);
    written += T_bits.serialize(os);

    written += A_bits_rs.serialize(os);
    written += C_bits_rs.serialize(os);
    written += G_bits_rs.serialize(os);
    written += T_bits_rs.serialize(os);

    write_log("MatrixRank bit vectors total " + to_string((double)written/A_bits.size()*8) + " bits total per node", LogLevel::MINOR);
    return written;
}

void load(istream& is){
    A_bits.load(is);
    C_bits.load(is);
    G_bits.load(is);
    T_bits.load(is);

    if(std::is_same<sdsl::rank_support_v5<>, rank_support_t>::value){
        // Special handling needed for rank_support_v5 because of a design flaw in sdsl
        A_bits_rs.load(is, &A_bits);
        C_bits_rs.load(is, &C_bits);
        G_bits_rs.load(is, &G_bits);
        T_bits_rs.load(is, &T_bits);
    } else{
        A_bits_rs.load(is);
        C_bits_rs.load(is);
        G_bits_rs.load(is);
        T_bits_rs.load(is);

        A_bits_rs.set_vector(&A_bits);
        C_bits_rs.set_vector(&C_bits);
        G_bits_rs.set_vector(&G_bits);
        T_bits_rs.set_vector(&T_bits);
    }
}

};

}