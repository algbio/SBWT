#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include "globals.hh"
#include "SubsetMatrixRank.hh"
#include <map>

namespace sbwt{

using namespace std;

/** A subset select support based on bitvector select support on indicator bitvectors for each character.
* This class does not own the bit vectors its pointing to. But it does own the select support data.
*/
template<typename bitvector_t>
class SubsetMatrixSelectSupport{

    public:

    // Select supports
    typename bitvector_t::select_1_type A_bits_ss;
    typename bitvector_t::select_1_type C_bits_ss;
    typename bitvector_t::select_1_type G_bits_ss;
    typename bitvector_t::select_1_type T_bits_ss;

    int64_t select(int64_t pos, char c) const{
        if(c == 'A') return A_bits_ss.select(pos);
        if(c == 'C') return C_bits_ss.select(pos);
        if(c == 'G') return G_bits_ss.select(pos);
        if(c == 'T') return T_bits_ss.select(pos);
        return 0;
    }

    SubsetMatrixSelectSupport(){}

    /** Warning: this select structure points to internal vectors of `mr`, so the select support
    * can be used only as long as those pointers are valid.
    */
    SubsetMatrixSelectSupport(const SubsetMatrixRank& mr){
        sdsl::util::init_support(this->A_bits_ss, mr.A_bits);
        sdsl::util::init_support(this->C_bits_ss, mr.C_bits);
        sdsl::util::init_support(this->G_bits_ss, mr.G_bits);
        sdsl::util::init_support(this->T_bits_ss, mr.T_bits);
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

    /** Warning: this select structure points to internal vectors of `mr`, so the select support
    * can be used only as long as those pointers are valid.
    */
    void load(istream& is, const SubsetMatrixRank& mr){

        A_bits_ss.load(is);
        C_bits_ss.load(is);
        G_bits_ss.load(is);
        T_bits_ss.load(is);

        A_bits_ss.set_vector(mr.A_bits);
        C_bits_ss.set_vector(mr.C_bits);
        G_bits_ss.set_vector(mr.G_bits);
        T_bits_ss.set_vector(mr.T_bits);
    }

};

}