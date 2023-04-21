#ifndef DEFINED_MOD_EF_VECTOR
#define DEFINED_MOD_EF_VECTOR

//An implementation of Elias-Fano that is faster and smaller than sd_vector in the SDSL.
//This code was written by: Bella Zhukova, Simon J. Puglisi, and Rajeev Raman.
//
//If you use this please cite the paper below and publish the URL where you downloaded the code.
//
//Danyang Ma, Simon J. Puglisi, Rajeev Raman, Bella Zhukova:
//On Elias-Fano for Rank Queries in FM-Indexes. DCC 2021: 223-232.
//
//SJP/BZ, Helsinki, 2020.

#include <sdsl/int_vector.hpp>
#include <sdsl/sfstream.hpp>
#include <sdsl/util.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/hyb_vector.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/bit_vector_il.hpp>
#include <string>
#include <cstdint>

#include <iostream>
#include <stdexcept>
#include <list>

#if defined(__BMI2__)
#include <immintrin.h>
#endif

#include <fstream>

namespace sbwt{

using namespace sdsl;
using namespace std;


template<
        class t_rank_1     = typename bit_vector::rank_1_type>
class rank_support_mod_ef;



template<
        class t_rank_1     = typename bit_vector::rank_1_type>
class mod_ef_vector {
public:
    typedef bit_vector::size_type                   size_type;
    typedef rank_support_mod_ef<t_rank_1>           rank_1_type;

private:
    bit_vector  m_upper;
    bit_vector  m_lower;

    t_rank_1    m_mef_upper_rank_1;          // rank support on m_upper
    t_rank_1    m_mef_lower_rank_1;          // rank support on m_lower

public:
    size_type      m_m;     // m_m — size (# of bits) in the bitvector
    uint8_t     m_wl = 0;   // width of the lower part, in bits

    const t_rank_1&     mef_upper_rank_1    = m_mef_upper_rank_1;
    const t_rank_1&     mef_lower_rank_1    = m_mef_lower_rank_1;

    friend class rank_support_mod_ef<>;


public:
    mod_ef_vector() {
    }

    mod_ef_vector(const mod_ef_vector& mef)
    {
        copy(mef);
    }

    mod_ef_vector(mod_ef_vector&& mef)
    {
        *this = std::move(mef);
    }

    mod_ef_vector(const bit_vector &b) {
        m_m = b.size();
        bit_vector copyOfB;
        optimize_w(b, m_wl);
        size_type m_bucket_size = 1 << m_wl;

        m_upper = bit_vector((m_m / m_bucket_size) + 1, 0);
        size_type count = 0;

        /* loop over all buckets, deciding which ones to copy and which to discard */
        /* we discard all buckets that are all 0 */
        for (size_type i = 0; i < m_m / m_bucket_size; i++) {
            size_type sum = 0;
            for (size_type j = 0; j < m_bucket_size; j++)
                sum += b[i * m_bucket_size + j];
            if (sum != 0) {
                /* this bucket will be copied since it is non-empty */
                count++;
                m_upper[i] = 1;
            }
        }

        m_upper[m_m / m_bucket_size] = 1; /* final bucket is always copied, even if it is all zeros */

        m_lower = bit_vector((count + 1) * m_bucket_size);

        /* go over all buckets except maybe the final one which may be of
           size < m_bucket_size */
        size_type next = 0;
        for (size_type i = 0; i < m_m / m_bucket_size; i++) {
            if (m_upper[i] == 1) { /* bucket is not all zeros, copy it */
                for (size_type j = 0; j < m_bucket_size; j++) {
                    size_type ind_bot = next * m_bucket_size + j;
                    size_type ind_b = i * m_bucket_size + j;
                    m_lower[ind_bot] = b[ind_b];
                }
                next++;
            }
        }

        /* Copy over last bucket */
        for (size_type i = 0; i < (m_m % m_bucket_size); i++)
            m_lower[next * m_bucket_size + i] = b[(m_m / m_bucket_size) * m_bucket_size + i];

        util::init_support(m_mef_upper_rank_1, &m_upper);
        util::init_support(m_mef_lower_rank_1, &m_lower);
    }

    mod_ef_vector(bit_vector b, uint8_t wl) {

        m_m = b.size();
//        size_type m_n = util::cnt_one_bits(b);       // never used

        m_wl = wl;
        size_type m_bucket_size = 1 << m_wl;

//#ifdef DEBUG
//        cout << "wl: " << std::to_string(m_wl) << " bucket size = " << m_bucket_size << endl;
//#endif

        m_upper = bit_vector((m_m / m_bucket_size) + 1, 0);

        size_type count = 0;

        /* loop over all buckets, deciding which ones to copy and which to discard */
        /* we discard all buckets that are all 0 */

        for (size_type i = 0; i < m_m / m_bucket_size; i++) {
            size_type sum = 0;
            for (size_type j = 0; j < m_bucket_size; j++)
                sum += b[i * m_bucket_size + j];
            if (sum != 0) {
                /* this bucket will be copied since it is non-empty */
                count++;
                m_upper[i] = 1;
            }
        }

        m_upper[m_m / m_bucket_size] = 1; /* final bucket is always copied, even if it is
                   all zeros */

        m_lower = bit_vector((count + 1) * m_bucket_size);

//#ifdef DEBUG
//        cout << "m_lower.size() = " << m_lower.size() << endl;
//#endif

        /* go over all buckets except maybe the final one which may be of
           size < m_bucket_size */

        size_type next = 0;

        for (size_type i = 0; i < m_m / m_bucket_size; i++) {
            if (m_upper[i] == 1) { /* bucket is not all zeros, copy it */
                for (size_type j = 0; j < m_bucket_size; j++) {
                    size_type ind_bot = next * m_bucket_size + j;
                    size_type ind_b = i * m_bucket_size + j;
                    m_lower[ind_bot] = b[ind_b];
                }
                next++;
            }
        }

        /* Copy over last bucket */

        for (size_type i = 0; i < (m_m % m_bucket_size); i++)
            m_lower[next * m_bucket_size + i] = b[(m_m / m_bucket_size) * m_bucket_size + i];

        util::init_support(m_mef_upper_rank_1, &m_upper);
        util::init_support(m_mef_lower_rank_1, &m_lower);

//#ifdef DEBUG
//        cout << "length of input " << m_m << endl;
////        cout << "#1s in input " << m_n << endl;
//        cout << "length of upper " << m_upper.size() << endl;
//        cout << "#1s in upper " << util::cnt_one_bits(m_upper) << endl;
//        cout << "length of lower " << m_lower.size() << endl;
//        cout << "#1s in lower " << util::cnt_one_bits(m_lower) << endl;
//#endif
    }

    //! Returns the size of the original bit vector.
    const size_type& size() const {
        return (m_m);
    }

    ~mod_ef_vector() = default;

    mod_ef_vector& operator=(const mod_ef_vector& v)
    {
        if (this != &v) {
            copy(v);
        }
        return *this;
    }

    mod_ef_vector& operator=(mod_ef_vector&& v)
    {
        if (this != &v) {
            m_m      = v.m_m;
            m_wl     = v.m_wl;
            m_upper  = std::move(v.m_upper);
            m_lower  = std::move(v.m_lower);

            m_mef_upper_rank_1 = std::move(v.m_mef_upper_rank_1);
            m_mef_upper_rank_1.set_vector(&m_upper);
            m_mef_lower_rank_1 = std::move(v.m_mef_lower_rank_1);
            m_mef_lower_rank_1.set_vector(&m_lower);
        }
        return *this;
    }

    //! Serializes the data structure into the given ostream
    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += write_member(m_m, out, child, "size");
        written_bytes += write_member(m_wl, out, child, "wl");
        written_bytes += m_upper.serialize(out, child, "upper");
        written_bytes += m_lower.serialize(out, child, "lower");

        // not necessary to be serialized, can be calculated from previous
        written_bytes += m_mef_upper_rank_1.serialize(out, child, "mef_upper_rank_1");
        written_bytes += m_mef_lower_rank_1.serialize(out, child, "mef_lower_rank_1");

        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Loads the data structure from the given istream.
    void load(std::istream& in)
    {
        read_member(m_m, in);
        read_member(m_wl, in);
        m_upper.load(in);
        m_lower.load(in);

        // can be calculated from previous
        m_mef_upper_rank_1.load(in, &m_upper);
        m_mef_upper_rank_1.set_vector(&m_upper);
        m_mef_lower_rank_1.load(in, &m_lower);
        m_mef_lower_rank_1.set_vector(&m_lower);
    }

private:
    void copy(const mod_ef_vector& v)
    {
        m_m      = v.m_m;
        m_wl     = v.m_wl;
        m_upper  = v.m_upper;
        m_lower  = v.m_lower;

        m_mef_upper_rank_1 = v.m_mef_upper_rank_1;
        m_mef_upper_rank_1.set_vector(&m_upper);
        m_mef_lower_rank_1 = v.m_mef_lower_rank_1;
        m_mef_lower_rank_1.set_vector(&m_lower);
    }

    void optimize_w(bit_vector b, uint8_t &wl) {
        size_type m = b.size();
        //size_type n = util::cnt_one_bits(b);

        size_type best_size = m;
        size_type top_size = m;
        size_type bot_size = 0;

        wl = 0;

        while(b.size() >= 64) {
            wl++;
            //cout << countemptyranges(b, (1<<wl)) << endl;
            shrink(b);
            top_size = b.size();
            bot_size = (util::cnt_one_bits(b)) * (1 << wl) ;
            if (top_size + bot_size < best_size)
                best_size = top_size + bot_size;
            else {
                wl--;
                return;
            }
//            cout << "wl: " << (int) wl;
//            cout << " top: " << std::fixed << std::setprecision(1) << (float) top_size/((1<<20)*8) << "MB";
//            cout << " bot: " << std::fixed << std::setprecision(1) << (float) bot_size/((1<<20)*8) << "MB";
//            cout << " tot: " << std::fixed << std::setprecision(1) << (float) best_size/((1<<20)*8) << "MB";
//            cout << endl;
        }
//        fprintf(stderr, "n = %lu, wl = %lu\n", n, (size_type)wl);
//        cout << "Shouldn't be here." << endl;
        return;
    }

    // Software implementation for the _pext_u64 instruction in the bmi2 instruction set
    uint64_t pext_u64_fallback(uint64_t x, uint64_t m) {
        uint64_t r, s, b;    // Result, shift, mask bit. 

        r = 0; 
        s = 0; 
        do {
            b = m & 1; 
            r = r | ((x & b) << s); 
            s = s + b; 
            x = x >> 1; 
            m = m >> 1; 
        } while (m != 0); 
        return r; 
    } 

    uint64_t shrink(unsigned long long x) {
        #if defined(__BMI2__)
        return((uint64_t) _pext_u64( ((x & 0xAAAAAAAAAAAAAAAAULL) >> 1 | x), 0x5555555555555555ULL));
        #else
        return(pext_u64_fallback( ((x & 0xAAAAAAAAAAAAAAAAULL) >> 1 | x), 0x5555555555555555ULL));
        #endif 
    }

    void shrink(bit_vector &b) {
        for(size_type i = 0, j = 0; i < b.size() - 64; i+=64, j+=32) {
            uint64_t nextword = b.get_int(i);
            uint64_t newword = shrink(nextword);
            //    if(__builtin_popcount(newword) < __builtin_popcount(nextword) / 2) {
            //  std::cout.fill('0');
            //  cout << std::setw(16) << std::hex << nextword << endl << std::setw(16) << newword << endl;
            //}

            b.set_int(j, (newword&0xFFFFFFFFULL), 32);
        }

        b.resize(b.size()/2);
    }

};


template<class t_rank_1>
class rank_support_mod_ef {
public:
    typedef bit_vector::size_type size_type;
    typedef mod_ef_vector<t_rank_1> bit_vector_type;

private:
    const bit_vector_type* m_v;
    size_type m_mask = 0;

public:

    explicit rank_support_mod_ef(const bit_vector_type* v=nullptr)
    {
        set_vector(v);
    }

    size_type rank(size_type i) const {
        #ifdef DEBUG
            if (i > initialBitVectorSize) fprintf(stderr, "rank(%lu / %lu)\n", i, m_v->size());
        #endif

        assert(i <= m_v->size());

        size_type nz_block_id = m_v->mef_upper_rank_1(i >> m_v->m_wl);

        size_type lob = (m_v->m_upper[i >> m_v->m_wl] == 0) ? 0 : (i & m_mask);
        //TODO: speed up by replacing branch with multiply?

        return (m_v->mef_lower_rank_1((nz_block_id << m_v->m_wl) + lob));
    }


    size_type operator()(size_type i)const
    {
        return rank(i);
    }

    size_type size()const
    {
        return m_v->size();
    }

    void set_vector(const bit_vector_type* v=nullptr)
    {
        m_v = v;
        if (v != nullptr) { m_mask = (((size_type) 1 << m_v->m_wl) - 1); }
    }

    rank_support_mod_ef& operator=(const rank_support_mod_ef& rs)
    {
        if (this != &rs) {
            set_vector(rs.m_v);
        }
        return *this;
    }

    void swap(rank_support_mod_ef&) { }

    void load(std::istream& in, const bit_vector_type* v=nullptr)
    {
        read_member(m_mask, in);
        set_vector(v);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += write_member(m_mask, out, child, "mask");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

};

}

#endif
