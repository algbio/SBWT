#pragma once

#include "Kmer.hh"
#include "BOSS.hh"
#include <unordered_map>
#include <stdexcept>

template<typename iter_t>
class Kmer_stream_from_iterator_pair{

public:

    iter_t it;
    iter_t end;

    Kmer_stream_from_iterator_pair(iter_t begin, iter_t end) : it(begin), end(end){}

    bool done(){
        return it == end;
    }

    Kmer<MAX_KMER_LENGTH> next(){
        Kmer<MAX_KMER_LENGTH> x = *it;
        it++;
        return x;
    }

};

class Kmer_stream_in_memory{
public:

    const vector<string>& reads;
    int64_t cur_read_id = 0;
    int64_t cur_read_pos = 0;
    int64_t k;

    /**
     * Assumers reads have characters only from alphabet {A,C,G,T}
     */
    Kmer_stream_in_memory(const vector<string>& reads, int64_t k) : reads(reads), k(k){
        
        cur_read_id = reads.size();

        // Find the first read that has length at least k
        for(int64_t i = 0; i < (int64_t)reads.size(); i++){
            if((int64_t)reads[i].size() >= k){
                cur_read_id = i;
                break;
            }
        }
    }

    bool done(){
        return cur_read_id == (int64_t)reads.size();
    }

    Kmer<MAX_KMER_LENGTH> next(){
        Kmer<MAX_KMER_LENGTH> x(reads[cur_read_id].substr(cur_read_pos++, k));

        // Get ready for the next read
        while(cur_read_id < (int64_t)reads.size() && cur_read_pos + k - 1 >= (int64_t)reads[cur_read_id].size()){
            cur_read_id++;
            cur_read_pos = 0;
        }

        return x;
    }
};

template<typename payload_t>
class Kmer_sorter_in_memory{

private:

vector<pair<Kmer<MAX_KMER_LENGTH>, payload_t> > vec;
bool sorted = false;
int64_t stream_idx = 0;

public:

    void add(Kmer<MAX_KMER_LENGTH> x, payload_t E){
        vec.push_back({x,E});
    }

    void sort(){
        std::sort(vec.begin(), vec.end(), [](const pair<Kmer<MAX_KMER_LENGTH>, payload_t>& A, const pair<Kmer<MAX_KMER_LENGTH>, payload_t>& B){
            return A.first < B.first; // Colexicographic comparison
        });
        sorted = true;
    }

    // This does nothing but it's here to have the same interface as the other sorters.
    void set_mem_budget(int64_t bytes){
        (void)bytes;
        return; 
    }

    void reset_stream(){
        stream_idx = 0;
    }

    bool stream_done() const{
        return stream_idx == (int64_t)vec.size();
    }

    pair<Kmer<MAX_KMER_LENGTH>, payload_t> stream_next(){
        assert(sorted);
        return vec[stream_idx++];
    }

    pair<Kmer<MAX_KMER_LENGTH>, payload_t> peek_next(){
        assert(sorted);
        return vec[stream_idx];
    }


};

// A single sorted stream out of two sorted streams
template<typename payload_t, typename sorter_t>
class Kmer_stream_merger{

    sorter_t& A;
    sorter_t& B;

public:

    Kmer_stream_merger(sorter_t& A, sorter_t& B) : A(A), B(B){}

    bool stream_done(){
        return A.stream_done() && B.stream_done();
    }

    pair<Kmer<MAX_KMER_LENGTH>, payload_t> stream_next(){
        if(A.stream_done()) return B.stream_next();
        if(B.stream_done()) return A.stream_next();
        if(A.peek_next().first < B.peek_next().first) return A.stream_next();
        else return B.stream_next();
    }

};

// Edgemer stream should have functions:
// Kmer next(); // Returns the next (k+1)-mer
// bool done(); // Returns whether all (k+1)-mers have been iterated
template <typename boss_t, typename edgemer_stream>
class BOSS_builder{
    
private:

    struct Edgeset{
        //std::bitset<8> edges;
        uint8_t data;

        Edgeset() : data(0) {}

        void set_bit(uint8_t bit, uint8_t value){
            data = (data & ~(((uint8_t)1) << bit)) | (value << bit);
        }

        bool get_bit(uint8_t bit) const{
            return (data >> bit) & ((uint8_t)1);
        }

        bool have_out(char c) const{
            switch(c){
                case 'A': return get_bit(0);
                case 'C': return get_bit(1);
                case 'G': return get_bit(2);
                case 'T': return get_bit(3);
            }
            assert(false); // Invalid character
            return false;
        }

        bool have_in(char c) const{
            switch(c){
                case 'A': return get_bit(4);
                case 'C': return get_bit(5);
                case 'G': return get_bit(6);
                case 'T': return get_bit(7);
            }
            assert(false); // Invalid character
            return false;
        }

        void set_have_out(char c, bool have){
            switch(c){
                case 'A': set_bit(0, have); return;
                case 'C': set_bit(1, have); return;
                case 'G': set_bit(2, have); return;
                case 'T': set_bit(3, have); return;
            }
            assert(false); // Invalid character
        }

        void set_have_in(char c, bool have){
            switch(c){
                case 'A': set_bit(4, have); return;
                case 'C': set_bit(5, have); return;
                case 'G': set_bit(6, have); return;
                case 'T': set_bit(7, have); return;
            }
            assert(false); // Invalid character
        }

        string to_string(){
            stringstream ss;
            ss << "In: {";
            for(char c : string("ACGT")){
                if(have_in(c)) ss << c;
            }
            ss << "} Out: {";
            for(char c : string("ACGT")){
                if(have_out(c)) ss << c;
            }
            ss << "}";
            return ss.str();
        }

        static int64_t size_in_bytes(){
            return sizeof(data);
        }

        void serialize(ostream& out) const{
            out.write((char*)(&data), sizeof(data));
        }

        void load(istream& in){
            in.read((char*)(&data), sizeof(data));
        }

        // out must have at least size_in_bytes() bytes of space
        void serialize(char* out) const{
            memcpy(out, (char*)(&data), sizeof(data));
        }

        void load(const char* in){
            memcpy((char*)(&data), in, sizeof(data));
        }

    };

    // Returns number of prefixes added
    int64_t process_prefixes(Kmer<MAX_KMER_LENGTH> x, Edgeset E, Kmer_sorter_in_memory<Edgeset>& new_sorter){
        int64_t k = x.get_k();
        int64_t count = 0;
        if(!E.have_in('A') && !E.have_in('C') && !E.have_in('G') && !E.have_in('T')){
            // Add prefixes
            char c = '\0';
            while(true){
                Edgeset E_new;
                if(x.get_k() > 0)
                    E_new.set_have_in('A', true); // This is actually "dollar" but we use 'A' as dollar because the edgeset does not have a bit for dollar
                if(x.get_k() < k)
                    E_new.set_have_out(c, true);
                new_sorter.add(x, E_new);
                count++;
                if(x.get_k() > 0){
                    c = x.last();
                    x.dropright();
                } else break;
            }
        }
        return count;
    }

    // Assumes the old sorter is already sorted and the other is just initialized.
    // Adds the new edges to the new sorter. Returns number of dummies added.
    int64_t add_dummies(Kmer_sorter_in_memory<Edgeset>& old_sorter, Kmer_sorter_in_memory<Edgeset>& new_sorter){
        string ACGT = "ACGT";
        Kmer<MAX_KMER_LENGTH> prev_kmer;
        Edgeset cur_edgeset;
        int64_t loopcount = 0;
        int64_t addcount = 0;

        while(!old_sorter.stream_done()){
            Kmer<MAX_KMER_LENGTH> cur_kmer; Edgeset E;
            tie(cur_kmer,E) = old_sorter.stream_next();

            if(loopcount > 0 && cur_kmer != prev_kmer){
                // New node -> process the collected edges of the previous node
                addcount += process_prefixes(prev_kmer, cur_edgeset, new_sorter);
                cur_edgeset = Edgeset();
            }

            // Add the edge to collection
            for(char c : ACGT){
                if(E.have_in(c)) cur_edgeset.set_have_in(c, true);
                if(E.have_out(c)) cur_edgeset.set_have_out(c, true);
            }

            // Special case for processing the last node
            if(old_sorter.stream_done()) 
                addcount += process_prefixes(cur_kmer, cur_edgeset, new_sorter);

            prev_kmer = cur_kmer;
            loopcount++;
        }
        return addcount;
    }

public:

    boss_t build(edgemer_stream& input){
        if(input.done()) return boss_t(); 

        Kmer_sorter_in_memory<Edgeset> sorter1;
        int64_t k = 0;
        int64_t n_records_written = 0;
        int64_t edge_count = 0;
        while(!input.done()){
            Kmer<MAX_KMER_LENGTH> edgemer = input.next();
            k = edgemer.get_k() - 1; // Edgemers are (k+1)-mers

            char first = edgemer.get(0);
            char last = edgemer.get(edgemer.get_k() - 1);

            Kmer<MAX_KMER_LENGTH> prefix = edgemer.copy().dropright();
            Edgeset prefixE; prefixE.set_have_out(last, true);
            sorter1.add(prefix, prefixE); n_records_written++;

            Kmer<MAX_KMER_LENGTH> suffix = edgemer.copy().dropleft();
            Edgeset suffixE; suffixE.set_have_in(first, true);
            sorter1.add(suffix, suffixE); n_records_written++;

            edge_count++;

        }
        write_log("Sorting " + to_string(edge_count) + " (k+1)-mers", LogLevel::MAJOR);
        sorter1.sort();

        // Dummies
        write_log("Adding dummy (k+1)-mers", LogLevel::MAJOR);
        Kmer_sorter_in_memory<Edgeset> sorter2;
        int64_t n_dummies_disk = add_dummies(sorter1, sorter2);
        sorter2.sort();
        sorter1.reset_stream(); // Rewind back to start because add_dummies iterates over this

        // Build BOSS vectors
        Kmer_stream_merger<Edgeset, Kmer_sorter_in_memory<Edgeset>> merger(sorter1, sorter2);
        
        vector<bool> I, O;
        string outlabels;
        string ACGT = "ACGT";
        Kmer<MAX_KMER_LENGTH> prev_kmer;
        Edgeset cur_edgeset;
        int64_t loopcount = 0;

        if(n_dummies_disk == 0){
            // There were no dummies. But the BOSS structure wants there always
            // the be a source node. So we need to add it separately.
            I.push_back(1); // No in-edges
            O.push_back(1); // No out-edges
        }

        auto add_node = [&](){
            I.push_back(1);
            O.push_back(1);
            for(char c : ACGT){
                if(cur_edgeset.have_in(c)) I.push_back(0);
                if(cur_edgeset.have_out(c)){
                    O.push_back(0);
                    outlabels += c;
                }
            }
        };

        write_log("Constructing Wheeler BOSS components.", LogLevel::MAJOR);
        while(!merger.stream_done()){
            Kmer<MAX_KMER_LENGTH> cur_kmer; Edgeset E;
            tie(cur_kmer,E) = merger.stream_next();

            if(loopcount > 0 && cur_kmer != prev_kmer){
                // New node -> process the collected edges of a the previous node
                add_node();
                cur_edgeset = Edgeset();
            }

            // Add the edge to collection
            for(char c : ACGT){
                if(E.have_in(c)) cur_edgeset.set_have_in(c, true);
                if(E.have_out(c)){
                    cur_edgeset.set_have_out(c, true);
                }
            }

            // Special case for processing the last node
            if(merger.stream_done()) add_node();

            prev_kmer = cur_kmer;
            loopcount++;
        }

        assert(I.size() == O.size());

        return boss_t(outlabels, I, O, k);

    }

};
