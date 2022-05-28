#include <set>
#include <vector>
#include <stdint.h>
#include <iostream>

namespace sbwt{

using namespace std;

template<typename T> 
vector<set<T>> generate_all_subsets(const vector<T>& vec){
    vector<set<T>> subsets;
    for(int64_t mask = 0; mask < (1LL << vec.size()); mask++){
        set<T> subset;
        for(int64_t i = 0; i < vec.size(); i++){
            if(mask & (1LL << i))
                subset.insert(vec[i]);
        }
        subsets.push_back(subset);
    }
    return subsets;
}

template<typename T> 
vector<T> concat_lists(const vector<T>& A, const vector<T>& B){
    vector<T> AB;
    for(T x : A) AB.push_back(x);
    for(T x : B) AB.push_back(x);
    return AB; 
}

// If S is the empty set, returns an empty list
// Otherwise returns a list of all set partitions of S
template<typename T> 
vector<vector<set<T>>> generate_all_set_partitions(const set<T>& S){
    // Take the first element and some subsets of other elements to go with it. Repeat.
    vector<vector<set<T>>> partitions;
    if(S.size() == 0) return partitions;

    vector<T> vec(S.begin(), S.end());
    vector<T> tail(vec.begin()+1, vec.end());

    for(set<T> tail_subset : generate_all_subsets(tail)){
        set<T> first_set = {vec[0]};
        for(T x : tail_subset) first_set.insert(x);

        set<T> not_used = S;
        for(T x : first_set) not_used.erase(x);

        if(not_used.size() == 0) partitions.push_back({first_set});
        else{
            for(vector<set<T>>& recursed : generate_all_set_partitions(not_used)){
                partitions.push_back(concat_lists({first_set}, recursed));
            }
        }
    }

    return partitions;

}

}