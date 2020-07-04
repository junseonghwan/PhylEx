//
//  compact_tssb_state.hpp
//  lib_simul
//
//  Created by Seong-Hwan Jun on 2019-09-11.
//

#ifndef compact_tssb_state_hpp
#define compact_tssb_state_hpp

#include <stdio.h>

#include "clone_node.hpp"
#include "data_util.hpp"
#include "tssb_state.hpp"

class CompactTSSBState
{
    //gsl_matrix *ancestral_matrix;
    //vector<unordered_set<S*> > clustering;
    vector<string> datum2node;
    vector<CloneTreeNodeParam *> datum2param;
    unordered_map<string, double> node2param;
    vector<unsigned int> cluster_labels;
    string newick;
public:
    CompactTSSBState(TSSBState &tssb_state);
    const vector<CloneTreeNodeParam *> get_param() const { return datum2param; };
    inline const vector<unsigned int> &get_cluster_labels() const { return cluster_labels; };
    inline const string get_newick() const { return newick; }
    inline const vector<string> &get_datum2node() const { return datum2node; }
    inline const unordered_map<string, double> &get_node2param() const { return node2param; }
    
    ~CompactTSSBState();
};

#endif /* compact_tssb_state_hpp */
