//
//  compact_tssb_state.hpp
//  lib_simul
//
//  Created by Seong-Hwan Jun on 2019-09-11.
//

#ifndef compact_tssb_state_hpp
#define compact_tssb_state_hpp

#include <stdio.h>

#include "tssb_state.hpp"

template <class S, class C, class P>
class CompactTSSBState
{
    //gsl_matrix *ancestral_matrix;
    //vector<unordered_set<S*> > clustering;
    vector<string> datum2node;
    vector<P *> datum2param;
    unordered_map<string, double> node2param;
    vector<unsigned int> cluster_labels;
    string newick;
public:
    CompactTSSBState(TSSBState<S,C,P> &tssb_state);
    const vector<P *> get_param() const { return datum2param; };
    inline const vector<unsigned int> &get_cluster_labels() const { return cluster_labels; };
    inline const string get_newick() const { return newick; }
    inline const vector<string> &get_datum2node() const { return datum2node; }
    inline const unordered_map<string, double> &get_node2param() const { return node2param; }
    
    ~CompactTSSBState();
};

template <class S, class C, class P>
CompactTSSBState<S,C,P>::CompactTSSBState(TSSBState<S,C,P> &tssb_state)
{
    newick = write_newick(tssb_state.get_root());
    fill_node_to_param(tssb_state.get_root(), node2param);

    // construct ancestral matrix
    //ancestral_matrix = TSSBState<S,C,P>::get_ancestral_matrix(tssb_state);

    // construct datum2param
    Node<S,P> *node;
    for (size_t i = 0; i < tssb_state.get_data().size(); i++) {
        node = tssb_state.get_node(tssb_state.get_data().at(i));
        datum2param.push_back(new P(node->get_node_parameter()));
        datum2node.push_back(node->get_name());
    }

    // get non-empty nodes
//    vector<Node<S,P> *> nodes;
//    TSSBState<S,C,P>::get_all_nodes(true, tssb_state.get_root(), nodes);
//
//    for (size_t i = 0; i < nodes.size(); i++) {
//        unordered_set<S*> data(nodes[i]->get_data()); // copy the elements
//        clustering.push_back(data);
//    }

    Node<S,P>::get_cluster_labels(tssb_state.get_root(), tssb_state.get_data(), cluster_labels);
}

template <class S, class C, class P>
CompactTSSBState<S,C,P>::~CompactTSSBState()
{
//    if (ancestral_matrix != 0)
//        gsl_matrix_free(ancestral_matrix);
    for (size_t i = 0; i < datum2param.size(); i++) {
        // delete the parameter since it is a copy
        // do not delete S*, since it is used globally in the program
        delete datum2param[i];
    }
}

#endif /* compact_tssb_state_hpp */
