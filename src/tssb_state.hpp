//
//  tssb_state.hpp
//  tssb
//
//  Created by Seong-Hwan Jun on 2018-07-17.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef tssb_state_hpp
#define tssb_state_hpp

#include <chrono>
#include <cmath>
#include <map>
#include <queue>
#include <stdio.h>
#include <set>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <eigen3/Eigen/Dense>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include "clone_node.hpp"
#include "node.hpp"
#include "model_params.hpp"
#include "numerical_utils.hpp"

using namespace std;

template <class S, class P> class Node;

template <class S, class C, class P>
class TSSBState
{
    Node<S,P> *root;

    double log_lik = DOUBLE_NEG_INF;
    double log_lik_bulk = DOUBLE_NEG_INF;
    double log_lik_sc = DOUBLE_NEG_INF;

    vector<S *> *bulk_data;
    vector<C *> *sc_data;

    unordered_map<const S *, Node<S,P> *> datum2node;

    // vector index for sc_cache is the same as the sc_data idx
    // for eaach cell, store the mapping from node to double (log likelihood of being assigned to that node)
    vector<unordered_map<Node<S,P> *, double> > sc_cache;

    // function pointer to compute likelihood of assigning bulk datum to node
    double (*log_lik_datum)(const Node<S,P> *v, const S *s, const ModelParams &params) = 0;
    // function pointer to compute likelihood of single cell at site s if it has the variant given by has_snv
    double (*log_lik_sc_at_site)(const S *s, const C *c, bool has_snv, const ModelParams &params) = 0;

    // private constructor -- used by simulator, where data has not been generated yet
    TSSBState(Node<S,P> *root);

    // helper functions for initializing assignment of SNV
    //void initialize_data_assignment(const gsl_rng *random, S *datum, const ModelParams &params);
    void initialize_data_assignment(const gsl_rng *random, size_t mut_id, const ModelParams &params);
    // helper functions for performing slice sampling on an SNV
    //void slice_sample_data_assignment(const gsl_rng *random, S *datum, const ModelParams &params);
    void slice_sample_data_assignment(const gsl_rng *random, size_t mut_id, const ModelParams &model_params);
    // assign data to a node
    //void assign_data_point(Node<S,P> *curr_node, Node<S,P> *new_node, S *datum, const ModelParams &model_params);
    void assign_data_point(Node<S,P> *curr_node, Node<S,P> *new_node, size_t mut_id, const ModelParams &model_params);

    void initialize_sc_cache(const ModelParams &model_params);
    void initialize_sc_cache(size_t c, Node<S,P> *v, unordered_set<const S*> &snvs, const ModelParams &model_params);
    //void update_sc_cache(Node<S,P> *curr_node, Node<S,P> *new_node, S *datum, const ModelParams &params);
    void update_sc_cache(Node<S,P> *curr_node, Node<S,P> *new_node, size_t mut_id, const ModelParams &params);
    void print_cache();

    //double compute_loglik_sc(unordered_set<const S *> &snvs, C *cell, const ModelParams &params);
    double compute_loglik_sc(unordered_set<const S *> &snvs, size_t cell_id);
    
    EigenMatrix sc_presence_matrix_, sc_absence_matrix_;
    void PreComputeScLikelihood(const ModelParams &model_params);

public:
    TSSBState(const gsl_rng *random,
              Node<S,P> *root,
              const ModelParams &params,
              double (*log_lik_datum)(const Node<S,P> *node, const S *datum, const ModelParams &params),
              double (*log_lik_sc_at_site)(const S *s, const C *c, bool has_snv, const ModelParams &params),
              vector<S *> *bulk_data,
              vector<C *> *sc_data);

    void insert_datum(const gsl_rng *random, S* datum, const ModelParams &params);
    void set_sc_data(vector<C *> *sc_data, const ModelParams &model_params);
    
    const vector<S*> &get_data() const;

    // getters
    size_t get_num_nodes();
    inline Node<S,P> *get_root() { return root; }
    Node<S,P> *get_node(const S *datum);
    static void get_mixture(Node<S,P> *node, double mass, vector<pair<double, Node<S,P> *> > &ret);
    static void get_all_nodes(bool non_empty, Node<S,P> *root_node, vector<Node<S,P> *> &ret);
    void get_all_nodes(bool non_empty, vector<Node<S,P> *> &ret) const;
    void get_all_nodes(vector<Node<S,P> *> &ret) const;
    inline bool is_root(const Node<S,P> *node) const { return root == node; }

    // operations relating to sticks
    void update_sticks(const gsl_rng *random, const ModelParams &params);
    pair<size_t, size_t> descend_and_sample_sticks(const gsl_rng *random, Node<S,P> *node, const ModelParams &params);
    void reorder_sticks(const gsl_rng *random, const ModelParams &model_params);

    void resample_data_assignment(const gsl_rng *random, const ModelParams &params);
    //void move_datum(Node<S,P> *node, S *datum, const ModelParams &model_params);
    void move_datum(Node<S,P> *new_node, size_t mut_id, const ModelParams &model_params);
    
    double compute_log_likelihood_bulk(const ModelParams &params);
    double compute_log_likelihood_sc_cached(const ModelParams &params, bool verbose=false);

    double get_log_lik();
    static double get_log_prior_assignment(Node<S,P> *root);
    inline double get_log_lik_bulk() { return log_lik_bulk; }
    inline double get_log_lik_sc() { return log_lik_sc; }
    //double compute_log_likelihood_sc_subset(const gsl_rng *random, const ModelParams &params);

    // for output and debugging
    string print();
    
    // static functions
    static gsl_matrix *get_ancestral_matrix(TSSBState<S,C,P> &state);
    static TSSBState<S,C,P> *construct_trivial_state(Node<S,P> *root,
                                                     double (*log_lik_datum)(const Node<S,P> *v, const S *s, const ModelParams &params),
                                                     double (*log_lik_sc_at_site)(const S *s, const C *c, bool has_snv, const ModelParams &params));
    static double compute_loglik_sc(const vector<S *> &bulk_data,
                                    unordered_set<const S *> &snvs,
                                    C *cell,
                                    const ModelParams &params,
                                    double (*log_lik_sc_at_site)(const S *s, const C *c, bool has_snv, const ModelParams &params));
    static double compute_log_likelihood_sc(Node<S,P> *root,
                                            const vector<S *> &bulk_data,
                                            const vector<C *> &sc_data,
                                            double (*log_lik_sc_at_site)(const S *s, const C *c, bool has_snv, const ModelParams &params),
                                            const ModelParams &params,
                                            bool verbose=false);
    static void sample_alpha0(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<Node<S,P> *> &nodes);
    static void sample_lambda(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<Node<S,P> *> &nodes);
    static void sample_gamma(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<double> psi_sticks);
    
    static void update_hyper_params(const gsl_rng *random,
                                     size_t n_mh_iter,
                                     TSSBState<S,C,P> &state,
                                     ModelParams &params);

    friend size_t descend_and_cull(Node<S,P> *node);
    friend void descend_and_update_names(Node<S,P> *node);
    void clear_cache();
};

template <class S, class C, class P>
TSSBState<S,C,P>::TSSBState(const gsl_rng *random,
                            Node<S,P> *root,
                            const ModelParams &params,
                            double (*log_lik_datum)(const Node<S,P> *v, const S *s, const ModelParams &params),
                            double (*log_lik_sc_at_site)(const S *s, const C *c, bool has_snv, const ModelParams &params),
                            vector<S*> *bulk_data,
                            vector<C*> *sc_data) :
root(root), bulk_data(bulk_data), sc_data(sc_data), log_lik_datum(log_lik_datum), log_lik_sc_at_site(log_lik_sc_at_site)
{
    if (sc_data != 0) {
      sc_cache.resize(sc_data->size());
    }

    // Pre-compute single cell likelihood.
    PreComputeScLikelihood(params);

    root->sample_node_parameters(random, params, 0);
    for (size_t idx = 0; idx < bulk_data->size(); idx++) {
        S *datum = (*bulk_data)[idx];
        //this->initialize_data_assignment(random, datum, params);
        this->initialize_data_assignment(random, idx, params);
    }
    
    // initialize the log likelihoods
    log_lik_bulk = compute_log_likelihood_bulk(params);
    cout << "Initializing cache..." << endl;
    log_lik_sc = compute_log_likelihood_sc_cached(params);

//    double sanity_check = compute_log_likelihood_sc(root, *bulk_data, *sc_data, log_lik_sc_at_site, params);
//    assert(abs(sanity_check - log_lik_sc) < 1e-3);
}

template <class S, class C, class P>
void TSSBState<S,C,P>::PreComputeScLikelihood(const ModelParams &model_params) {
    size_t cell_count = sc_data->size();
    size_t mutation_count = bulk_data->size();
    sc_presence_matrix_ = EigenMatrix::Zero(cell_count, mutation_count);
    sc_absence_matrix_ = EigenMatrix::Zero(cell_count, mutation_count);
    // Evaluate single cell log likelihoods.
    for (size_t c = 0; c < cell_count; c++) {
        cout << "Cell name: " << sc_data->at(c)->get_name() << endl;
        for (size_t n = 0; n < mutation_count; n++) {
            sc_presence_matrix_(c,n) = log_lik_sc_at_site(bulk_data->at(n), sc_data->at(c), true, model_params);
            sc_absence_matrix_(c,n) = log_lik_sc_at_site(bulk_data->at(n), sc_data->at(c), false, model_params);
        }
    }
}

// the nu stick for the root node will have been sampled or set to 0 depending on the use case
template <class S, class C, class P>
TSSBState<S,C,P>::TSSBState(Node<S,P> *root)
{
    sc_data = 0;
    bulk_data = new vector<S*>();
    this->root = root;
}

/*******************
 private functions
 *******************/
//template <class S, class C, class P>
//void TSSBState<S,C,P>::initialize_data_assignment(const gsl_rng *random, S *datum, const ModelParams &params)
//{
//    // sample an initial tree
////    double u = uniform(random, 0.0, 1.0);
////    Node<S,P> *node = Node<S,P>::find_node(random, u, root, params); // this call may generate new nodes
////    assign_data_point(0, node, datum, params);
//
//    // assign to the root
//    //assign_data_point(0, root, datum, params);
//
//    // sample once then, assign to first child of root
//    // root denotes healthy population and exist only for convenience regarding parameter sampling/computation
//    // as none of the SNVs should be assigned to the root, create a single child of this root and initialize the tree by
//    // assigning all of the SNVs there
//    if (root->get_num_children() == 0) {
//        // Spawn one child.
//        root->InitializeChild(random, params);
//        // sample an initial tree
//        //double u = uniform(random, 0.0, 1.0);
//        //Node<S,P>::find_node(random, u, root, params); // this call may generate new nodes
//    }
//
//    // assign the datum to the first child of root
//    assert(root->get_num_children() > 0);
//    Node<S,P> *first_child_node = root->get_child(0).second;
//    assign_data_point(0, first_child_node, datum, params);
//}

template <class S, class C, class P>
void TSSBState<S,C,P>::initialize_data_assignment(const gsl_rng *random, size_t mut_id, const ModelParams &params)
{
    if (root->get_num_children() == 0) {
        // Spawn one child.
        root->InitializeChild(random, params);
    }
    
    // assign the datum to the first child of root
    assert(root->get_num_children() > 0);
    Node<S,P> *first_child_node = root->get_child(0).second;
    assign_data_point(0, first_child_node, mut_id, params);
}


template <class S, class C, class P>
void TSSBState<S,C,P>::initialize_sc_cache(size_t c, Node<S,P> *v, unordered_set<const S*> &snvs, const ModelParams &model_params)
{
    unordered_map<Node<S,P> *, double> &cache_c = sc_cache[c];
    //cache_c[v] = v->compute_log_likelihood_of_sc(sc_data->at(c), snas, model_params);
    //cache_c[v] = compute_loglik_sc(snvs, sc_data->at(c), model_params);
    cache_c[v] = compute_loglik_sc(snvs, c);
}

//template <class S, class C, class P>
//void TSSBState<S,C,P>::update_sc_cache(Node<S,P> *curr_node, Node<S,P> *new_node, S *datum, const ModelParams &params)
//{
//    //cout << "===== Update sc cache =====" << endl;
//    if (sc_data == 0) {
//        return;
//    }
//
//    vector<Node<S,P> *> subtree_curr;
//    vector<Node<S,P> *> subtree_new;
//    get_all_nodes(false, curr_node, subtree_curr);
//    get_all_nodes(false, new_node, subtree_new);
//    unordered_set<Node<S,P> *> subtree_curr_set;
//    subtree_curr_set.insert(subtree_curr.begin(), subtree_curr.end());
//    unordered_set<Node<S,P> *> subtree_new_set;
//    subtree_new_set.insert(subtree_new.begin(), subtree_new.end());
//
//    // retrieve set of SNVs for each node in subtree_curr and subtree_new
//    unordered_map<Node<S,P> *, unordered_set<const S *> *> node2snvs;
//    for (Node<S,P> *v : subtree_curr)
//    {
//        unordered_set<const S *> *snvs = new unordered_set<const S *>();
//        Node<S,P>::get_dataset(v, *snvs);
//        node2snvs[v] = snvs;
//    }
//    for (Node<S,P> *v : subtree_new)
//    {
//        if (node2snvs.count(v) == 0) {
//            unordered_set<const S *> *snvs = new unordered_set<const S *>();
//            Node<S,P>::get_dataset(v, *snvs);
//            node2snvs[v] = snvs;
//        }
//    }
//
//    bool exp_mut_status;
//    for (size_t c = 0; c < sc_data->size(); c++) {
//        // update the entries that belongs to subtrees of curr_node and new_node
//        // code is being repeated: facor this out
//        //cout << "Updating subtree of current node: " << curr_node->get_name() << endl;
//        for (Node<S,P> *v : subtree_curr) {
//            if (sc_cache[c].count(v) == 0) {
//                initialize_sc_cache(c, v, *node2snvs[v], params);
//                //cout << "New init: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
//                continue;
//            }
//
//            //cout << "Before: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
//            double x = (*log_lik_sc_at_site)(datum, sc_data->at(c), true, params);
//            sc_cache[c][v] -= x;
//            exp_mut_status = (subtree_new_set.count(v) > 0) ? true : false; // check if v is expected to have the mutation or not
//            double y = (*log_lik_sc_at_site)(datum, sc_data->at(c), exp_mut_status, params);
//            sc_cache[c][v] += y;
////            cout << "After: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
////            cout << "x: " << x << ", y: " << y << ", has mutation: " << exp_mut_status << endl;
//            // check sc_cache[c][v] is correct
////            double exp_val = compute_loglik_sc(*node2snvs[v], sc_data->at(c), params);
////            double sanity_check = sc_cache[c][v];
////            assert(abs(sanity_check - exp_val) < 1e-3);
////            if (abs(sanity_check - exp_val) > 1e-3) {
////                cout << "sanity check failed." << endl;
////            }
//
//        }
//
//        //cout << "Updating subtree of new node: " << new_node->get_name() << endl;
//        for (Node<S,P> *v : subtree_new) {
//            if (sc_cache[c].count(v) == 0) {
//                initialize_sc_cache(c, v, *node2snvs[v], params);
//                //cout << "New init: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
//                continue;
//            }
//            exp_mut_status = (subtree_curr_set.count(v) > 0) ? true : false; // check if v had the mutation or not
//            //cout << "Before: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
//            double x = (*log_lik_sc_at_site)(datum, sc_data->at(c), exp_mut_status, params);
//            sc_cache[c][v] -= x;
//            double y = (*log_lik_sc_at_site)(datum, sc_data->at(c), true, params);
//            sc_cache[c][v] += y;
////            cout << "After: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
////            cout << "x: " << x << ", y: " << y << ", has mutation: " << exp_mut_status << endl;
//            // check sc_cache[c][v] is correct
////            double exp_val = compute_loglik_sc(*node2snvs[v], sc_data->at(c), params);
////            double sanity_check = sc_cache[c][v];
////            assert(abs(sanity_check - exp_val) < 1e-3);
////            if (abs(sanity_check - exp_val) > 1e-3) {
////                cout << "sanity check failed." << endl;
////            }
//        }
//    }
//
//    for (auto it = node2snvs.begin(); it != node2snvs.end(); ++it) {
//        delete it->second;
//    }
//}

template <class S, class C, class P>
void TSSBState<S,C,P>::update_sc_cache(Node<S,P> *curr_node, Node<S,P> *new_node, size_t mut_id, const ModelParams &params)
{
    //cout << "===== Update sc cache =====" << endl;
    if (sc_data == 0) {
        return;
    }

    vector<Node<S,P> *> subtree_curr;
    vector<Node<S,P> *> subtree_new;
    get_all_nodes(false, curr_node, subtree_curr);
    get_all_nodes(false, new_node, subtree_new);
    unordered_set<Node<S,P> *> subtree_curr_set;
    subtree_curr_set.insert(subtree_curr.begin(), subtree_curr.end());
    unordered_set<Node<S,P> *> subtree_new_set;
    subtree_new_set.insert(subtree_new.begin(), subtree_new.end());

    // retrieve set of SNVs for each node in subtree_curr and subtree_new
    unordered_map<Node<S,P> *, unordered_set<const S *> *> node2snvs;
    for (Node<S,P> *v : subtree_curr)
    {
        unordered_set<const S *> *snvs = new unordered_set<const S *>();
        // This loops until it reaches the root.
        // We can write better code for getting the mutations without traversing to the root each time.
        Node<S,P>::get_dataset(v, *snvs);
        node2snvs[v] = snvs;
    }
    for (Node<S,P> *v : subtree_new)
    {
        if (node2snvs.count(v) == 0) {
            unordered_set<const S *> *snvs = new unordered_set<const S *>();
            Node<S,P>::get_dataset(v, *snvs);
            node2snvs[v] = snvs;
        }
    }

    bool exp_mut_status;
    for (size_t c = 0; c < sc_data->size(); c++) {
        
        // update the entries that belong to subtrees of curr_node and new_node
        // TODO: code is being repeated: facor this out
        //cout << "Updating subtree of current node: " << curr_node->get_name() << endl;
        for (Node<S,P> *v : subtree_curr) {
            if (sc_cache[c].count(v) == 0) {
                initialize_sc_cache(c, v, *node2snvs[v], params);
                //cout << "New init: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
                continue;
            }
            
            //cout << "Before: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
            //double x = (*log_lik_sc_at_site)(datum, sc_data->at(c), true, params);
            double x = sc_presence_matrix_(c, mut_id);
            sc_cache[c][v] -= x;
            exp_mut_status = (subtree_new_set.count(v) > 0) ? true : false; // check if v is expected to have the mutation or not
            //double y = (*log_lik_sc_at_site)(datum, sc_data->at(c), exp_mut_status, params);
            double y = exp_mut_status ? sc_presence_matrix_(c, mut_id) : sc_absence_matrix_(c, mut_id);
            sc_cache[c][v] += y;
            //            cout << "After: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
            //            cout << "x: " << x << ", y: " << y << ", has mutation: " << exp_mut_status << endl;
            // check sc_cache[c][v] is correct
            //            double exp_val = compute_loglik_sc(*node2snvs[v], sc_data->at(c), params);
            //            double sanity_check = sc_cache[c][v];
            //            assert(abs(sanity_check - exp_val) < 1e-3);
            //            if (abs(sanity_check - exp_val) > 1e-3) {
            //                cout << "sanity check failed." << endl;
            //            }
            
        }
        
        //cout << "Updating subtree of new node: " << new_node->get_name() << endl;
        for (Node<S,P> *v : subtree_new) {
            if (sc_cache[c].count(v) == 0) {
                initialize_sc_cache(c, v, *node2snvs[v], params);
                //cout << "New init: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
                continue;
            }
            exp_mut_status = (subtree_curr_set.count(v) > 0) ? true : false; // check if v had the mutation or not
            //cout << "Before: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
            //double x = (*log_lik_sc_at_site)(datum, sc_data->at(c), exp_mut_status, params);
            double x = exp_mut_status ? sc_presence_matrix_(c, mut_id) : sc_absence_matrix_(c, mut_id);
            sc_cache[c][v] -= x;
            //double y = (*log_lik_sc_at_site)(datum, sc_data->at(c), true, params);
            double y = sc_presence_matrix_(c, mut_id);
            sc_cache[c][v] += y;
            //            cout << "After: l[" << c << ", " << v << "] = " << sc_cache[c][v] << endl;
            //            cout << "x: " << x << ", y: " << y << ", has mutation: " << exp_mut_status << endl;
            // check sc_cache[c][v] is correct
            //            double exp_val = compute_loglik_sc(*node2snvs[v], sc_data->at(c), params);
            //            double sanity_check = sc_cache[c][v];
            //            assert(abs(sanity_check - exp_val) < 1e-3);
            //            if (abs(sanity_check - exp_val) > 1e-3) {
            //                cout << "sanity check failed." << endl;
            //            }
        }
    }
    
    for (auto it = node2snvs.begin(); it != node2snvs.end(); ++it) {
        delete it->second;
    }
}


template <class S, class C, class P>
void TSSBState<S,C,P>::print_cache()
{
    // print cache
    //cout << "*****In sc cache.*****" << endl;
    vector<Node<S,P> *> all_nodes;
    get_all_nodes(false, all_nodes);
    for (size_t c = 0; c < sc_data->size(); c++) {
        for (size_t i = 0; i < all_nodes.size(); i++) {
            Node<S,P> *node = all_nodes[i];
            cout << "l[" << c << ", " << node << ", " << node->get_name() <<  "]=" << sc_cache[c][node] << endl;
        }
    }
}

//template <class S, class C, class P>
//void TSSBState<S,C,P>::slice_sample_data_assignment(const gsl_rng *random, S *datum, const ModelParams &model_params)
//{
//    double u_min = 0.0, u_max = 1.0;
//    double u;
//    Node<S,P> *curr_node = datum2node[datum];
//    Node<S,P> *new_node = 0;
//
//    double curr_log_lik_bulk = log_lik_datum(curr_node, datum, model_params);
//    double curr_log_lik_sc = compute_log_likelihood_sc_cached(model_params);
////    double sanity_check = compute_log_likelihood_sc(root, *bulk_data, *sc_data, log_lik_sc_at_site, model_params);
////    assert(abs(curr_log_lik_sc - sanity_check) < 1e-3);
//    double curr_log_lik = curr_log_lik_bulk + curr_log_lik_sc;
//
//    double log_slice = log(uniform(random, 0.0, 1.0)); // sample the slice
//    double new_log_lik_bulk = 0.0, new_log_lik_sc = 0.0, new_log_lik = 0.0;
//    size_t iter = 0;
//
//    while ((abs(u_max - u_min) > 1e-3)) {
//        //cout << print() << endl;
//
//        iter++;
//        u = uniform(random, u_min, u_max);
//        new_node = Node<S,P>::find_node(random, u, root, model_params); // this call may generate new nodes
//
//        if (new_node == curr_node) {
//            break; // no need to carry out further computation
//        }
//
//        // assign data point from curr_node to new_node
//        assign_data_point(curr_node, new_node, datum, model_params);
//
//        // compute the log likelihood
//        new_log_lik_bulk = log_lik_datum(new_node, datum, model_params);
//        new_log_lik_sc = compute_log_likelihood_sc_cached(model_params);
////        sanity_check = compute_log_likelihood_sc(root, *bulk_data, *sc_data, log_lik_sc_at_site, model_params);
////        assert(abs(new_log_lik_sc - sanity_check) < 1e-3);
////        if (abs(new_log_lik_sc - sanity_check) > 1e-3)
////        {
////            cout << print() << endl;
////            new_log_lik_sc = compute_log_likelihood_sc_cached(model_params, true);
////            sanity_check = compute_log_likelihood_sc(root, *bulk_data, *sc_data, log_lik_sc_at_site, model_params, true);
////            unordered_set<const S*> snvs;
////            Node<S,P>::get_dataset(new_node, snvs);
////            for (const S* snv : snvs)
////            {
////                cout << snv << endl;
////            }
////        }
//        new_log_lik = new_log_lik_bulk + new_log_lik_sc;
//        if (new_log_lik > (log_slice + curr_log_lik)) {
//            // update log_lik_bulk, log_lik_sc, log_lik
//            log_lik_bulk = compute_log_likelihood_bulk(model_params);
//            log_lik_sc = new_log_lik_sc;
//            log_lik = log_lik_bulk + log_lik_sc;
//            break;
//        } else {
//            // revert the changes
//            assign_data_point(new_node, curr_node, datum, model_params);
//            new_log_lik_sc = compute_log_likelihood_sc_cached(model_params);
////            sanity_check = compute_log_likelihood_sc(root, *bulk_data, *sc_data, log_lik_sc_at_site, model_params);
////            assert(abs(new_log_lik_sc - sanity_check) < 1e-3);
//            // this value may not equal to curr_log_lik_sc because new nodes may have been created
//            // which could alter the sc log likelihood, same goes for bulk
//            if (Node<S,P>::less(new_node, curr_node)) {
//                u_min = u;
//            } else {
//                u_max = u;
//            }
//            if (abs(u_max - u_min) < 1e-6) {
//                cout << "Slice sampler shrank!" << endl;
//                break;
//            }
//        }
//
//        if (iter == 20) { // slice sampler can degenerate, keep the max iteration to 20
//            break;
//        }
//    }
//
//    // note: if datum is not assigned to a new node (reverted), then
//    // any new node created does not modify previous bulk log likelihood nor the single cells
//
//    if (abs(u_max - u_min) <= 1e-3) {
//        // slice sampler degenerated: this occurs far too often, wasteful
//        //cout << "Slice sampler degenerated. TODO: figure this out at some point!" << endl;
//    }
//}

template <class S, class C, class P>
void TSSBState<S,C,P>::slice_sample_data_assignment(const gsl_rng *random, size_t mut_id, const ModelParams &model_params)
{
    auto datum = bulk_data->at(mut_id);
    double u_min = 0.0, u_max = 1.0;
    double u;
    Node<S,P> *curr_node = datum2node[datum];
    Node<S,P> *new_node = 0;
    
    double curr_log_lik_bulk = log_lik_datum(curr_node, datum, model_params);
    double curr_log_lik_sc = compute_log_likelihood_sc_cached(model_params);
    //    double sanity_check = compute_log_likelihood_sc(root, *bulk_data, *sc_data, log_lik_sc_at_site, model_params);
    //    assert(abs(curr_log_lik_sc - sanity_check) < 1e-3);
    double curr_log_lik = curr_log_lik_bulk + curr_log_lik_sc;
    
    double log_slice = log(uniform(random, 0.0, 1.0)); // sample the slice
    double new_log_lik_bulk = 0.0, new_log_lik_sc = 0.0, new_log_lik = 0.0;
    size_t iter = 0;
    
    while ((abs(u_max - u_min) > 1e-3)) {
        //cout << print() << endl;
        
        iter++;
        u = uniform(random, u_min, u_max);
        new_node = Node<S,P>::find_node(random, u, root, model_params); // this call may generate new nodes
        
        if (new_node == curr_node) {
            break; // no need to carry out further computation
        }
        
        // assign data point from curr_node to new_node
        assign_data_point(curr_node, new_node, mut_id, model_params);
        
        // compute the log likelihood
        new_log_lik_bulk = log_lik_datum(new_node, datum, model_params);
        new_log_lik_sc = compute_log_likelihood_sc_cached(model_params);
        //        sanity_check = compute_log_likelihood_sc(root, *bulk_data, *sc_data, log_lik_sc_at_site, model_params);
        //        assert(abs(new_log_lik_sc - sanity_check) < 1e-3);
        //        if (abs(new_log_lik_sc - sanity_check) > 1e-3)
        //        {
        //            cout << print() << endl;
        //            new_log_lik_sc = compute_log_likelihood_sc_cached(model_params, true);
        //            sanity_check = compute_log_likelihood_sc(root, *bulk_data, *sc_data, log_lik_sc_at_site, model_params, true);
        //            unordered_set<const S*> snvs;
        //            Node<S,P>::get_dataset(new_node, snvs);
        //            for (const S* snv : snvs)
        //            {
        //                cout << snv << endl;
        //            }
        //        }
        new_log_lik = new_log_lik_bulk + new_log_lik_sc;
        if (new_log_lik > (log_slice + curr_log_lik)) {
            // update log_lik_bulk, log_lik_sc, log_lik
            log_lik_bulk = compute_log_likelihood_bulk(model_params);
            log_lik_sc = new_log_lik_sc;
            log_lik = log_lik_bulk + log_lik_sc;
            break;
        } else {
            // revert the changes
            assign_data_point(new_node, curr_node, mut_id, model_params);
            new_log_lik_sc = compute_log_likelihood_sc_cached(model_params);
            //            sanity_check = compute_log_likelihood_sc(root, *bulk_data, *sc_data, log_lik_sc_at_site, model_params);
            //            assert(abs(new_log_lik_sc - sanity_check) < 1e-3);
            // this value may not equal to curr_log_lik_sc because new nodes may have been created
            // which could alter the sc log likelihood, same goes for bulk
            if (Node<S,P>::less(new_node, curr_node)) {
                u_min = u;
            } else {
                u_max = u;
            }
            if (abs(u_max - u_min) < 1e-6) {
                cout << "Slice sampler shrank!" << endl;
                break;
            }
        }
        
        if (iter == 20) { // slice sampler can degenerate, keep the max iteration to 20
            break;
        }
    }
    
    // note: if datum is not assigned to a new node (reverted), then
    // any new node created does not modify previous bulk log likelihood nor the single cells
    
    if (abs(u_max - u_min) <= 1e-3) {
        // Slice sampler degenerated: this seems to occur quite a bit.
        //cout << "Slice sampler degenerated. TODO: figure this out at some point!" << endl;
    }
}


template <class S, class C, class P>
const vector<S*> &TSSBState<S,C,P>::get_data() const
{
    return *bulk_data;
}

//template <class S, class C, class P>
//void TSSBState<S,C,P>::assign_data_point(Node<S,P> *curr_node, Node<S,P> *new_node, S *datum, const ModelParams &model_params)
//{
//    if (curr_node != 0) {
//        curr_node->remove_datum(datum);
//    }
//
//    // update the map : datum -> node
//    datum2node[datum] = new_node;
//    // add the datum to the new node
//    new_node->add_datum(datum);
//
//    // update single cell cache
//    if (sc_data != 0 && sc_data->size() > 0)
//        update_sc_cache(curr_node, new_node, datum, model_params);
//}

template <class S, class C, class P>
void TSSBState<S,C,P>::assign_data_point(Node<S,P> *curr_node, Node<S,P> *new_node, size_t mut_id, const ModelParams &model_params)
{
    auto datum = bulk_data->at(mut_id);
    if (curr_node != 0) {
        curr_node->remove_datum(datum);
    }

    // update the map : datum -> node
    datum2node[datum] = new_node;
    // add the datum to the new node
    new_node->add_datum(datum);
    
    // update single cell cache.
    // TODO: do it only if mut_id has single cell coverage.
    if (sc_data != 0 && sc_data->size() > 0)
        update_sc_cache(curr_node, new_node, mut_id, model_params);
}


/**********
 public functions
 **********/
template <class S, class C, class P>
double TSSBState<S,C,P>::get_log_prior_assignment(Node<S,P> *root)
{
    // log prior of the data assignment requires computing the probability of
    // assignment to a node and the number of data points at that node.
    vector<pair<double, Node<S,P> *> > mixture;
    get_mixture(root, 1.0, mixture);
    
    double log_prior_assignment = 0.0;
    for (size_t i = 0; i < mixture.size(); i++) {
        size_t n_data = mixture[i].second->get_num_data();
        if (n_data > 0) {
            log_prior_assignment += n_data * log(mixture[i].first);
        }
    }
    return log_prior_assignment;
}

template <class S, class C, class P>
double TSSBState<S,C,P>::get_log_lik()
{
    double log_prior_assignment = get_log_prior_assignment(get_root());
    log_lik = log_lik_bulk + log_lik_sc + log_prior_assignment;
    return log_lik;
}

template <class S, class C, class P>
double TSSBState<S,C,P>::compute_log_likelihood_sc(Node<S,P> *root,
                                                   const vector<S *> &bulk_data,
                                                   const vector<C *> &sc_data,
                                                   double (*log_lik_sc_at_site)(const S *s, const C *c, bool has_snv, const ModelParams &params),
                                                   const ModelParams &params,
                                                   bool verbose)
{
    if (sc_data.size() == 0) {
        return 0.0;
    }
    
    double log_lik_sc = 0.0;
    
    // marginalize over assignment to the nodes
    vector<Node<S,P> *> all_nodes;
    get_all_nodes(true, root, all_nodes);
    unordered_map<Node<S,P> *, unordered_set<const S *> > node2snvs;
    for (Node<S,P> *v : all_nodes)
    {
        unordered_set<const S *> snvs;
        Node<S,P>::get_dataset(v, snvs);
        node2snvs[v] = snvs;
    }
    
    double log_prior_assignment = -log(all_nodes.size());
    for (size_t c = 0; c < sc_data.size(); c++) {
        double log_lik_cell = DOUBLE_NEG_INF;
        C *cell = sc_data.at(c);
        for (size_t i = 0; i < all_nodes.size(); i++) {
            Node<S,P> *node = all_nodes[i];
            auto snvs = node2snvs.at(node);
            double log_val = TSSBState<S,C,P>::compute_loglik_sc(bulk_data, snvs, cell, params, log_lik_sc_at_site);
            if (verbose) {
                cout << "[" << node->get_name() << "]: " << log_val << endl;
            }
            log_lik_cell = log_add(log_val, log_lik_cell);
        }
        log_lik_cell += log_prior_assignment;
        log_lik_sc += log_lik_cell;
        //cout << "Cell " << c << ": " << log_lik_cell << endl;
    }
    
    //    for (auto it = node2snvs.begin(); it != node2snvs.end(); ++it) {
    //        delete it->second;
    //    }
    
    return log_lik_sc;
    
    //    // get all nodes -- find SNVs for each node
    //    vector<Node<S,P> *> all_nodes;
    //    get_all_nodes(false, all_nodes);
    //    double log_prior_assignment = -log(all_nodes.size());
    //    unordered_map<Node<S,P> *, unordered_set<S*> > node2snas;
    //    for (size_t i = 0; i < all_nodes.size(); i++) {
    //        Node<S,P> *node = all_nodes[i];
    //        unordered_set<S*> parent_snas = node2snas[node->get_parent_node()];
    //        unordered_set<S*> node_snas = node->get_data();
    //        node_snas.insert(parent_snas.begin(), parent_snas.end());
    //        node2snas[node] = node_snas;
    //    }
    //
    //    // get all non-empty nodes, these are possible candidates for assignment of single cells
    //    //vector<Node<S,P> *> non_empty_nodes;
    //    //get_all_nodes(false, non_empty_nodes);
    //
    ////    cout << "*****In compute sc log lik.*****" << endl;
    //    double log_lik_sc = 0.0;
    //    // TODO: this can be parallelized over cells
    //    for (size_t c = 0; c < sc_data->size(); c++) {
    //        // marginalize over the nodes
    //        double log_lik_cell = DOUBLE_NEG_INF;
    //        for (size_t i = 0; i < all_nodes.size(); i++) {
    //            Node<S,P> *node = all_nodes[i];
    //            double log_val = node->compute_log_likelihood_of_sc(sc_data->at(c), node2snas[node], params);
    //            log_lik_cell = log_add(log_val, log_lik_cell);
    //        }
    //        log_lik_cell += log_prior_assignment;
    //        log_lik_sc += log_lik_cell;
    //    }
    //    return log_lik_sc;
}

template <class S, class C, class P>
double TSSBState<S,C,P>::compute_log_likelihood_sc_cached(const ModelParams &params, bool verbose)
{
    if (sc_data == 0 || sc_data->size() == 0) {
        return 0.0;
    }

    // get all non-empty nodes, these are possible candidates for assignment of single cells
    vector<Node<S,P> *> all_nodes;
    get_all_nodes(true, all_nodes);
    unordered_map<Node<S,P> *, unordered_set<const S *> *> node2snvs;
    // TODO: This can be improved by passing down SNVs.
    for (Node<S,P> *v : all_nodes)
    {
        unordered_set<const S *> *snvs = new unordered_set<const S *>();
        Node<S,P>::get_dataset(v, *snvs);
        node2snvs[v] = snvs;
    }

    double log_prior_assignment = -log(all_nodes.size());

    double log_lik_sc = 0.0;
    double log_val;
    // TODO: this can be parallelized over cells
    for (size_t c = 0; c < sc_data->size(); c++) {
        // marginalize over the nodes
        unordered_map<Node<S,P> *, double> &cache_c = sc_cache[c];
        double log_lik_cell = DOUBLE_NEG_INF;
        for (size_t i = 0; i < all_nodes.size(); i++) {
            Node<S,P> *v = all_nodes[i];
            if (cache_c.count(v) == 0) {
                initialize_sc_cache(c, v, *node2snvs[v], params);
            }
            log_val = cache_c[v];
            if (verbose) {
                cout << "[" << v->get_name() << "]: " << log_val << endl;
            }
            log_lik_cell = log_add(log_val, log_lik_cell);
        }
        log_val = (log_lik_cell + log_prior_assignment);
        //cout << "Cell " << c << ": " << log_val << endl;
        log_lik_sc += log_val;
    }

    for (auto it = node2snvs.begin(); it != node2snvs.end(); ++it)
    {
        delete it->second;
    }

    //print_cache();
    return log_lik_sc;
}

//template <class S, class C, class P>
//double TSSBState<S,C,P>::compute_loglik_sc(unordered_set<const S *> &snvs, C *cell, const ModelParams &params)
//{
//    double log_lik = 0.0;
//    bool has_snv;
//    for (size_t i = 0; i < bulk_data->size(); i++) {
//        const S *snv = bulk_data->at(i);
//        has_snv = (snvs.count(snv) > 0);
//        double log_val = log_lik_sc_at_site(snv, cell, has_snv, params);
//        log_lik += log_val;
//        //cout << "Has SNV: " << has_snv << ", " << log_val << endl;
//    }
//
//    return log_lik;
//}

template <class S, class C, class P>
double TSSBState<S,C,P>::compute_loglik_sc(unordered_set<const S *> &snvs, size_t cell_id)
{
    double log_lik = 0.0;
    bool has_snv;
    // TODO: change this loop can be made faster by looping over only the sites with reads.
    for (size_t i = 0; i < bulk_data->size(); i++) {
        const S *snv = bulk_data->at(i);
        has_snv = (snvs.count(snv) > 0);
        double log_val = has_snv ? sc_presence_matrix_(cell_id, i) :
                            sc_absence_matrix_(cell_id, i);
        log_lik += log_val;
        //cout << "Has SNV: " << has_snv << ", " << log_val << endl;
    }
    return log_lik;
}


template <class S, class C, class P>
double TSSBState<S,C,P>::compute_loglik_sc(const vector<S *> &bulk_data,
                                           unordered_set<const S *> &snvs,
                                           C *cell,
                                           const ModelParams &params,
                                           double (*log_lik_sc_at_site)(const S *s, const C *c, bool has_snv, const ModelParams &params))
{
    double log_lik = 0.0;
    bool has_snv;
    for (size_t i = 0; i < bulk_data.size(); i++) {
        const S *snv = bulk_data.at(i);
        has_snv = (snvs.count(snv) > 0);
        log_lik += log_lik_sc_at_site(snv, cell, has_snv, params);
    }

    return log_lik;
}

//template <class S, class C, class P>
//double TSSBState<S,C,P>::compute_log_likelihood_sc_subset(const gsl_rng *random, const ModelParams &params)
//{
//    if (sc_data == 0 || sc_data->size() == 0) {
//        return 0.0;
//    }
//
//    // get all nodes -- find SNAs for each node
//    vector<Node<S,P> *> all_nodes;
//    get_all_nodes(false, all_nodes);
//    unordered_map<Node<S,P> *, unordered_set<S*> > node2snas;
//    for (size_t i = 0; i < all_nodes.size(); i++) {
//        Node<S,P> *node = all_nodes[i];
//        unordered_set<S*> parent_snas = node2snas[node->get_parent_node()];
//        unordered_set<S*> node_snas = node->get_data();
//        node_snas.insert(parent_snas.begin(), parent_snas.end());
//        node2snas[node] = node_snas;
//    }
//
//    // get all non-empty nodes, these are possible candidates for assignment of single cells
//    vector<Node<S,P> *> non_empty_nodes;
//    get_all_nodes(true, non_empty_nodes);
//
//    double log_lik_sc = 0.0;
//    // select subset
//    size_t M = 50;
//    M = sc_data->size() < M ? sc_data->size() : M;
//    size_t *n = new size_t[sc_data->size()];
//    for (size_t c = 0; c < sc_data->size(); c++) {
//        n[c] = c;
//    }
//    gsl_ran_shuffle(random, n, sc_data->size(), sizeof(size_t));
//    for (size_t c = 0; c < M; c++) {
//        size_t cell = n[c];
//        // marginalize over the nodes
//        double log_lik_cell = DOUBLE_NEG_INF;
//        for (size_t i = 0; i < non_empty_nodes.size(); i++) {
//            Node<S,P> *node = non_empty_nodes[i];
//            double log_val = node->compute_log_likelihood_of_sc(sc_data->at(cell), node2snas[node], params);
//            log_lik_cell = log_add(log_val, log_lik_cell);
//        }
//
//        log_lik_sc += log_lik_cell;
//    }
//
//    delete [] n;
//
//    return log_lik_sc;
//}

template <class S, class C, class P>
void TSSBState<S,C,P>::resample_data_assignment(const gsl_rng *random, const ModelParams &params)
{
    for (size_t mut_id = 0; mut_id < bulk_data->size(); mut_id++)
    {
        //slice_sample_data_assignment(random, bulk_data->at(i), params);
        slice_sample_data_assignment(random, mut_id, params);
    }
//    size_t i = gsl_rng_uniform_int(random, bulk_data->size());
    //const S *s = bulk_data->at(i);
    //cout << "Assigning: " << s << endl;
//    slice_sample_data_assignment(random, bulk_data->at(i), params);
}

//template <class S, class C, class P>
//void TSSBState<S,C,P>::move_datum(Node<S,P> *new_node, S *datum, const ModelParams &model_params)
//{
//    // should check if this node is part of the tree
//    Node<S,P> *curr_node = new_node;
//    // if curr_node has the root as ancestor, it is a valid node
//    while (curr_node != root) {
//        if (curr_node == 0) {
//            cerr << "Error: Node is not part of the tree!" << endl;
//        }
//        curr_node = curr_node->get_parent_node();
//    }
//
//    curr_node = datum2node[datum];
//    assign_data_point(curr_node, new_node, datum, model_params);
//}

template <class S, class C, class P>
void TSSBState<S,C,P>::move_datum(Node<S,P> *new_node, size_t mut_id, const ModelParams &model_params)
{
    // should check if this node is part of the tree
    Node<S,P> *curr_node = new_node;
    // if curr_node has the root as ancestor, it is a valid node
    while (curr_node != root) {
        if (curr_node == 0) {
            cerr << "Error: Node is not part of the tree!" << endl;
        }
        curr_node = curr_node->get_parent_node();
    }
    
    curr_node = datum2node[bulk_data->at(mut_id)];
    assign_data_point(curr_node, new_node, mut_id, model_params);
}


// helper functions for updating the sticks
template <class S, class C, class P>
void TSSBState<S,C,P>::reorder_sticks(const gsl_rng *random, const ModelParams &model_params)
{
    queue<Node<S,P> *> q;
    q.push(root);
    while (!q.empty()) {
        Node<S,P> *node = q.front();
        q.pop();

        node->reorder_sticks(random, model_params);
        // enqueue children nodes to the queue
        for (auto it = node->get_idx2child().begin(); it != node->get_idx2child().end(); ++it) {
            q.push(it->second.second);
        }
    }
}

//template <class S, class C, class P>
//void TSSBState<S,C,P>::cull()
//{
//    descend_and_cull(root);
//    descend_and_update_names(root);
//}

template <class S, class P>
void descend_and_update_names(Node<S,P> *node)
//void TSSBState<S,C,P>::descend_and_update_names(Node<S,P> *node)
{
    node->reset_children_names();
    Node<S,P> *child = 0;
    
    unordered_map<size_t, pair<double, Node<S,P> *> > &children = node->get_idx2child();
    for (auto it = children.begin(); it != children.end(); it++) {
        child = it->second.second;
        descend_and_update_names(child);
    }
}


template <class S, class P>
size_t descend_and_cull(Node<S,P> *node)
{
    size_t n_data = node->get_num_data();
    size_t n_data_desc = 0;
    Node<S,P> *child = 0;

    unordered_map<size_t, pair<double, Node<S,P> *> > &children = node->get_idx2child();
    unordered_set<size_t> cull_list;
    bool stop = false;
    for (int i = children.size() - 1; i >= 0; i--) {

        child = children[i].second;

        size_t ret = descend_and_cull(child);
        n_data_desc += ret;

        if (ret > 0) {
            stop = true;
        } else if (ret == 0 && !stop) {
            // cull child node with index i
            cull_list.insert(i);
        }

    }

    node->cull(cull_list);

    return n_data + n_data_desc;
}

// returns pair<n_data at node, n_data at descendants of node>
template <class S, class C, class P>
pair<size_t, size_t> TSSBState<S,C,P>::descend_and_sample_sticks(const gsl_rng *random, Node<S,P> *node, const ModelParams &params)
{
    size_t n_data = node->get_num_data();
    size_t total_num_data_at_desc = 0; // number of data points below this node (over all descendant nodes)

    unordered_map<size_t, pair<double, Node<S,P> *> > &children = node->get_idx2child();
    vector<pair<size_t, size_t> > n_data_desc; // n_data_desc[i] stores <n_data at i, total_num_desc of i> (note: the second count does not include n_data at i)
    for (size_t i = 0; i < children.size(); i++)
    {
        Node<S,P> *child = children.at(i).second;
        pair<size_t, size_t> ret = descend_and_sample_sticks(random, child, params);
        n_data_desc.push_back(ret);
        total_num_data_at_desc += (ret.first + ret.second);
    }

    // update nu-stick except for the root, which stays at 0
    if (node == root && node->get_nu_stick() == 0) {
        // do not sample nu-stick for root
    } else {
        double temp = bounded_beta(random, n_data + 1, total_num_data_at_desc + params.alpha(node->get_name_vec()));
        node->set_nu_stick(temp);
    }

    // update psi-sticks
    size_t cumulative_num_data = 0;
    for (int i = children.size() - 1; i >= 0; i--)
    {
        size_t num_data_at_child = n_data_desc[i].first;
        size_t num_data_at_desc_of_child = n_data_desc[i].second;
        size_t n_data_i = num_data_at_child + num_data_at_desc_of_child;
        //double temp = bounded_beta(random, n_data_i + 1, params.get_gamma() + (total_num_data_at_desc - n_data_i - cumulative_num_data));
        double temp = bounded_beta(random, n_data_i + 1, params.get_gamma() + cumulative_num_data);
        node->set_psi_stick(i, temp);
        cumulative_num_data += n_data_i;
    }
    assert(total_num_data_at_desc == cumulative_num_data);

    return make_pair(n_data, total_num_data_at_desc);
}

// public functions

// getters and setters

// get the total number of nodes in this tree (state)
template <class S, class C, class P>
size_t TSSBState<S,C,P>::get_num_nodes()
{
    vector<Node<S,P> *> nodes;
    get_all_nodes(nodes);
    return nodes.size();
}

// get the node where this SomaticSNV (ssm) first appears in
template <class S, class C, class P>
Node<S,P> *TSSBState<S,C,P>::get_node(const S *datum)
{
    return datum2node[datum];
}

template <class S, class C, class P>
void TSSBState<S,C,P>::get_all_nodes(bool non_empty, vector<Node<S,P> *> &ret) const
{
    get_all_nodes(non_empty, root, ret);
}

template <class S, class C, class P>
void TSSBState<S,C,P>::get_all_nodes(bool non_empty, Node<S,P> *root_node, vector<Node<S,P> *> &ret)
{
    // get nodes in the subtree rooted at root
    if (ret.size() > 0) {
        cout << "Warning: get_all_nodes() is going to clear vector ret." << endl;
        ret.clear();
    }
    if (root_node == 0)
        return;
    
    set<Node<S,P> *> nodes;
    Node<S,P> *node = 0;
    
    queue<Node<S,P> *> q;
    q.push(root_node);
    while (!q.empty())
    {
        node = q.front();
        q.pop();
        if (!nodes.count(node)) {
            if (!non_empty) {
                ret.push_back(node);
                nodes.insert(node);
            } else {
                if (node->get_num_data() > 0) {
                    ret.push_back(node);
                    nodes.insert(node);
                }
            }
        }
        unordered_map<size_t, pair<double, Node<S,P> *> > &children = node->get_idx2child();
        for (size_t i = 0; i < children.size(); i++)
        {
            Node<S,P> *child = children[i].second;
            q.push(child);
        }
    }
}

// get the nodes along with their mass
// the nodes are ordered by depth first traversal
template <class S, class C, class P>
void TSSBState<S,C,P>::get_mixture(Node<S,P> *node, double mass, vector<pair<double, Node<S,P> *> > &ret)
{
    double node_mass = mass * node->get_nu_stick();
    ret.push_back(make_pair(node_mass, node));
    unordered_map<size_t, pair<double, Node<S,P> *> > &children = node->get_idx2child();
    double cum_prod = 1.0;
    vector<double> weights;
    for (size_t i = 0; i < children.size(); i++)
    {
        Node<S,P> *child = children[i].second;
        double val = cum_prod * (1 - children[i].first);
        weights.push_back(1 - val);
        cum_prod *= (1 - children[i].first);
        get_mixture(child, (1 - node->get_nu_stick()) * mass * weights[i], ret);
    }
}

template <class S, class C, class P>
void TSSBState<S,C,P>::insert_datum(const gsl_rng *random, S* datum, const ModelParams &params)
{
    double u = gsl_ran_flat(random, 0, 1);
    Node<S,P> *node = Node<S,P>::find_node(random, u, root, params);
    bulk_data->push_back(datum);
    assign_data_point(0, node, datum, params);
    
    // do not update sc_cache as it will initialized and called in compute_likelihood_sc_cached
}

template <class S, class C, class P>
void TSSBState<S,C,P>::initialize_sc_cache(const ModelParams &model_params)
{
    sc_cache.resize(sc_data->size());
    vector<Node<S,P> *> all_nodes;
    get_all_nodes(false, all_nodes);
    unordered_map<Node<S,P> *, unordered_set<const S *> > node2snvs;
    for (Node<S,P> *v : all_nodes)
    {
        unordered_set<const S *> snvs;
        Node<S,P>::get_dataset(v, snvs);
        node2snvs[v] = snvs;
    }
    for (size_t c = 0; c < sc_data->size(); c++) {
        for (auto v : all_nodes) {
            initialize_sc_cache(c, v, node2snvs[v], model_params);
        }
    }
}

template <class S, class C, class P>
void TSSBState<S,C,P>::set_sc_data(vector<C *> *sc_data, const ModelParams &model_params)
{
    this->sc_data = sc_data;
    initialize_sc_cache(model_params);
}

// get all nodes in this tree with at least one datum
// via breadth-first traversal
// populate ret
template <class S, class C, class P>
void TSSBState<S,C,P>::get_all_nodes(vector<Node<S,P> *> &ret) const
{
    get_all_nodes(false, root, ret);
}

template <class S, class C, class P>
gsl_matrix *TSSBState<S,C,P>::get_ancestral_matrix(TSSBState<S,C,P> &state)
{
    size_t N = state.bulk_data->size();
    vector<unordered_set<S *> > ancestors(N);
    for (size_t i = 0; i < N; i++) {
        const S *datum = state.bulk_data->at(i);
        Node<S,P> *v = state.get_node(datum);
        // trace up to the root and set the row of A
        while (v != state.root) {
            v = v->get_parent_node();
            unordered_set<S *> data = v->get_data();
            ancestors[i].insert(data.begin(), data.end());
        }
    }

    // A_ij = 1 if i occurs before j
    gsl_matrix *A = gsl_matrix_alloc(N, N);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            if (ancestors[j].count(state.bulk_data->at(i))) {
                gsl_matrix_set(A, i, j, 1);
            } else {
                gsl_matrix_set(A, i, j, 0);
            }
        }
    }
    
    return A;
}

/**********************
 Functions for sampling
 **********************/

// compute (or get) the log likelihood of the current state
template <class S, class C, class P>
double TSSBState<S,C,P>::compute_log_likelihood_bulk(const ModelParams &model_params)
{
    double log_lik_bulk = 0.0;
    // compute the log likelihood over data S
    for (size_t i = 0; i < bulk_data->size(); i++) {
        Node<S,P> *node = datum2node[(*bulk_data)[i]];
        //double log_val = node->compute_log_likelihood_of_datum(*(*bulk_data)[i], model_params);
        double log_val = log_lik_datum(node, (*bulk_data)[i], model_params);
        log_lik_bulk += log_val;
    }
    return log_lik_bulk;
}

// update stick lengths, and exchange stick ordering
template <class S, class C, class P>
void TSSBState<S,C,P>::update_sticks(const gsl_rng *random, const ModelParams &params)
{
    descend_and_sample_sticks(random, root, params);
    reorder_sticks(random, params);
    descend_and_sample_sticks(random, root, params);
}

template <class S, class P>
double log_prod_beta(vector<Node<S,P> *> &nodes, const ModelParams &params)
{
    double log_sum = 0.0;
    double alpha;
    for (Node<S,P> *node : nodes) {
        if (node->get_depth() == 0 && node->get_nu_stick() == 0) // root is assigned 0 nu-stick as it denotes healthy population that should not have any mutation assigned to it
            continue;
        double x = node->get_nu_stick();
        alpha = params.alpha(node->get_name_vec());
        double y = log_beta_pdf(x, 1.0, alpha);
        log_sum += y;
        // below: for a small value of alpha, beta(1, alpha) is essentially a point mass at 1
        // otherwise, the value of y may get too large and dominate
        //log_sum += y > 20 ? 20 : y;
    }
    return log_sum;
}

template <class S, class C, class P>
void TSSBState<S,C,P>::sample_alpha0(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<Node<S,P> *> &nodes)
{
    //double curr_log_lik = log_prod_beta(nodes, params);
    double curr_log_lik_slice = log_prod_beta(nodes, params) + log(gsl_ran_flat(random, 0, 1));
    double new_log_lik;
    double alpha_0_ll = params.get_alpha0_bound(false);
    double alpha_0_uu = params.get_alpha0_bound(true);
    double curr_alpha0 = params.get_alpha0();
    double new_alpha0;
    double u;
    //for (size_t i = 0; i < n_mh_iter; i++)
    while (true)
    {
        // propose a value
        //new_alpha0 = params.get_alpha0() + gsl_ran_gaussian(random, params.get_alpha0_sigma());
        u = gsl_ran_flat(random, 0, 1);
        new_alpha0 = (alpha_0_uu - alpha_0_ll) * u + alpha_0_ll;

//        if (new_alpha0 <= alpha_0_ll || new_alpha0 > alpha_0_uu) {
//            continue;
//        }
        params.set_alpha0(new_alpha0);
        new_log_lik = log_prod_beta(nodes, params);
        //double log_unif = log(gsl_ran_flat(random, 0.0, 1.0));
        //if (log_unif < (new_log_lik - curr_log_lik)) {
        if (new_log_lik > curr_log_lik_slice) {
            // accept
            //curr_log_lik = new_log_lik;
            //curr_alpha0 = new_alpha0;
            break;
        } else {
            // revert
            //params.set_alpha0(curr_alpha0);
            if (new_alpha0 < curr_alpha0) {
                alpha_0_ll = new_alpha0;
            } else if (new_alpha0 > curr_alpha0) {
                alpha_0_uu = new_alpha0;
            }
        }
        if (abs(alpha_0_uu - alpha_0_ll) < 1e-6) {
            // revert
            params.set_alpha0(curr_alpha0);
            cout << "alpha0 slice sampler shrank to 0!" << endl;
            break;
        }
    }
}

template <class S, class C, class P>
void TSSBState<S,C,P>::sample_lambda(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<Node<S,P> *> &nodes)
{
    //double curr_log_lik = log_prod_beta(nodes, params);
    double curr_log_slice = log_prod_beta(nodes, params) + log(gsl_ran_flat(random, 0, 1));
    double new_log_lik;
    double ll = params.get_lambda_bound(false);
    double uu = params.get_lambda_bound(true);
    double curr_val = params.get_lambda();
    double new_val;
    double u;
    //for (size_t i = 0; i < n_mh_iter; i++)
    while (true)
    {
        // propose a value
        //new_val = params.get_lambda() + gsl_ran_gaussian(random, params.get_lambda_sigma());
        u = gsl_ran_flat(random, 0, 1);
        new_val = (uu - ll) * u + ll;
        
//        if (new_val <= ll || new_val > uu) {
//            continue;
//        }
        params.set_lambda(new_val);
        new_log_lik = log_prod_beta(nodes, params);
        //double log_unif = log(gsl_ran_flat(random, 0.0, 1.0));
        //if (log_unif < (new_log_lik - curr_log_lik)) {
        if (new_log_lik > curr_log_slice) {
            // accept
            //curr_log_lik = new_log_lik;
            //curr_val = new_val;
            break;
        } else {
            // revert
            //params.set_lambda(curr_val);
            if (new_val < curr_val) {
                ll = new_val;
            } else if (new_val > curr_val) {
                uu = new_val;
            }
        }
        if (abs(uu - ll) < 1e-6) {
            // revert
            params.set_lambda(curr_val);
            cout << "lambda slice sampler shrank to 0!" << endl;
            break;
        }
    }
}

template <class S, class C, class P>
void TSSBState<S,C,P>::sample_gamma(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<double> psi_sticks)
{
    //double curr_log_lik = log_prod_beta(psi_sticks, params.get_gamma());
    double curr_log_slice = log_prod_beta(psi_sticks, params.get_gamma()) + log(gsl_ran_flat(random, 0, 1));
    double new_log_lik;
    double ll = params.get_gamma_bound(false);
    double uu = params.get_gamma_bound(true);
    double curr_gamma = params.get_gamma();
    double new_gamma;
    //for (size_t i = 0; i < n_mh_iter; i++)
    while (true)
    {
        // propose a value
        new_gamma = (uu - ll) * gsl_ran_flat(random, 0, 1) + ll;
        //new_gamma = gsl_ran_gaussian(random, params.get_gamma_sigma()) + curr_gamma;
//        if (new_gamma < ll || new_gamma > uu)
//            continue;
        params.set_gamma(new_gamma);
        new_log_lik = log_prod_beta(psi_sticks, params.get_gamma());
        //double log_unif = log(gsl_ran_flat(random, 0.0, 1.0));
        //if (log_unif < (new_log_lik - curr_log_lik)) {
        if (new_log_lik > curr_log_slice) {
            // accept
            //cout << "new log_lik: " << new_log_lik << ", " << curr_log_lik << endl;
            //curr_gamma = new_gamma;
            //curr_log_lik = new_log_lik;
            break;
        } else {
            //params.set_gamma(curr_gamma);
            if (new_gamma < curr_gamma) {
                ll = new_gamma;
            } else if (new_gamma > curr_gamma) {
                uu = new_gamma;
            }
        }
        if (abs(uu - ll) < 1e-6) {
            // revert
            params.set_gamma(curr_gamma);
            cout << "gamma slice sampler shrank to 0!" << endl;
            break;
        }
    }
}

template <class S, class C, class P>
void TSSBState<S,C,P>::update_hyper_params(const gsl_rng *random,
                         size_t n_mh_iter,
                         TSSBState<S,C,P> &state,
                         ModelParams &params)
{
    vector<Node<S,P> *> all_nodes;
    state.get_all_nodes(all_nodes);
    vector<double> psi_sticks;
    for (Node<S,P> *node : all_nodes) {
        for (size_t i = 0; i < node->get_idx2child().size(); i++) {
            psi_sticks.push_back(node->get_idx2child()[i].first);
        }
    }
    
    sample_alpha0(random, n_mh_iter, params, all_nodes);
    sample_lambda(random, n_mh_iter, params, all_nodes);
    sample_gamma(random, n_mh_iter, params, psi_sticks);
}

template <class S, class C, class P>
string TSSBState<S,C,P>::print()
{
    string ret = "";
    queue<Node<S,P> *> q;
    q.push(root);
    while (!q.empty()) {
        Node<S,P> *node = q.front();
        q.pop();

        ret += node->print() + "\n";
        // enqueue children nodes to the queue
        unordered_map<size_t, pair<double, Node<S,P> *> > idx2child = node->get_idx2child();
        for (size_t i = 0; i < idx2child.size(); i++) {
            q.push(idx2child[i].second);
        }
    }

    ret += "\n";
    // print data assignments
    //double log_lik = 0.0, log_lik_datum = 0.0;
//    for (size_t idx = 0; idx < bulk_data->size(); idx++) {
//        Node<S,P> *node = datum2node[(*bulk_data)[idx]];
//        ret += "Data point " + to_string(idx) + " assigned to <" + (node->is_root() ? "root" : node->get_name()) + ">\n";
        //S *datum = (*data)[idx];
//        log_lik_datum = node->compute_log_likelihood_of_datum(*datum, params);
//        ret += "Log likelihood: " + to_string(log_lik_datum) + "\n";
//        log_lik += log_lik_datum;
//    }
//    ret += "Complete data log lik: " + to_string(log_lik) + "\n";
    return ret;
}

template <class S, class C, class P>
TSSBState<S,C,P> *TSSBState<S,C,P>::construct_trivial_state(Node<S,P> *root,
                                                            double (*log_lik_datum)(const Node<S,P> *v, const S *s, const ModelParams &params),
                                                            double (*log_lik_sc_at_site)(const S *s, const C *c, bool has_snv, const ModelParams &params))
{
    auto state = new TSSBState<S,C,P>(root);
    state->log_lik_datum = log_lik_datum;
    state->log_lik_sc_at_site = log_lik_sc_at_site;
    return state;
}

template <class S, class C, class P>
void TSSBState<S,C,P>::clear_cache()
{
    if (sc_cache.size() == 0)
        return;
    
    vector<Node<S,P> *> active_nodes_vec;
    get_all_nodes(false, root, active_nodes_vec);
    unordered_set<Node<S,P> *> active_nodes(active_nodes_vec.begin(), active_nodes_vec.end());
    
    unordered_set<Node<S,P> *> purge_set;
    for (auto it = sc_cache[0].begin(); it != sc_cache[0].end(); ++it)
    {
        Node<S,P> *v = it->first;
        if (active_nodes.count(v) == 0) {
            purge_set.insert(v);
        }
    }
    for (size_t c = 0; c < sc_data->size(); c++) {
        unordered_map<Node<S,P> *, double> &cache_c = sc_cache[c];
        for (auto v : purge_set)
            cache_c.erase(v);
    }
}

#endif /* tssb_state_hpp */
