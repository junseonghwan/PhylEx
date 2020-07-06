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
#include "model_params.hpp"
#include "numerical_utils.hpp"

using namespace std;

class TSSBState
{
    CloneTreeNode *root;

    double log_lik = DOUBLE_NEG_INF;
    double log_lik_bulk = DOUBLE_NEG_INF;
    double log_lik_sc = DOUBLE_NEG_INF;

    vector<BulkDatum *> *bulk_data;
    vector<SingleCellData *> *sc_data;

    unordered_map<const BulkDatum *, CloneTreeNode *> datum2node;

    // vector index for sc_cache is the same as the sc_data idx
    // for eaach cell, store the mapping from node to double (log likelihood of being assigned to that node)
    vector<unordered_map<CloneTreeNode *, double> > sc_cache;

    // function pointer to compute likelihood of assigning bulk datum to node
    double (*log_lik_datum)(size_t region,
                            const CloneTreeNode *v,
                            const BulkDatum *s,
                            const ModelParams &params) = 0;
    // function pointer to compute likelihood of single cell at site s if it has the variant given by has_snv
    double (*log_lik_sc_at_site)(const BulkDatum *s, const SingleCellData *c, bool has_snv, const ModelParams &params) = 0;

    // private constructor -- used by simulator, where data has not been generated yet
    TSSBState(CloneTreeNode *root);

    // helper functions for initializing assignment of SNV
    //void initialize_data_assignment(const gsl_rng *random, BulkDatum *datum, const ModelParams &params);
    void initialize_data_assignment(const gsl_rng *random, size_t mut_id, const ModelParams &params);
    // helper functions for performing slice sampling on an SNV
    //void slice_sample_data_assignment(const gsl_rng *random, BulkDatum *datum, const ModelParams &params);
    void slice_sample_data_assignment(const gsl_rng *random, size_t mut_id, const ModelParams &model_params);
    // assign data to a node
    //void assign_data_point(CloneTreeNode *curr_node, CloneTreeNode *new_node, BulkDatum *datum, const ModelParams &model_params);
    void assign_data_point(CloneTreeNode *curr_node, CloneTreeNode *new_node, size_t mut_id, const ModelParams &model_params);

    void initialize_sc_cache(const ModelParams &model_params);
    void initialize_sc_cache(size_t c, CloneTreeNode *v, unordered_set<const BulkDatum *> &snvs, const ModelParams &model_params);
    //void update_sc_cache(CloneTreeNode *curr_node, CloneTreeNode *new_node, BulkDatum *datum, const ModelParams &params);
    void update_sc_cache(CloneTreeNode *curr_node, CloneTreeNode *new_node, size_t mut_id, const ModelParams &params);
    void print_cache();

    //double compute_loglik_sc(unordered_set<const BulkDatum *> &snvs, SingleCellData *cell, const ModelParams &params);
    double compute_loglik_sc(unordered_set<const BulkDatum *> &snvs, size_t cell_id);
    
    EigenMatrix sc_presence_matrix_, sc_absence_matrix_;
    void PreComputeScLikelihood(const ModelParams &model_params);

    double LogLikDatum(CloneTreeNode *node,
                       BulkDatum *datum,
                       const ModelParams &model_params);
public:
    TSSBState(const gsl_rng *random,
              CloneTreeNode *root,
              const ModelParams &params,
              double (*log_lik_datum)(size_t region,
                                      const CloneTreeNode *node,
                                      const BulkDatum *datum,
                                      const ModelParams &params),
              double (*log_lik_sc_at_site)(const BulkDatum *s,
                                           const SingleCellData *c,
                                           bool has_snv,
                                           const ModelParams &params),
              vector<BulkDatum *> *bulk_data,
              vector<SingleCellData *> *sc_data);

    //void insert_datum(const gsl_rng *random, BulkDatum * datum, const ModelParams &params);
    void set_sc_data(vector<SingleCellData *> *sc_data, const ModelParams &model_params);
    
const vector<BulkDatum *> &get_data() const;

    // getters
    size_t get_num_nodes();
    inline CloneTreeNode *get_root() { return root; }
    CloneTreeNode *get_node(const BulkDatum *datum);
    static void get_mixture(CloneTreeNode *node, double mass, vector<pair<double, CloneTreeNode *> > &ret);
    static void get_all_nodes(bool non_empty, CloneTreeNode *root_node, vector<CloneTreeNode *> &ret);
    void get_all_nodes(bool non_empty, vector<CloneTreeNode *> &ret) const;
    void get_all_nodes(vector<CloneTreeNode *> &ret) const;
    inline bool is_root(const CloneTreeNode *node) const { return root == node; }

    // operations relating to sticks
    void update_sticks(const gsl_rng *random, const ModelParams &params);
    pair<size_t, size_t> descend_and_sample_sticks(const gsl_rng *random, CloneTreeNode *node, const ModelParams &params);
    void reorder_sticks(const gsl_rng *random, const ModelParams &model_params);

    void resample_data_assignment(const gsl_rng *random, const ModelParams &params);
    //void move_datum(CloneTreeNode *node, BulkDatum *datum, const ModelParams &model_params);
    void move_datum(CloneTreeNode *new_node, size_t mut_id, const ModelParams &model_params);
    
    double compute_log_likelihood_bulk(size_t region,
                                       const ModelParams &model_params);
    double compute_log_likelihood_bulk(const ModelParams &params);
    double compute_log_likelihood_sc_cached(const ModelParams &params, bool verbose=false);

    double get_log_lik();
    static double get_log_prior_assignment(CloneTreeNode *root);
    inline double get_log_lik_bulk() { return log_lik_bulk; }
    inline double get_log_lik_sc() { return log_lik_sc; }
    //double compute_log_likelihood_sc_subset(const gsl_rng *random, const ModelParams &params);

    // for output and debugging
    string print();
    
    // static functions
    static gsl_matrix *get_ancestral_matrix(TSSBState &state);
    static TSSBState *construct_trivial_state(CloneTreeNode *root,
                            double (*log_lik_datum)(size_t region,
                                                    const CloneTreeNode *v,
                                                    const BulkDatum *s,
                                                    const ModelParams &params),
                            double (*log_lik_sc_at_site)(const BulkDatum *s,
                                                         const SingleCellData *c,
                                                         bool has_snv,
                                                         const ModelParams &params));
    static double compute_loglik_sc(const vector<BulkDatum *> &bulk_data,
                                    unordered_set<const BulkDatum *> &snvs,
                                    SingleCellData *cell,
                                    const ModelParams &params,
                                    double (*log_lik_sc_at_site)(const BulkDatum *s, const SingleCellData *c, bool has_snv, const ModelParams &params));
    static double compute_log_likelihood_sc(CloneTreeNode *root,
                                            const vector<BulkDatum *> &bulk_data,
                                            const vector<SingleCellData *> &sc_data,
                                            double (*log_lik_sc_at_site)(const BulkDatum *s, const SingleCellData *c, bool has_snv, const ModelParams &params),
                                            const ModelParams &params,
                                            bool verbose=false);
    static void sample_alpha0(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<CloneTreeNode *> &nodes);
    static void sample_lambda(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<CloneTreeNode *> &nodes);
    static void sample_gamma(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<double> psi_sticks);
    
    static void update_hyper_params(const gsl_rng *random,
                                     size_t n_mh_iter,
                                     TSSBState &state,
                                     ModelParams &params);

    friend size_t descend_and_cull(CloneTreeNode *node);
    friend void descend_and_update_names(CloneTreeNode *node);
    void clear_cache();
};

void update_params(size_t region,
                   vector<CloneTreeNode *> &nodes,
                   double *new_clone_freq);
double update_cellular_prev_recursive(size_t region, CloneTreeNode *node);
void get_clone_freqs(size_t region,
                     TSSBState &state,
                     double *clone_freqs);
void sample_params_bottom_up(const gsl_rng *random,
                             size_t n_mh_iter,
                             TSSBState &tree,
                             const ModelParams &params);
double sample_params_dirichlet(const gsl_rng *random,
                               size_t n_mh_iter,
                               TSSBState &tree,
                               const ModelParams &params);
void cull(CloneTreeNode *root);
bool check_clone_freq(size_t region, CloneTreeNode *root);

#endif /* tssb_state_hpp */
