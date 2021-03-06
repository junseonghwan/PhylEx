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

#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include "clone_node.hpp"
#include "model_params.hpp"
#include "numerical_utils.hpp"

using namespace std;

class TSSBState
{
    // We will have two roots:
    // 1. A root that represents healthy cells and
    // 2. A child of the root that represents the root of the cancer clones.
    CloneTreeNode *root;
    CloneTreeNode *ancestral_clone;

    //double log_lik = DOUBLE_NEG_INF;
    //double log_lik_bulk = DOUBLE_NEG_INF;
    double log_lik_sc = DOUBLE_NEG_INF;
    bool use_geometric_mean_;

    vector<BulkDatum *> *bulk_data_;
    vector<bool> has_sc_coverage_;
    vector<SingleCellData *> *sc_data;

    unordered_map<const BulkDatum *, CloneTreeNode *> datum2node;

    // function pointer to compute likelihood of assigning bulk datum to node
    double (*log_lik_datum)(size_t region,
                            const CloneTreeNode *v,
                            const BulkDatum *s,
                            const ModelParams &params) = 0;
    // function pointer to compute likelihood of single cell at site s if it has the variant given by has_snv
    double (*log_lik_sc_at_site)(size_t loci_idx,
                                 const BulkDatum *s,
                                 const SingleCellData *c,
                                 bool has_snv,
                                 const ModelParams &params) = 0;

    // Helper functions for initializing assignment of SNV.
    // It assign all SNVs to the first child of root.
    void InitializeDataAssignment(const gsl_rng *random, size_t mut_id, const ModelParams &params);
    // helper functions for performing slice sampling on an SNV
    void slice_sample_data_assignment(const gsl_rng *random,
                                      size_t mut_id,
                                      const ModelParams &model_params);
    void slice_sample_data_assignment_with_sc(const gsl_rng *random,
                                              size_t mut_id,
                                              const ModelParams &model_params);
    void AssignDatum(CloneTreeNode *curr_node, CloneTreeNode *new_node, size_t mut_id, const ModelParams &model_params, bool update_cache = true);

    void InitializeCacheForNode(CloneTreeNode *v);
    void UpdateSingleCellCache(CloneTreeNode *curr_node, CloneTreeNode *new_node, size_t mut_id, const ModelParams &params);

    double compute_loglik_sc(CloneTreeNode *v, size_t cell_id);

    // N x C matrix storing likelihood of single cell reads for mutation n, cell c
    // when the cell carries the mutation (presnce) and when it doesn't (absence).
    DoubleMatrix sc_presence_matrix_, sc_absence_matrix_;
    // Pre-compute single cell likelihoods and determine SNVs that don't have
    // any single cell coverage.
    void ProcessSingleCellData(const ModelParams &model_params);

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
              double (*log_lik_sc_at_site)(size_t loci_idx,
                                           const BulkDatum *s,
                                           const SingleCellData *c,
                                           bool has_snv,
                                           const ModelParams &params),
              vector<BulkDatum *> *bulk_data,
              vector<SingleCellData *> *sc_data,
              bool use_geometric_mean = false);

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

    void resample_data_assignment(const gsl_rng *random,
                                  const ModelParams &params,
                                  bool use_sc = true);
    //void move_datum(CloneTreeNode *node, BulkDatum *datum, const ModelParams &model_params);
    void move_datum(CloneTreeNode *new_node, size_t mut_id, const ModelParams &model_params);
    
    double compute_log_likelihood_bulk(size_t region,
                                       const ModelParams &model_params);
    double compute_log_likelihood_bulk(const ModelParams &params);
    double compute_log_likelihood_sc_cached(bool verbose=false);

    //double get_log_lik();
    static double get_log_prior_assignment(CloneTreeNode *root);
    //inline double get_log_lik_bulk() { return log_lik_bulk; }
    inline double get_log_lik_sc() { return log_lik_sc; }

    // for output and debugging
    string print();
    
    // static functions
    static gsl_matrix *get_ancestral_matrix(TSSBState &state);
    
    double compute_log_likelihood_sc(bool verbose=false);
    
    static void sample_alpha0(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<CloneTreeNode *> &nodes);
    static void sample_lambda(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<CloneTreeNode *> &nodes);
    static void sample_gamma(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<double> psi_sticks);
    
    static void update_hyper_params(const gsl_rng *random,
                                     size_t n_mh_iter,
                                     TSSBState &state,
                                     ModelParams &params);

    friend size_t descend_and_cull(CloneTreeNode *node);
    friend void descend_and_update_names(CloneTreeNode *node);
    //void clear_cache();
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
double ComputeSingleCellLikelihood(CloneTreeNode *root,
                    vector<BulkDatum *> &bulk_data,
                    vector<SingleCellData *> &sc_data,
                    const ModelParams &model_params);

#endif /* tssb_state_hpp */
