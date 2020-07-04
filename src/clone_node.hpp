//
//  clone_node.hpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-06-27.
//

#ifndef clone_node_hpp
#define clone_node_hpp

#include <string>

#include "bulk_datum.hpp"
#include "node.hpp"
#include "single_cell.hpp"

using namespace std;

class CloneTreeNodeParam
{
    double clone_freq = 0.0;
    double cellular_prev = 0.0;
public:
    inline double get_clone_freq() const { return clone_freq; }
    inline double get_cellular_prev() const { return cellular_prev; }
    void set_clone_freq(double val);
    void set_cellular_prev(double val);
    bool is_consistent();
};

class CloneTreeNode : public Node<BulkDatum,CloneTreeNodeParam>
{
public:
    CloneTreeNode(size_t child_idx, Node<BulkDatum,CloneTreeNodeParam> *parent);

    const CloneTreeNodeParam &get_node_parameter() const override;
    void sample_node_parameters(const gsl_rng *random,
                                const ModelParams &params,
                                Node<BulkDatum,CloneTreeNodeParam> *parent) override;
    Node<BulkDatum,CloneTreeNodeParam> *spawn_child(double psi) override;
    string print() override;

    void change_clone_freq(double diff);
    void set_clone_freq(double new_val);
    void set_cellular_prev(double new_val);
    //double compute_expectation_of_xi(const ModelParams &model_params);

    static CloneTreeNode *create_root_node();
    static void get_snvs(Node<BulkDatum,CloneTreeNodeParam> *node, unordered_set<Locus> &ret);
};

double ScLikelihood(const BulkDatum *s,
                     const SingleCellData *sc,
                     bool has_snv,
                     const ModelParams &model_params);
double BulkLogLikWithTotalCopyNumber(const Node<BulkDatum, CloneTreeNodeParam> *node,
                       const BulkDatum *s,
                       const ModelParams &model_params);
double zero_bulk_likelihood(const Node<BulkDatum, CloneTreeNodeParam> *node,
                            const BulkDatum *datum,
                            const ModelParams &model_params);
double BulkLogLikWithCopyNumberProfile(const Node<BulkDatum, CloneTreeNodeParam> *node,
                                       const BulkDatum *datum,
                                       const ModelParams &model_params);
double BulkLogLikWithGenotype(const Node<BulkDatum, CloneTreeNodeParam> *node,
                              const BulkDatum *datum,
                              const ModelParams &model_params);
void update_params(vector<Node<BulkDatum,CloneTreeNodeParam> *> &nodes, double *new_clone_freq);
double update_cellular_prev_recursive(Node<BulkDatum,CloneTreeNodeParam> *node);
void get_clone_freqs(TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &state,
                     double *clone_freqs);
void sample_params_bottom_up(const gsl_rng *random,
                             size_t n_mh_iter,
                             TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &tree,
                             const ModelParams &params);
double sample_params_dirichlet(const gsl_rng *random,
                             size_t n_mh_iter,
                             TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &tree,
                             const ModelParams &params);
void cull(Node<BulkDatum,CloneTreeNodeParam> *root);

#endif /* clone_node_hpp */
