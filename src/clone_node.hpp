//
//  clone_node.hpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-06-27.
//

#ifndef clone_node_hpp
#define clone_node_hpp

#include <iostream>
#include <limits>
#include <map>
#include <memory.h>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <gsl/gsl_randist.h>

#include "bulk_datum.hpp"
#include "model_params.hpp"
#include "numerical_utils.hpp"
#include "sampling_utils.hpp"
#include "single_cell.hpp"

using namespace std;

class CloneTreeNodeParam
{
    vector<double> clone_freqs;
    vector<double> cellular_prevs;
public:
    CloneTreeNodeParam(size_t region_count);
    inline const vector<double> &get_clone_freqs() const { return clone_freqs; }
    inline const vector<double> &get_cellular_prevs() const { return cellular_prevs; }
    inline double get_clone_freqs(size_t region) const { return clone_freqs[region]; }
    inline double get_cellular_prevs(size_t region) const { return cellular_prevs[region]; }
    void set_clone_freq(size_t idx, double val);
    void set_cellular_prev(size_t idx, double val);
    void set_clone_freq(vector<double> &vec);
    void set_cellular_prev(vector<double> &vec);
    bool is_consistent();
    void SetRootParameters();
    size_t GetRegionCount() const {
        return clone_freqs.size();
    }
    string GetCloneFreqsAsString();
    string GetCellularPrevsAsString();
};

class CloneTreeNode
{
    CloneTreeNodeParam param;
    
    vector<size_t> name;
    double nu = 0.0;
    CloneTreeNode *parent_node = 0;
    unordered_set<const BulkDatum *> data;
    vector<double> sc_cache;

    // map: child branch idx -> pair<psi_stick, Node *>
    unordered_map<size_t, pair<double, CloneTreeNode *> > idx2child;
    bool operator==(const CloneTreeNode &other) const
    {
        if (this->parent_node != other.parent_node)
            return false;
        if (name.size() != other.name.size())
            return false;
        size_t n = name.size();
        return (name[n-1] == other.name[n-1]);
    }

    // change the last part of the node's name to j -- used extensively by reorder_sticks
    void edit_name(size_t j);
    
    // For creating the root node.
    CloneTreeNode(size_t region_count);
    // For creating non-root node.
    // Assumption: parent != 0.
    CloneTreeNode(size_t child_idx, CloneTreeNode *parent);
public:
    ~CloneTreeNode();

    bool IsCacheAllocated(size_t cell_count);
    void AllocateCache(size_t cell_count);
    void UpdateCache(size_t c, double val);
    double GetScCache(size_t c);
    string get_name() const;
    inline const vector<size_t> &get_name_vec() const { return name; };
    inline size_t get_depth() const { return name.size() - 1; }
    size_t get_num_data() const;
    double get_nu_stick() const;
    CloneTreeNode *get_parent_node() const;
    const unordered_set<const BulkDatum *> &get_data() const;
    const pair<double, CloneTreeNode *> &get_child(size_t child_idx) const;
    size_t get_num_children() const;
    bool is_root() const;
    bool is_leaf() const;
    unordered_map<size_t, pair<double, CloneTreeNode *> > &get_idx2child();

    // setters
    void set_nu_stick(double nu);
    void set_psi_stick(size_t child_idx, double psi);
    
    // add/remove datum
    void add_datum(BulkDatum *datum);
    void remove_datum(BulkDatum *datum);
    bool contains_datum(BulkDatum *datum) const;
    
    // operations on sticks
    void cull(unordered_set<size_t> &cull_list);
    void reorder_sticks(const gsl_rng *random, const ModelParams &params);
    void reset_children_names();
    
    // Returns true if this node is ancestor of other or other == this.
    bool IsAncestorOf(CloneTreeNode *other);
    bool IsDescendantOf(CloneTreeNode *other);
    
    // identify the branch that contains u, return the corresponding child Node
    CloneTreeNode *locate_child(const gsl_rng *random, double &u, const ModelParams &hyper_params);
    void InitializeChild(const gsl_rng *random,
                         const ModelParams &params);
    
    static string form_node_string(string curr_node_str, size_t branch);
    static string get_parent_string(string curr_node_str);
    
    static void breadth_first_traversal(CloneTreeNode *root, vector<CloneTreeNode *> &ret, bool non_empty = false);
    static CloneTreeNode *find_node(const gsl_rng *random, double u, CloneTreeNode *root, const ModelParams &params);
    // return all data from ancestors including node itself
    static void GetDataset(CloneTreeNode *node, unordered_set<const BulkDatum *> &dataset);
    static void get_cluster_labels(CloneTreeNode *root,
                                   const vector<BulkDatum *> &data,
                                   vector<unsigned int> &cluster_labels);
    static void construct_datum2node(vector<CloneTreeNode *> &all_nodes,
                                     unordered_map<const BulkDatum *, CloneTreeNode *> &datum2node);
    static gsl_matrix *GetAncestralMatrix(CloneTreeNode *root,
                                          const vector<BulkDatum *> &bulk_data,
                                          const unordered_map<const BulkDatum *, CloneTreeNode *> &datum2node);

    static bool less(CloneTreeNode *l, CloneTreeNode *r);
    friend bool operator<(const CloneTreeNode& lhs, const CloneTreeNode& rhs) {
        vector<string> lhs_arr, rhs_arr;
        boost::split(lhs_arr, lhs.name, boost::is_any_of("_"));
        boost::split(rhs_arr, rhs.name, boost::is_any_of("_"));
        if (lhs_arr.size() < rhs_arr.size())
            return true;
        else if (lhs_arr.size() > rhs_arr.size())
            return false;
        else {
            // same length, compare one by one
            for (size_t i = 0; i < lhs_arr.size(); i++) {
                if (stoi(lhs_arr[i]) < stoi(rhs_arr[i])) {
                    return true;
                } else if (stoi(lhs_arr[i]) > stoi(rhs_arr[i])) {
                    return false;
                }
            }
            return false;
        }
        //return tie(lhs.name) < tie(rhs.name);
    }
    friend bool operator> (const CloneTreeNode& lhs, const CloneTreeNode& rhs){ return rhs < lhs; }
    friend bool operator<=(const CloneTreeNode& lhs, const CloneTreeNode& rhs){ return !(lhs > rhs); }
    friend bool operator>=(const CloneTreeNode& lhs, const CloneTreeNode& rhs){ return !(lhs < rhs); }

//    inline EigenVectorRef GetCellularPrevs() const {
//        return param.get_cellular_prevs();
//    }
//    inline EigenVectorRef GetCloneFreqs() const {
//        return param.get_clone_freqs();
//    }
    inline double GetCellularPrevs(size_t region) const {
        return param.get_cellular_prevs(region);
    }
    inline double GetCloneFreqs(size_t region) const {
        return param.get_clone_freqs(region);
    }
    CloneTreeNodeParam &get_node_parameter();
    void sample_node_parameters(const gsl_rng *random,
                                const ModelParams &params,
                                CloneTreeNode *parent);
    CloneTreeNode *spawn_child(double psi);
    string print();

    void set_clone_freq(size_t idx, double new_val);
    void set_cellular_prev(size_t idx, double new_val);
    void set_clone_freq(vector<double> &vec);
    void set_cellular_prev(vector<double> &vec);

    static CloneTreeNode *create_root_node(size_t region_count);
    static void RetrieveLoci(CloneTreeNode *node, unordered_set<Locus> &ret);
};

double ScLikelihood(size_t loci_idx,
                    const BulkDatum *s,
                    const SingleCellData *sc,
                    bool has_snv,
                    const ModelParams &model_params);
double BulkLogLikWithTotalCopyNumber(size_t region,
                                     const CloneTreeNode *node,
                                     const BulkDatum *s,
                                     const ModelParams &model_params);
double BulkLogLikWithCopyNumberProfile(size_t region,
                                       const CloneTreeNode *node,
                                       const BulkDatum *datum,
                                       const ModelParams &model_params);
double BulkLogLikWithGenotype(size_t region,
                              const CloneTreeNode *node,
                              const BulkDatum *datum,
                              const ModelParams &model_params);
double ZeroBulkLikelihood(const CloneTreeNode *node,
                          const BulkDatum *datum,
                          const ModelParams &model_params);

#endif /* clone_node_hpp */
