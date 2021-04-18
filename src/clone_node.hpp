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
    vector<double> clone_freqs_;
    vector<double> cellular_prevs_;

public:
    CloneTreeNodeParam(size_t region_count);
    inline const vector<double> &GetCloneFreqs() const { return clone_freqs_; }
    inline const vector<double> &GetCellularPrevalences() const { return cellular_prevs_; }
    inline double GetCloneFreqAtRegion(size_t region) const { return clone_freqs_[region]; }
    inline double GetCellularPrevalenceAtRegion(size_t region) const { return cellular_prevs_[region]; }
    void SetCloneFrequencyAtRegion(size_t idx, double val);
    void SetCellularPrevalenceAtRegion(size_t idx, double val);
    void SetRootParameters(const gsl_rng *random);
    void SetRootParameters();
    bool IsConsistent() const;
    size_t GetRegionCount() const {
        return clone_freqs_.size();
    }
    string GetCloneFreqsAsString() const;
    string GetCellularPrevsAsString() const;
};

class CloneTreeNode
{
    CloneTreeNodeParam param_;

    vector<size_t> name_;
    double nu_ = 0.0;
    CloneTreeNode *parent_node_ = 0;
    unordered_set<const BulkDatum *> data_;
    vector<double> sc_cache_;

    // map: child branch idx -> pair<psi_stick, Node *>
    unordered_map<size_t, pair<double, CloneTreeNode *> > idx2child_;

    // clone-specific copy number per each gene
    vector<size_t> cn_profile;

    bool operator==(const CloneTreeNode &other) const;

    // change the last part of the node's name to j -- used extensively by reorder_sticks
    void EditName(size_t j);
    
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
    double GetCache(size_t c);
    string GetName() const;
    inline const vector<size_t> &GetNameVector() const { return name_; };
    inline size_t GetDepth() const { return name_.size() - 1; }
    size_t DataCount() const;
    double GetNuStick() const;
    CloneTreeNode *GetParentNode() const;
    const unordered_set<const BulkDatum *> &GetData() const;
    const pair<double, CloneTreeNode *> &GetChild(size_t child_idx) const;

    const vector<size_t> &getCnProfile() const;

    void setCnProfile(const vector<size_t> &cnProfile);

    size_t GetChildrenCount() const;
    bool IsRoot() const;
    bool IsLeaf() const;
    unordered_map<size_t, pair<double, CloneTreeNode *> > &GetIdx2Child();

    // setters
    void SetNuStick(double nu);
    void SetPsiStick(size_t child_idx, double psi);
    
    // add/remove datum
    void AddDatum(BulkDatum *datum);
    void RemoveDatum(BulkDatum *datum);
    bool ContainsDatum(BulkDatum *datum) const;
    
    // operations on sticks
    void Cull(unordered_set<size_t> &cull_list);
    void ResampleStickOrder(const gsl_rng *random, const ModelParams &params);
    void ResetChildrenNames();
    
    // Returns true if this node is ancestor of other or other == this.
    bool IsAncestorOf(CloneTreeNode *other);
    bool IsDescendantOf(CloneTreeNode *other);
    
    // identify the branch that contains u, return the corresponding child Node
    CloneTreeNode *LocateChild(const gsl_rng *random, double &u, const ModelParams &hyper_params);
    void InitializeChild(const gsl_rng *random,
                         const ModelParams &params);
    
    static string ConstructNodeString(string curr_node_str, size_t branch);
    static string RetrieveParentString(string curr_node_str);
    
    static void BreadthFirstTraversal(CloneTreeNode *root, vector<CloneTreeNode *> &ret, bool non_empty = false);
    static CloneTreeNode *FindNode(const gsl_rng *random, double u, CloneTreeNode *root, const ModelParams &params);
    // Return all data from ancestors including node itself.
    static void GetDataset(CloneTreeNode *node, unordered_set<const BulkDatum *> &dataset);
    static void GetClusterLabels(CloneTreeNode *root,
                                 const vector<BulkDatum *> &data,
                                 vector<unsigned int> &cluster_labels);
    static void Datum2Node(vector<CloneTreeNode *> &all_nodes,
                           unordered_map<const BulkDatum *, CloneTreeNode *> &datum2node);
    static gsl_matrix *GetAncestralMatrix(CloneTreeNode *root,
                                          const vector<BulkDatum *> &bulk_data,
                                          const unordered_map<const BulkDatum *, CloneTreeNode *> &datum2node);

    static bool Less(CloneTreeNode *l, CloneTreeNode *r);
    friend bool operator<(const CloneTreeNode& lhs, const CloneTreeNode& rhs);
    friend bool operator> (const CloneTreeNode& lhs, const CloneTreeNode& rhs){ return rhs < lhs; }
    friend bool operator<=(const CloneTreeNode& lhs, const CloneTreeNode& rhs){ return !(lhs > rhs); }
    friend bool operator>=(const CloneTreeNode& lhs, const CloneTreeNode& rhs){ return !(lhs < rhs); }

    inline double GetCellularPrevs(size_t region) const {
        return param_.GetCellularPrevalenceAtRegion(region);
    }
    inline double GetCloneFreqs(size_t region) const {
        return param_.GetCloneFreqAtRegion(region);
    }
    CloneTreeNodeParam &NodeParameter();
    void SampleNuStick(const gsl_rng *random,
                       const ModelParams &params);
    void SampleNodeParameters(const gsl_rng *random,
                                const ModelParams &params);
    CloneTreeNode *SpawnChild(double psi);
    string Print();

    void SetCloneFrequencyAtRegion(size_t region, double new_val);
    void SetCellularPrevalenceAtRegion(size_t region, double new_val);
    void SetCloneFrequencyVector(vector<double> &vec);
    void SetCellularPrevalenceVector(vector<double> &vec);

    static CloneTreeNode *CreateRootNode(size_t region_count);
    static void RetrieveLoci(CloneTreeNode *node, unordered_set<Locus> &ret);
};

double ScLikelihood(size_t loci_idx,
                    const BulkDatum *s,
                    const SingleCellData *sc,
                    bool has_snv,
                    const ModelParams &model_params);
double ScLikelihoodWithDropout(size_t loci_idx,
                               const BulkDatum *bulk,
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
