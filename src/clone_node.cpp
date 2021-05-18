//
//  clone_node.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-06-27.
//

#include "clone_node.hpp"

#include <iomanip>
#include <unordered_set>

#include "numerical_utils.hpp"

CloneTreeNodeParam::CloneTreeNodeParam(size_t region_count)
{
    clone_freqs_.resize(region_count);
    cellular_prevs_.resize(region_count);
}

void CloneTreeNodeParam::SetCloneFrequencyAtRegion(size_t idx, double val)
{
    this->clone_freqs_[idx] = val;
}
void CloneTreeNodeParam::SetCellularPrevalenceAtRegion(size_t idx, double val)
{
    this->cellular_prevs_[idx] = val;
}
bool CloneTreeNodeParam::IsConsistent() const
{
    for (size_t i = 0; i < clone_freqs_.size(); i++) {
        if (clone_freqs_[i] > cellular_prevs_[i]) {
            return false;
        }
    }
    return true;
}

void CloneTreeNodeParam::SetRootParameters(const gsl_rng *random)
{
    // FIXME not actually random, this implementation is identical
    //  to `CloneTreeNodeParam::SetRootParameters()`
    double curr_cellular_prev = 1.0;
    for (size_t i = 0; i < GetRegionCount(); i++) {
        SetCellularPrevalenceAtRegion(i, curr_cellular_prev);
        SetCloneFrequencyAtRegion(i, curr_cellular_prev);
    }
}

// Initialize root cellular prevalence and clone frequency both to 1.
void CloneTreeNodeParam::SetRootParameters()
{
    for (size_t i = 0; i < GetRegionCount(); i++) {
        SetCellularPrevalenceAtRegion(i, 1.0);
        SetCloneFrequencyAtRegion(i, 1.0);
    }
}


string CloneTreeNodeParam::GetCloneFreqsAsString() const
{
    string str = "(";
    for (size_t i = 0; i < clone_freqs_.size() - 1; i++) {
        str += to_string(clone_freqs_[i]);
        str += ",";
    }
    str += to_string(clone_freqs_[clone_freqs_.size() - 1]);
    str += ")";
    return str;
}

string CloneTreeNodeParam::GetCellularPrevsAsString() const
{
    string str = "(";
    for (size_t i = 0; i < cellular_prevs_.size() - 1; i++) {
        str += to_string(cellular_prevs_[i]);
        str += ",";
    }
    str += to_string(cellular_prevs_[cellular_prevs_.size() - 1]);
    str += ")";
    return str;
}

CloneTreeNode::CloneTreeNode(size_t region_count) :
param_(region_count)
{
    parent_node_ = 0;
    this->name_.push_back(0);
}

CloneTreeNode::CloneTreeNode(size_t child_idx,
                             CloneTreeNode *parent) :
param_(parent->param_.GetRegionCount()),
parent_node_(parent)
{
    if (parent == 0) {
        cerr << "Error: parent node cannot be null.\n";
        exit(-1);
    }
    this->name_ = parent->name_; // copy the name vector
    this->name_.push_back(child_idx);
}

CloneTreeNode::~CloneTreeNode() {
}

CloneTreeNodeParam &CloneTreeNode::NodeParameter()
{
    return param_;
}

void CloneTreeNode::SampleNodeParameters(const gsl_rng *random, const ModelParams &params)
{
    if (parent_node_ == 0) {
        this->param_.SetRootParameters(random);
    } else {
        CloneTreeNode *parent_node = parent_node_;
        for (size_t i = 0; i < param_.GetRegionCount(); i++) {
            double parent_clone_freq = parent_node->param_.GetCloneFreqAtRegion(i);
            double curr_cellular_prev = uniform(random) * parent_clone_freq;
            this->param_.SetCellularPrevalenceAtRegion(i, curr_cellular_prev);
            this->param_.SetCloneFrequencyAtRegion(i, curr_cellular_prev);
            parent_node->param_.SetCloneFrequencyAtRegion(i, parent_clone_freq - curr_cellular_prev);
        }
    }
}

void CloneTreeNode::RetrieveLoci(CloneTreeNode *node, unordered_set<Locus> &ret)
{
    // traverse to the root to retrieve all SNVs
    CloneTreeNode *v = node;
    while (v != 0) {
        unordered_set<const BulkDatum *> data = v->GetData();
        for (const BulkDatum *datum : data) {
            ret.insert(datum->GetLocus());
        }
        v = v->GetParentNode();
    }
}

double get_branch_length(CloneTreeNode *child, CloneTreeNode *parent)
{
    // compute the sum of branch length from parent to child
    double len = 0.0;
    while (child != parent) {
        child = child->GetParentNode();
        len++;
        if (child == 0) {
            // reached root, child and parent are not in a ancestral relationship
            return -1;
        }
    }
    return len;
}

CloneTreeNode *CloneTreeNode::SpawnChild(double psi)
{
    size_t j = GetIdx2Child().size();
    //CloneTreeNode *child = new CloneTreeNode(j, this, tssb);
    auto child = new CloneTreeNode(j, this);
    idx2child_[j] = make_pair(psi, child);
    return child;
}

string CloneTreeNode::Print()
{
    string psi_stick_str = "(";
    unordered_map<size_t, pair<double, CloneTreeNode *> > &children = GetIdx2Child();
    for (size_t i = 0; i < GetChildrenCount(); i++) {
        psi_stick_str += to_string(children[i].first);
        if (i < GetChildrenCount() - 1) {
            psi_stick_str += ", ";
        }
    }
    psi_stick_str += ")";
    // print the CloneTreeNode name, nu-stick, CloneTreeNode params
    string ret = "[" + (parent_node_ == 0 ? "root" : this->GetName()) + ", ";
    ret += "phi=" + param_.GetCellularPrevsAsString() + ", ";
    ret += "eta=" + param_.GetCloneFreqsAsString() + ", ";
    ret += "nu=" + to_string(nu_) + ", ";
    ret += "psi=" + psi_stick_str + ", ";
    
    ret += "data=( ";
    for (const BulkDatum *datum : GetData()) {
        ret += datum->GetId() + " ";
    }
    ret += ")";
    ret += "]";

    return ret;
}

CloneTreeNode *CloneTreeNode::CreateRootNode(size_t region_count)
{
    auto node = new CloneTreeNode(region_count);
    return node;
}

void CloneTreeNode::SetCloneFrequencyAtRegion(size_t region, double new_val)
{
    param_.SetCloneFrequencyAtRegion(region, new_val);
}

void CloneTreeNode::SetCellularPrevalenceAtRegion(size_t region, double new_val)
{
    param_.SetCellularPrevalenceAtRegion(region, new_val);
}

string CloneTreeNode::RetrieveParentString(string curr_node_str)
{
    // split the string using "_"
    vector<string> fields;
    boost::split(fields, curr_node_str, boost::is_any_of("_"));
    string ret = "";
    size_t num_iter = fields.size() - 1;
    for (size_t i = 0; i < num_iter; i++)
    {
        ret += fields[i];
        if (i < num_iter - 1) {
            ret += "_";
        }
    }
    return ret;
}

string CloneTreeNode::ConstructNodeString(string curr_node_str, size_t branch)
{
    if (curr_node_str == "")
        return to_string(branch);
    
    return curr_node_str + "_" + to_string(branch);
}

CloneTreeNode *CloneTreeNode::FindNode(const gsl_rng *random, double u, CloneTreeNode *root, const ModelParams &params)
{
    CloneTreeNode *node = root;
    double nu = 0.0;
    while (true) {
        nu = node->GetNuStick();
        
        if (u < nu) {
            // found the node!!
            break;
        }
        
        // Shrink u relative to the remaining stick (1 - nu)
        u = (u - nu) / (1 - nu);
        
        // find sub branch by enumerating over the branching sticks
        node = node->LocateChild(random, u, params);
    }
    return node; // TODO: consider returning the path for the purposes of better estimation of the parameters
}

string CloneTreeNode::GetName() const
{
    // concatenate vector into string
    string str = "";
    for (size_t i = 0; i < name_.size(); i++) {
        str = this->ConstructNodeString(str, name_[i]);
    }
    return str;
}

void CloneTreeNode::EditName(size_t j)
{
    size_t n = this->name_.size();
    for (size_t i = 0; i < n - 1; i++) {
        name_[i] = parent_node_->name_[i];
    }
    name_[n-1] = j;
}

// getters
CloneTreeNode *CloneTreeNode::GetParentNode() const
{
    return parent_node_;
}

const unordered_set<const BulkDatum *> &CloneTreeNode::GetData() const
{
    return data_;
}

void CloneTreeNode::GetDataset(CloneTreeNode *node,
                            unordered_set<const BulkDatum *> &dataset)
{
    // trace up to the root node to get all SNVs
    while (node != 0) {
        dataset.insert(node->GetData().begin(), node->GetData().end());
        node = node->GetParentNode();
    }
}

size_t CloneTreeNode::DataCount() const
{
    return data_.size();
}

unordered_map<size_t, pair<double, CloneTreeNode *> > &CloneTreeNode::GetIdx2Child()
{
    return idx2child_;
}

double CloneTreeNode::GetNuStick() const
{
    return this->nu_;
}

size_t CloneTreeNode::GetChildrenCount() const
{
    return idx2child_.size();
}

bool CloneTreeNode::IsRoot() const
{
    return parent_node_ == 0;
}

bool CloneTreeNode::IsLeaf() const
{
    return (GetChildrenCount() == 0);
}

const pair<double, CloneTreeNode *> &CloneTreeNode::GetChild(size_t child_idx) const
{
    if (child_idx < idx2child_.size()) {
        return idx2child_.at(child_idx);
    }
    cerr << "Error: child index out of bounds." << endl;
    exit(-1);
}

void CloneTreeNode::SetNuStick(double nu)
{
    this->nu_ = nu;
}

void CloneTreeNode::SetPsiStick(size_t child_idx, double psi)
{
    pair<double, CloneTreeNode *> &child = idx2child_.at(child_idx);
    child.first = psi;
}

void CloneTreeNode::AddDatum(BulkDatum *datum)
{
    data_.insert(datum);
}

void CloneTreeNode::RemoveDatum(BulkDatum *datum)
{
    data_.erase(datum);
}

bool CloneTreeNode::ContainsDatum(BulkDatum *datum) const
{
    return (data_.count(datum) > 0);
}

void CloneTreeNode::Cull(unordered_set<size_t> &cull_list)
{
    unordered_map<size_t, pair<double, CloneTreeNode *> > new_idx2child;
    size_t idx = 0;
    for (size_t j = 0; j < idx2child_.size(); j++)
    {
        CloneTreeNode *node = idx2child_[j].second;
        if (cull_list.count(j) == 0) {
            new_idx2child[idx] = idx2child_[j];
            idx++;
        } else {
            //cout << get_name() << " cull child " << j << endl;
            idx2child_[j].second = 0;
            delete node;
        }
    }
    idx2child_ = new_idx2child;
}

void CloneTreeNode::ResetChildrenNames()
{
    // clear the children info
    unordered_map<size_t, pair<double, CloneTreeNode *> > new_map;
    
    CloneTreeNode *child = 0;
    size_t n_children = idx2child_.size();
    size_t idx = 0;
    for (size_t i = 0; i < n_children; i++) {
        child = idx2child_[i].second;
        if (child == 0) {
            //idx2child.erase(i); // this is not necessarily since we will overwrite idx2child
        } else {
            child->EditName(idx);
            new_map[idx] = idx2child_[i];
            idx++;
        }
    }
    idx2child_ = new_map;
}

void CloneTreeNode::ResampleStickOrder(const gsl_rng *random, const ModelParams &params)
{
    if (idx2child_.size() == 0)
        return;

    vector<double> unnorm_w(idx2child_.size());
    vector<double> intervals(idx2child_.size());
    double cum_prod = 1.0;
    for (size_t i = 0; i < idx2child_.size(); i++)
    {
        double ww = idx2child_[i].first;
        unnorm_w[i] = ww * cum_prod;
        cum_prod *= (1 - ww);
        intervals[i] = (1 - cum_prod);
    }
    size_t N = unnorm_w.size();
    double norm = intervals[N - 1];
    vector<size_t> new_order;
    vector<bool> represented(N, true);
    int idx;
    while (new_order.size() < N)
    {
        idx = multinomial(random, unnorm_w, norm);
        if (represented[idx]) {
            new_order.push_back(idx);
            represented[idx] = false;
            norm -= unnorm_w[idx];
            unnorm_w[idx] = 0.0;
        }
    }
    unordered_map<size_t, pair<double, CloneTreeNode *> > new_map;
    for (size_t i = 0; i < new_order.size(); i++) {
        idx2child_.at(new_order.at(i)).second->EditName(i);
        new_map[i] = idx2child_.at(new_order.at(i));
    }
    idx2child_ = new_map;
}

void CloneTreeNode::SampleNuStick(const gsl_rng *random,
                                  const ModelParams &params) {
    double nu_stick = bounded_beta(random, 1.0, params.ComputeAlpha(GetNameVector()));
    SetNuStick(nu_stick);

}

void CloneTreeNode::InitializeChild(const gsl_rng *random,
                                const ModelParams &params)
{
    double psi_j = bounded_beta(random, 1, params.GetGamma());
    CloneTreeNode *child = this->SpawnChild(psi_j);
    child->SampleNuStick(random, params);
    child->SampleNodeParameters(random, params);
}

CloneTreeNode *CloneTreeNode::LocateChild(const gsl_rng *random, double &u, const ModelParams &params)
{
    double cum_prod = 1;
    for (size_t i = 0; i < idx2child_.size(); i++) {
        cum_prod *= (1 - idx2child_[i].first);
    }
    while (u > (1 - cum_prod)) {
        InitializeChild(random, params);
        cum_prod *= (1 - idx2child_[idx2child_.size()-1].first);
    }

    cum_prod = 1.0;
    vector<double> intervals(idx2child_.size() + 1);
    intervals[0] = 0.0;
    for (size_t j = 0; j < idx2child_.size(); j++) {
        cum_prod *= (1 - idx2child_.at(j).first);
        intervals[j+1] = (1 - cum_prod);
    }
    for (size_t j = 0; j < idx2child_.size(); j++) {
        if (u < intervals[j+1]) {
            u = (u - intervals[j])/(intervals[j+1] - intervals[j]);
            return idx2child_[j].second;
        }
    }

    cerr << "Error in locate_child.\n";
    exit(-1);
}

bool CloneTreeNode::Less(CloneTreeNode *lhs, CloneTreeNode *rhs)
{
    size_t lhs_size = lhs->name_.size();
    size_t rhs_size = rhs->name_.size();
    if (rhs_size == 0 && lhs_size == 0) {
        return false;
    } else if (lhs_size == 0) {
        return true;
    } else if (rhs_size == 0) {
        return false;
    }
    
    size_t len = min(lhs_size, rhs_size);
    for (size_t i = 0; i < len; i++) {
        if (lhs->name_[i] < rhs->name_[i]) {
            return true;
        } else if (lhs->name_[i] > rhs->name_[i]) {
            return false;
        }
    }
    
    // if we are here, the two strings are same up to len
    // then, we determine the order one is shorter than the other
    if (lhs_size < rhs_size) {
        return true;
    } else {
        return false;
    }
}

/**
 * Breadth first traversal starting from `root_node`.
 *
 * @param[in] root_node
 * @param[out] ret sequence vector of traversed nodes
 * @param non_empty if true, only nodes with SNVs assigned to it are selected
 */
void CloneTreeNode::BreadthFirstTraversal(CloneTreeNode *root_node, vector<CloneTreeNode *> &ret, bool non_empty)
{
    // get nodes in the subtree rooted at root
    if (ret.size() > 0) {
        cout << "Warning: get_all_nodes() is going to clear vector ret." << endl;
        ret.clear();
    }
    if (root_node == 0) {
        cout << "Warning: root node is null." << endl;
        return;
    }
    
    set<CloneTreeNode *> nodes;
    CloneTreeNode *node = 0;
    
    queue<CloneTreeNode *> q;
    q.push(root_node);
    while (!q.empty())
    {
        node = q.front();
        q.pop();
        if (!nodes.count(node)) {
            if (!non_empty) {
                // complete bft
                ret.push_back(node);
                nodes.insert(node);
            } else {
                // select only nodes with SNVs assigned to it
                if (node->DataCount() > 0) {
                    ret.push_back(node);
                    nodes.insert(node);
                }
            }
        }
        unordered_map<size_t, pair<double, CloneTreeNode *> > &children = node->GetIdx2Child();
        for (size_t i = 0; i < children.size(); i++)
        {
            CloneTreeNode *child = children[i].second;
            q.push(child);
        }
    }
}

void CloneTreeNode::GetClusterLabels(CloneTreeNode *root,
                                   const vector<BulkDatum *> &data,
                                   vector<unsigned int> &cluster_labels)
{
    vector<CloneTreeNode *> all_nodes;
    CloneTreeNode::BreadthFirstTraversal(root, all_nodes, false);
    unordered_map<const BulkDatum *, CloneTreeNode *> datum2node;
    Datum2Node(all_nodes, datum2node);
    
    // determine the classes and clusters
    // assign class label from 0, ..., N = num nodes
    unordered_map<CloneTreeNode *, size_t> node2class;
    size_t num_classes = all_nodes.size();
    for (size_t c = 0; c < num_classes; c++) {
        node2class[all_nodes[c]] = c;
    }
    for (size_t i = 0; i < data.size(); i++) {
        CloneTreeNode *v = datum2node[data[i]];
        size_t c = node2class[v];
        cluster_labels.push_back(c);
    }
}

bool CloneTreeNode::IsAncestorOf(CloneTreeNode *other) {
    if (this->name_.size() > other->name_.size()) {
        return false;
    }
    for (size_t i = 0; i < name_.size(); i++) {
        if (name_.at(i) != other->name_.at(i)) {
            return false;
        }
    }
    return true;
}

bool CloneTreeNode::IsDescendantOf(CloneTreeNode *other) {
    if (this->name_.size() < other->name_.size()) {
        return false;
    }
    for (size_t i = 0; i < other->name_.size(); i++) {
        if (name_.at(i) != other->name_.at(i)) {
            return false;
        }
    }
    return true;
}

bool CloneTreeNode::IsCacheAllocated(size_t cell_count)
{
    return (sc_cache_.size() == cell_count);
}

void CloneTreeNode::AllocateCache(size_t cell_count)
{
    if (sc_cache_.size() == 0) {
        sc_cache_.resize(cell_count);
    } else {
        if (cell_count != sc_cache_.size()) {
            cerr << "Cache size does not match the required size of " << cell_count << ".\n";
            exit(-1);
        }
    }
}

void CloneTreeNode::UpdateCache(size_t c, double val)
{
    sc_cache_[c] += val;
}

double CloneTreeNode::GetCache(size_t c)
{
    return sc_cache_[c];
}

void CloneTreeNode::Datum2Node(vector<CloneTreeNode *> &all_nodes,
                                     unordered_map<const BulkDatum *, CloneTreeNode *> &datum2node)
{
    for (auto it = all_nodes.begin(); it != all_nodes.end(); ++it) {
        for (const BulkDatum *datum : (*it)->GetData())
        {
            datum2node[datum] = (*it);
        }
    }
}

gsl_matrix *CloneTreeNode::GetAncestralMatrix(CloneTreeNode *root,
                                          const vector<BulkDatum *> &data,
                                          const unordered_map<const BulkDatum *, CloneTreeNode *> &datum2node)
{
    size_t N = data.size();
    vector<unordered_set<const BulkDatum *> > ancestors(N);
    for (size_t i = 0; i < N; i++) {
        BulkDatum *datum = data.at(i);
        CloneTreeNode *v = datum2node.at(datum);
        // trace up to the root and set the row of A
        while (v != root) {
            v = v->GetParentNode();
            unordered_set<const BulkDatum *> d = v->GetData();
            ancestors[i].insert(d.begin(), d.end());
        }
    }
    
    // A_ij = 1 if i occurs before j
    gsl_matrix *A = gsl_matrix_alloc(N, N);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            if (ancestors[j].count(data.at(i))) {
                gsl_matrix_set(A, i, j, 1);
            } else {
                gsl_matrix_set(A, i, j, 0);
            }
        }
    }
    
    return A;
}

bool CloneTreeNode::operator==(const CloneTreeNode &other) const {
    if (this->parent_node_ != other.parent_node_)
        return false;
    if (name_.size() != other.name_.size())
        return false;
    size_t n = name_.size();
    return (name_[n-1] == other.name_[n-1]);
}

bool operator<(const CloneTreeNode &lhs, const CloneTreeNode &rhs) {
    vector<string> lhs_arr, rhs_arr;
    boost::split(lhs_arr, lhs.name_, boost::is_any_of("_"));
    boost::split(rhs_arr, rhs.name_, boost::is_any_of("_"));
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

const vector<size_t> &CloneTreeNode::getCnProfile() const {
    return cn_profile;
}

void CloneTreeNode::setCnProfile(const vector<size_t> &cnProfile) {
    cn_profile = cnProfile;
}

const shared_ptr<vector<Bin>> &CloneTreeNode::getBins() const {
    return bins;
}

void CloneTreeNode::setBins(const shared_ptr<vector<Bin>> &bins) {
    CloneTreeNode::bins = bins;
}

const vector<size_t> &CloneTreeNode::getGeneCnProfile() const {
    return geneCnProfile;
}

void CloneTreeNode::setGeneCnProfile(const vector<size_t> &geneCnProfile) {
    CloneTreeNode::geneCnProfile = geneCnProfile;
}

// TODO incorporate gene expression data in likelihood functions

double ScLikelihoodWithDropout(size_t loci_idx,
                    const BulkDatum *bulk,
                    const SingleCellData *sc,
                    bool has_snv,
                    const ModelParams &model_params) {
    double log_lik = 0.0;

    // if sc has mutation s, then there are 3 cases
    // 1. Bi-allelic,
    // 2. Bursty for variant,
    // 3. Bursty for reference (dropout).
    size_t var_reads = sc->GetVariantReads(loci_idx);
    size_t total_reads = sc->GetTotalReads(loci_idx);
    if (total_reads == 0) {
        return 0.0;
    }
    if (has_snv) {
        auto locus = bulk->GetLocus();
        double alpha = locus.get_alpha();
        double beta = locus.get_beta();
        double log_lik_biallelic = log_beta_binomial_pdf(var_reads,
                                                         total_reads,
                                                         alpha,
                                                         beta);
        log_lik_biallelic += log(model_params.GetScBiallelicProportion());
        double log_lik_bursty_variant = log_beta_binomial_pdf(var_reads,
                                                              total_reads,
                                                              model_params.GetScBurstyDistributionAlpha0(),
                                                              model_params.GetScBurstyDistributionBeta0());
        log_lik_bursty_variant += log(model_params.GetScBurstyVariantProportion());
        double log_lik_dropout = log_beta_binomial_pdf(var_reads,
                                                      total_reads,
                                                      model_params.GetScDropoutDistributionAlpha0(),
                                                      model_params.GetScDropoutDistributionBeta0());
        log_lik_dropout += log(model_params.GetScDropoutProportion());
        log_lik = log_add(log_lik_dropout, log_lik_biallelic);
        log_lik = log_add(log_lik, log_lik_bursty_variant);
    } else {
        log_lik = log_beta_binomial_pdf(var_reads,
                                        total_reads,
                                        model_params.GetSequencingError(),
                                        1 - model_params.GetSequencingError());
    }
    
    return log_lik;

}

double ScLikelihood(size_t loci_idx,
                    const BulkDatum *bulk,
                    const SingleCellData *sc,
                    bool has_snv,
                    const ModelParams &model_params) {
    double log_lik = 0.0;

    // if sc has mutation s, then there are 2 cases
    // 1. Bi-allelic,
    // 2. Bursty distribution.
    size_t var_reads = sc->GetVariantReads(loci_idx);
    size_t total_reads = sc->GetTotalReads(loci_idx);
    if (total_reads == 0) {
        return 0.0;
    }
    if (has_snv) {
        auto locus = bulk->GetLocus();
        double alpha = locus.get_alpha();
        double beta = locus.get_beta();
        double log_lik_biallelic = log_beta_binomial_pdf(var_reads,
                                                          total_reads,
                                                          alpha,
                                                          beta);
        log_lik_biallelic += log(1-locus.GetBurstyProbability());
        double log_lik_dropout = log_beta_binomial_pdf(var_reads,
                                                       total_reads,
                                                       model_params.GetScBurstyDistributionAlpha0(),
                                                       model_params.GetScBurstyDistributionBeta0());
        log_lik_dropout += log(locus.GetBurstyProbability());
        log_lik = log_add(log_lik_biallelic, log_lik_dropout);
    } else {
        log_lik = log_beta_binomial_pdf(var_reads,
                                        total_reads,
                                        model_params.GetSequencingError(),
                                        1 - model_params.GetSequencingError());
    }

    return log_lik;
}

/**
 * Log likelihood of a single cell for all genes and SNVs in a bin
 *
 * @param bin_idx
 * @param bulk_data
 * @param sc
 * @return
 */
double SCLogLikWithCopyNumber(size_t bin_idx, const vector<BulkDatum *> &bulk_data, const SingleCellData *sc,
                              const CloneTreeNode *node, const vector<Gene *> &geneSet, const ModelParams &modelParams) {
    // precompute the normalization factor for the mean
    double normFactor = computeNormFactor(node);

    // for each gene in the bin
    double binLogLik = 0.0;
    for (auto geneIdx: node->getBins()->at(bin_idx).getGeneIdxs()) {
        double geneLogLik = DOUBLE_NEG_INF;
        int total_cn = node->getCnProfile()[bin_idx];
        for (int e = 1; e <= total_cn; ++e) {
            // marginalize the number of expressed copies `e`
            double allImbLogLik = 0.0; // when no snv is present in the gene
            for (auto loc_idx: sc->GetLoci()) {
                // find snv in gene (if any)
                // TODO can be optimized saving SNVs in bins/genes
                if (geneSet[geneIdx]->containsSNV(bulk_data[loc_idx])) {
                    // if SNV is present in the gene, compute the allelic imbalance prob
                    // marginalizing over variant copy number
                    size_t var_reads = sc->GetVariantReads(loc_idx);
                    size_t total_reads = sc->GetTotalReads(loc_idx);
                    allImbLogLik = DOUBLE_NEG_INF;
                    for (int v = 0; v < e; ++v) {
                        // compute sc snv probability
                        double binomLogLik = gsl_sf_lnchoose(e, v) +
                                gsl_ran_binomial_pdf(var_reads, (double) v / e, total_reads);
                        allImbLogLik = log_add(binLogLik, allImbLogLik);
                    }
                    break;
                }
            }
            // compute gene expression prob (ZINB dist) multiplied by the prior for `e` (binomial)
            double exprLogLik = log(gsl_ran_binomial_pdf(e, geneSet[geneIdx]->getGeneCopyProb(), total_cn)); // binom
            // clonealign formula
            double mean = sc->getDepthSize() * geneSet[geneIdx]->getPerCopyExpr() * e / normFactor;
            switch (modelParams.getExprModel()) {
                case POISSON: {
                    exprLogLik += log_poisson_pdf(sc->getExprReads()[geneIdx], mean); // gene expr likelihood
                    break;
                }
                case NEG_BINOM: {
                    exprLogLik += log_negative_binomial_pdf(sc->getExprReads()[geneIdx], mean,
                                               geneSet[geneIdx]->getNbInvDispersion());
                    break;
                }
                case ZIP: {
                    exprLogLik += log_zip_pdf(sc->getExprReads()[geneIdx], mean,
                                                            sc->getZeroInflationProbs()[geneIdx]);
                    break;
                }
                case ZINB: {
                    exprLogLik += log_zinb_pdf(sc->getExprReads()[geneIdx], mean,
                                               geneSet[geneIdx]->getNbInvDispersion(),
                                               sc->getZeroInflationProbs()[geneIdx]);
                    break;
                }
            }
            exprLogLik += allImbLogLik; // add allelic imbalance component (can be 0)

            geneLogLik = log_add(geneLogLik, exprLogLik); // addend of the outer sum
        }
        // sum all the log likelihoods of each gene in the bin
        binLogLik += geneLogLik;
    }
    return binLogLik;
}

/**
 * Compute the normalization factor required for mean computation
 * in the clonealign formula for single cell $c$, node $v$
 *  $\sum_{g'=1}^G \mu_{g'} D_{b(g')v} \delta_{g'}$
 *
 * @param node node to which the single cell is assigned
 * @return normalization factor (just the denominator)
 */
double computeNormFactor(const CloneTreeNode *node) {
    double sum = 0.0;
    for (int b = 0; b < node->getBins()->size(); ++b) {
        for (auto gene: node->getBins()->at(b).getGenes()) {
            sum += gene->getPerCopyExpr() * node->getCnProfile()[b] * gene->getGeneCopyProb();
        }
    }
    return sum;
}

double BulkLogLikWithTotalCopyNumberByRegion(const CloneTreeNode *node,
                                             const ModelParams &model_params,
                                             size_t var_reads,
                                             size_t total_reads,
                                             size_t total_cn,
                                             double cell_prev) {
    if (total_reads == 0) {
        return 0.0;
    }
    
    double seq_err = model_params.GetSequencingError();
    double phi = cell_prev;
    double xi = 0.0;
    
    if (total_cn == 0) {
        xi = seq_err;
        return log(gsl_ran_binomial_pdf(var_reads,
                                        xi,
                                        total_reads));
    }

    double log_val;
    double log_prior_genotype;
    double log_lik = DOUBLE_NEG_INF;
    double log_prior_norm = log(1 - gsl_ran_binomial_pdf(0, model_params.GetVariantCopyProbability(), total_cn));
    // We marginalize over the number of variant copies
    for (size_t var_cn = 1; var_cn <= total_cn; var_cn++) {
        xi = (1 - phi) * seq_err;
        if (var_cn == total_cn) {
            xi += phi * (1 - seq_err);
        } else {
            xi += phi * var_cn / total_cn;
        }
        log_val = log_binomial_pdf(var_reads, xi, total_reads);
        log_prior_genotype = log_binomial_pdf(var_cn, model_params.GetVariantCopyProbability(), total_cn);
        log_prior_genotype -= log_prior_norm;
        log_lik = log_add(log_lik, log_val + log_prior_genotype);
    }
    return log_lik;
}

double BulkLogLikWithTotalCopyNumber(size_t region,
                                     const CloneTreeNode *node,
                                     const BulkDatum *datum,
                                     const ModelParams &model_params)
{
    auto var_reads = datum->GetVariantReadCount(region);
    auto total_reads = datum->GetReadCount(region);
    auto total_cn = datum->GetTotalCN(region);
    auto cell_prev = node->GetCellularPrevs(region);
    double log_lik = 0.0;
    log_lik += BulkLogLikWithTotalCopyNumberByRegion(node,
                                                     model_params,
                                                     var_reads,
                                                     total_reads,
                                                     total_cn,
                                                     cell_prev);
    return log_lik;
}

double BulkLogLikWithCopyNumberProfileByRegion(const CloneTreeNode *node,
                                               const ModelParams &model_params,
                                               size_t var_reads,
                                               size_t total_reads,
                                               const vector<double> &cn_probs,
                                               double cell_prev) {
    if (total_reads == 0) {
        return 0.0;
    }
    
    double seq_err = model_params.GetSequencingError();
    double xi = 0.0;
    double log_val;
    double log_prior_genotype;
    
    // Initialize log_lik for total_cn = 0.
    // When total_cn = 0, variant is only seen (1 - phi) of the cells as a result
    // of sequencing error.
    xi = (1 - cell_prev) * seq_err;
    double log_lik = log_binomial_pdf(var_reads,
                                      xi,
                                      total_reads);
    
    for (size_t total_cn = 0; total_cn < cn_probs.size(); total_cn++) {
        double log_prior_norm = log(1 - gsl_ran_binomial_pdf(0, model_params.GetVariantCopyProbability(), total_cn));
        double log_lik_cn = DOUBLE_NEG_INF;
        for (size_t var_cn = 1; var_cn <= total_cn; var_cn++) {
            xi = (1 - cell_prev) * seq_err;
            if (var_cn == total_cn) {
                xi += cell_prev * (1 - seq_err);
            } else if (var_cn == 0) {
                xi += cell_prev * seq_err;
            } else {
                xi += cell_prev * var_cn / total_cn;
            }
            log_val = log_binomial_pdf(var_reads,
                                       xi,
                                       total_reads);
            log_prior_genotype = log(gsl_ran_binomial_pdf(var_cn,
                                                          model_params.GetVariantCopyProbability(), total_cn));
            log_prior_genotype -= log_prior_norm;
            log_lik_cn = log_add(log_lik_cn, log_val + log_prior_genotype);
        }
        log_lik = log_add(log_lik, log_lik_cn + log(cn_probs[total_cn]));
    }
    return log_lik;
}

double BulkLogLikWithCopyNumberProfile(size_t region,
                                       const CloneTreeNode *node,
                                       const BulkDatum *datum,
                                       const ModelParams &model_params)
{
    auto var_reads = datum->GetVariantReadCount(region);
    auto total_reads = datum->GetReadCount(region);
    auto cell_prev = node->GetCellularPrevs(region);
    double log_lik = 0.0;
    auto cn_probs = datum->GetCNProbs(region);
    log_lik += BulkLogLikWithCopyNumberProfileByRegion(node,
                                                       model_params,
                                                       var_reads,
                                                       total_reads,
                                                       cn_probs,
                                                       cell_prev);
    return log_lik;
}

double BulkLogLikWithGenotypeByRegion(const CloneTreeNode *node,
                                      const ModelParams &model_params,
                                      size_t var_reads,
                                      size_t total_reads,
                                      size_t major_cn,
                                      size_t minor_cn,
                                      double cell_prev)
{
    if (major_cn < minor_cn) {
        cerr << "Error: major copy number < minor copy number.\n";
        exit(-1);
    }
    if (total_reads == 0) {
        if (var_reads > 0) {
            cerr << "Error in the data.\n";
            cerr << "Num reads: " << total_reads << ".\n";
            cerr << "Num variants: " << var_reads << ".\n";
            exit(-1);
        } else {
            return 0.0;
        }
    }

    double seq_err = model_params.GetSequencingError();
    
    size_t total_cn = major_cn + minor_cn;
    double xi = 0.0;
    
    if (total_cn == 0) {
        xi = seq_err;
        return log(gsl_ran_binomial_pdf(var_reads,
                                        xi,
                                        total_reads));
    }

    double log_val;
    double log_lik = DOUBLE_NEG_INF;
    
    // Marginalize over the var_cn over minor_cn.
    for (size_t var_cn = 1; var_cn <= minor_cn; var_cn++) {
        xi = (1 - cell_prev) * seq_err;
        if (var_cn == total_cn) {
            xi += cell_prev * (1 - seq_err);
        } else {
            xi += cell_prev * var_cn / total_cn;
        }
        log_val = log_binomial_pdf(var_reads, xi, total_reads);
        log_lik = log_add(log_lik, log_val);
    }
    // Marginalize over the var_cn over major_cn.
    for (size_t var_cn = 1; var_cn <= major_cn; var_cn++) {
        xi = (1 - cell_prev) * seq_err;
        if (var_cn == total_cn) {
            xi += cell_prev * (1 - seq_err);
        } else {
            xi += cell_prev * var_cn / total_cn;
        }
        log_val = log_binomial_pdf(var_reads, xi, total_reads);
        log_lik = log_add(log_lik, log_val);
    }
    log_lik -= log(total_cn);

    return log_lik;
}

double BulkLogLikWithGenotype(size_t region,
                              const CloneTreeNode *node,
                              const BulkDatum *datum,
                              const ModelParams &model_params)
{
    auto var_reads = datum->GetVariantReadCount(region);
    auto total_reads = datum->GetReadCount(region);
    auto cell_prevs = node->GetCellularPrevs(region);
    auto major_cns = datum->GetMajorCN(region);
    auto minor_cns = datum->GetMinorCN(region);
    double log_lik = BulkLogLikWithGenotypeByRegion(node,
                                                    model_params,
                                                    var_reads,
                                                    total_reads,
                                                    major_cns,
                                                    minor_cns,
                                                    cell_prevs);
    return log_lik;
}

double ZeroBulkLikelihood(const CloneTreeNode *node,
                          const BulkDatum *datum,
                          const ModelParams &model_params)
{
    return 0.0;
}

