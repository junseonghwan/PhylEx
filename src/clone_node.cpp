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
    clone_freqs.resize(region_count);
    cellular_prevs.resize(region_count);
}

void CloneTreeNodeParam::set_clone_freq(size_t idx, double val)
{
    this->clone_freqs[idx] = val;
}
void CloneTreeNodeParam::set_cellular_prev(size_t idx, double val)
{
    this->cellular_prevs[idx] = val;
}
void CloneTreeNodeParam::set_clone_freq(vector<double> &vec)
{
    this->clone_freqs.insert(clone_freqs.begin(), vec.begin(), vec.end());
}
void CloneTreeNodeParam::set_cellular_prev(vector<double> &vec)
{
    this->cellular_prevs.insert(cellular_prevs.begin(), vec.begin(), vec.end());
}

bool CloneTreeNodeParam::is_consistent()
{
    for (size_t i = 0; i < clone_freqs.size(); i++) {
        if (clone_freqs[i] > cellular_prevs[i]) {
            return false;
        }
    }
    return true;
}

void CloneTreeNodeParam::SetRootParameters()
{
    for (size_t i = 0; i < GetRegionCount(); i++) {
        clone_freqs[i] = 1.0;
        cellular_prevs[i] = 1.0;
    }
}

string CloneTreeNodeParam::GetCloneFreqsAsString()
{
    string str = "(";
    for (size_t i = 0; i < clone_freqs.size() - 1; i++) {
        str += to_string(clone_freqs[i]);
        str += ",";
    }
    str += to_string(clone_freqs[clone_freqs.size() - 1]);
    str += ")";
    return str;
}

string CloneTreeNodeParam::GetCellularPrevsAsString()
{
    string str = "(";
    for (size_t i = 0; i < cellular_prevs.size() - 1; i++) {
        str += to_string(cellular_prevs[i]);
        str += ",";
    }
    str += to_string(cellular_prevs[cellular_prevs.size() - 1]);
    str += ")";
    return str;
}

CloneTreeNode::CloneTreeNode(size_t region_count) :
param(region_count)
{
    parent_node = 0;
    this->name.push_back(0);
}

CloneTreeNode::CloneTreeNode(size_t child_idx,
                             CloneTreeNode *parent) :
param(parent->param.GetRegionCount()),
parent_node(parent)
{
    if (parent == 0) {
        cerr << "Error: parent node cannot be null.\n";
        exit(-1);
    }
    this->name = parent->name; // copy the name vector
    this->name.push_back(child_idx);
}

CloneTreeNode::~CloneTreeNode() {
}

CloneTreeNodeParam &CloneTreeNode::get_node_parameter()
{
    return param;
}

void CloneTreeNode::sample_node_parameters(const gsl_rng *random, const ModelParams &params,  CloneTreeNode *parent)
{
    if (parent_node == 0) {
        this->param.SetRootParameters();
    } else {
        CloneTreeNode *parent_node = (CloneTreeNode *)parent;
        for (size_t i = 0; i < param.GetRegionCount(); i++) {
            double parent_clone_freq = parent_node->param.get_clone_freqs(i);
            double curr_cellular_prev = uniform(random) * parent_clone_freq;
            this->param.set_cellular_prev(i, curr_cellular_prev);
            this->param.set_clone_freq(i, curr_cellular_prev);
            parent_node->param.set_clone_freq(i, parent_clone_freq - curr_cellular_prev);
        }
    }
}

void CloneTreeNode::RetrieveLoci(CloneTreeNode *node, unordered_set<Locus> &ret)
{
    // traverse to the root to retrieve all SNVs
    CloneTreeNode *v = node;
    while (v != 0) {
        unordered_set<const BulkDatum *> data = v->get_data();
        for (const BulkDatum *datum : data) {
            ret.insert(datum->GetLocus());
        }
        v = v->get_parent_node();
    }
}

double get_branch_length(CloneTreeNode *child, CloneTreeNode *parent)
{
    // compute the sum of branch length from parent to child
    double len = 0.0;
    while (child != parent) {
        child = child->get_parent_node();
        len++;
        if (child == 0) {
            // reached root, child and parent are not in a ancestral relationship
            return -1;
        }
    }
    return len;
}

CloneTreeNode *CloneTreeNode::spawn_child(double psi)
{
    size_t j = get_idx2child().size();
    //CloneTreeNode *child = new CloneTreeNode(j, this, tssb);
    CloneTreeNode *child = new CloneTreeNode(j, this);
    idx2child[j] = make_pair(psi, child);
    return child;
}

string CloneTreeNode::print()
{
    string psi_stick_str = "(";
    unordered_map<size_t, pair<double, CloneTreeNode *> > &children = get_idx2child();
    for (size_t i = 0; i < get_num_children(); i++) {
        psi_stick_str += to_string(children[i].first);
        if (i < get_num_children() - 1) {
            psi_stick_str += ", ";
        }
    }
    psi_stick_str += ")";
    // print the CloneTreeNode name, nu-stick, CloneTreeNode params
    string ret = "[" + (parent_node == 0 ? "root" : this->get_name()) + ", ";
    ret += "phi=" + param.GetCellularPrevsAsString() + ", ";
    ret += "eta=" + param.GetCloneFreqsAsString() + ", ";
    ret += "nu=" + to_string(nu) + ", ";
    ret += "psi=" + psi_stick_str + ", ";
    
    ret += "data=( ";
    for (const BulkDatum *datum : get_data()) {
        ret += datum->GetId() + " ";
    }
    ret += ")";
    ret += "]";

    return ret;
}

CloneTreeNode *CloneTreeNode::create_root_node(size_t region_count)
{
    auto node = new CloneTreeNode(region_count);
    return node;
}

void CloneTreeNode::set_clone_freq(size_t idx, double new_val)
{
    param.set_clone_freq(idx, new_val);
}

void CloneTreeNode::set_cellular_prev(size_t idx, double new_val)
{
    param.set_cellular_prev(idx, new_val);
}

string CloneTreeNode::get_parent_string(string curr_node_str)
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

string CloneTreeNode::form_node_string(string curr_node_str, size_t branch)
{
    if (curr_node_str == "")
        return to_string(branch);
    
    return curr_node_str + "_" + to_string(branch);
}

CloneTreeNode *CloneTreeNode::find_node(const gsl_rng *random, double u, CloneTreeNode *root, const ModelParams &params)
{
    CloneTreeNode *node = root;
    double nu = 0.0;
    while (true) {
        nu = node->get_nu_stick();
        
        if (u < nu) {
            // found the node!!
            break;
        }
        
        // Shrink u relative to the remaining stick (1 - nu)
        u = (u - nu) / (1 - nu);
        
        // find sub branch by enumerating over the branching sticks
        node = node->locate_child(random, u, params);
    }
    return node; // TODO: consider returning the path for the purposes of better estimation of the parameters
}

string CloneTreeNode::get_name() const
{
    // concatenate vector into string
    string str = "";
    for (size_t i = 0; i < name.size(); i++) {
        str = this->form_node_string(str, name[i]);
    }
    return str;
}

void CloneTreeNode::edit_name(size_t j)
{
    size_t n = this->name.size();
    for (size_t i = 0; i < n - 1; i++) {
        name[i] = parent_node->name[i];
    }
    name[n-1] = j;
}

// getters
CloneTreeNode *CloneTreeNode::get_parent_node() const
{
    return parent_node;
}

const unordered_set<const BulkDatum *> &CloneTreeNode::get_data() const
{
    return data;
}

void CloneTreeNode::GetDataset(CloneTreeNode *node,
                            unordered_set<const BulkDatum *> &dataset)
{
    // trace up to the root node to get all SNVs
    while (node != 0) {
        dataset.insert(node->get_data().begin(), node->get_data().end());
        node = node->get_parent_node();
    }
}

size_t CloneTreeNode::get_num_data() const
{
    return data.size();
}

unordered_map<size_t, pair<double, CloneTreeNode *> > &CloneTreeNode::get_idx2child()
{
    return idx2child;
}

double CloneTreeNode::get_nu_stick() const
{
    return this->nu;
}

size_t CloneTreeNode::get_num_children() const
{
    return idx2child.size();
}

bool CloneTreeNode::is_root() const
{
    return parent_node == 0;
}

bool CloneTreeNode::is_leaf() const
{
    return (get_num_children() == 0);
}

const pair<double, CloneTreeNode *> &CloneTreeNode::get_child(size_t child_idx) const
{
    if (child_idx < idx2child.size()) {
        return idx2child.at(child_idx);
    }
    cerr << "Error: child index out of bounds." << endl;
    exit(-1);
}

void CloneTreeNode::set_nu_stick(double nu)
{
    this->nu = nu;
}

void CloneTreeNode::set_psi_stick(size_t child_idx, double psi)
{
    pair<double, CloneTreeNode *> &child = idx2child.at(child_idx);
    child.first = psi;
}

void CloneTreeNode::add_datum(BulkDatum *datum)
{
    data.insert(datum);
}

void CloneTreeNode::remove_datum(BulkDatum *datum)
{
    data.erase(datum);
}

bool CloneTreeNode::contains_datum(BulkDatum *datum) const
{
    return (data.count(datum) > 0);
}

void CloneTreeNode::cull(unordered_set<size_t> &cull_list)
{
    unordered_map<size_t, pair<double, CloneTreeNode *> > new_idx2child;
    size_t idx = 0;
    for (size_t j = 0; j < idx2child.size(); j++)
    {
        CloneTreeNode *node = idx2child[j].second;
        if (cull_list.count(j) == 0) {
            new_idx2child[idx] = idx2child[j];
            idx++;
        } else {
            //cout << get_name() << " cull child " << j << endl;
            idx2child[j].second = 0;
            delete node;
        }
    }
    idx2child = new_idx2child;
}

void CloneTreeNode::reset_children_names()
{
    // clear the children info
    unordered_map<size_t, pair<double, CloneTreeNode *> > new_map;
    
    CloneTreeNode *child = 0;
    size_t n_children = idx2child.size();
    size_t idx = 0;
    for (size_t i = 0; i < n_children; i++) {
        child = idx2child[i].second;
        if (child == 0) {
            //idx2child.erase(i); // this is not necessarily since we will overwrite idx2child
        } else {
            child->edit_name(idx);
            new_map[idx] = idx2child[i];
            idx++;
        }
    }
    idx2child = new_map;
}

void CloneTreeNode::set_clone_freq(vector<double> &vec)
{
    param.set_clone_freq(vec);
}
void CloneTreeNode::set_cellular_prev(vector<double> &vec)
{
    param.set_cellular_prev(vec);
}


//void CloneTreeNode::reorder_sticks(const gsl_rng *random, const ModelParams &params)
//{
//    if (idx2child.size() == 0)
//        return;
//
//    vector<double> weights(idx2child.size());
//    unordered_set<int> represented;
//    double cum_prod = 1.0;
//    double represented_stick_length = 0.0;
//    for (size_t i = 0; i < idx2child.size(); i++)
//    {
//        weights[i] = idx2child[i].first * cum_prod;
//        represented_stick_length += weights[i];
//        cum_prod *= (1 - idx2child[i].first);
//        represented.insert(i);
//    }
//
//    // throw uniform darts until all represented sticks are sampled
//    size_t n_children = weights.size();
//    vector<size_t> new_order;
//    while (represented.size() > 0) {
//        double u = gsl_ran_flat(random, 0, 1);
//        while (u > represented_stick_length) {
//            // need to represent new children: i.e., draw psi sticks
//            double psi_j = bounded_beta(random, 1, params.get_gamma());
//            CloneTreeNode *child = this->spawn_child(psi_j);
//            double nu_stick = bounded_beta(random, 1.0, params.alpha(child->get_name_vec()));
//            child->set_nu_stick(nu_stick);
//            child->sample_node_parameters(random, params, this);
//
//            double ww = psi_j * cum_prod;
//            cum_prod *= (1 - psi_j);
//            weights.push_back(ww);
//            represented.insert(n_children);
//            represented_stick_length += ww;
//            n_children++;
//        }
//
//        vector<double> sub_weights;
//        vector<size_t> sub_indices;
//        double sum = 0.0;
//        for (size_t i = 0; i < n_children; i++) {
//            if (represented.count(i) > 0) {
//                sub_indices.push_back(i);
//                sub_weights.push_back(weights[i]);
//                sum += weights[i];
//            }
//        }
//
//        double cum = 0.0;
//        for (size_t i = 0; i < sub_weights.size(); i++) {
//            double norm_w = sub_weights.at(i)/sum;
//            if (u < cum + norm_w) {
//                represented.erase(sub_indices.at(i));
//                new_order.push_back(sub_indices.at(i));
//                break;
//            }
//            cum += norm_w;
//        }
//    }
//
//    unordered_map<size_t, pair<double, CloneTreeNode *> > new_map;
//    for (size_t i = 0; i < new_order.size(); i++) {
//        pair<double, CloneTreeNode *> &ret = idx2child[new_order.at(i)];
//        ret.second->edit_name(i);
//        new_map[i] = ret;
//    }
//    idx2child = new_map;
//}

void CloneTreeNode::reorder_sticks(const gsl_rng *random, const ModelParams &params)
{
    if (idx2child.size() == 0)
        return;

    vector<double> unnorm_w(idx2child.size());
    vector<double> intervals(idx2child.size());
    double cum_prod = 1.0;
    for (size_t i = 0; i < idx2child.size(); i++)
    {
        double ww = idx2child[i].first;
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
        idx2child.at(new_order.at(i)).second->edit_name(i);
        new_map[i] = idx2child.at(new_order.at(i));
    }
    idx2child = new_map;
}

void CloneTreeNode::InitializeChild(const gsl_rng *random,
                                const ModelParams &params)
{
    double psi_j = bounded_beta(random, 1, params.GetGamma());
    CloneTreeNode *child = this->spawn_child(psi_j);
    double nu_stick = bounded_beta(random, 1.0, params.ComputeAlpha(child->get_name_vec()));
    child->set_nu_stick(nu_stick);
    child->sample_node_parameters(random, params, this);
}

CloneTreeNode *CloneTreeNode::locate_child(const gsl_rng *random, double &u, const ModelParams &params)
{
    double cum_prod = 1;
    for (size_t i = 0; i < idx2child.size(); i++) {
        cum_prod *= (1 - idx2child[i].first);
    }
    while (u > (1 - cum_prod)) {
        InitializeChild(random, params);
        cum_prod *= (1 - idx2child[idx2child.size()-1].first);
    }

    cum_prod = 1.0;
    vector<double> intervals(idx2child.size() + 1);
    intervals[0] = 0.0;
    for (size_t j = 0; j < idx2child.size(); j++) {
        cum_prod *= (1 - idx2child.at(j).first);
        intervals[j+1] = (1 - cum_prod);
    }
    for (size_t j = 0; j < idx2child.size(); j++) {
        if (u < intervals[j+1]) {
            u = (u - intervals[j])/(intervals[j+1] - intervals[j]);
            return idx2child[j].second;
        }
    }

    cerr << "Error in locate_child.\n";
    exit(-1);
}

bool CloneTreeNode::less(CloneTreeNode *lhs, CloneTreeNode *rhs)
{
    size_t lhs_size = lhs->name.size();
    size_t rhs_size = rhs->name.size();
    if (rhs_size == 0 && lhs_size == 0) {
        return false;
    } else if (lhs_size == 0) {
        return true;
    } else if (rhs_size == 0) {
        return false;
    }
    
    size_t len = min(lhs_size, rhs_size);
    for (size_t i = 0; i < len; i++) {
        if (lhs->name[i] < rhs->name[i]) {
            return true;
        } else if (lhs->name[i] > rhs->name[i]) {
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

void CloneTreeNode::breadth_first_traversal(CloneTreeNode *root_node, vector<CloneTreeNode *> &ret, bool non_empty)
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
                ret.push_back(node);
                nodes.insert(node);
            } else {
                if (node->get_num_data() > 0) {
                    ret.push_back(node);
                    nodes.insert(node);
                }
            }
        }
        unordered_map<size_t, pair<double, CloneTreeNode *> > &children = node->get_idx2child();
        for (size_t i = 0; i < children.size(); i++)
        {
            CloneTreeNode *child = children[i].second;
            q.push(child);
        }
    }
}

void CloneTreeNode::get_cluster_labels(CloneTreeNode *root,
                                   const vector<BulkDatum *> &data,
                                   vector<unsigned int> &cluster_labels)
{
    vector<CloneTreeNode *> all_nodes;
    CloneTreeNode::breadth_first_traversal(root, all_nodes, false);
    unordered_map<const BulkDatum *, CloneTreeNode *> datum2node;
    construct_datum2node(all_nodes, datum2node);
    
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
    if (this->name.size() > other->name.size()) {
        return false;
    }
    for (size_t i = 0; i < name.size(); i++) {
        if (name.at(i) != other->name.at(i)) {
            return false;
        }
    }
    return true;
}

bool CloneTreeNode::IsDescendantOf(CloneTreeNode *other) {
    if (this->name.size() < other->name.size()) {
        return false;
    }
    for (size_t i = 0; i < other->name.size(); i++) {
        if (name.at(i) != other->name.at(i)) {
            return false;
        }
    }
    return true;
}

bool CloneTreeNode::IsCacheAllocated(size_t cell_count)
{
    return (sc_cache.size() == cell_count);
}

void CloneTreeNode::AllocateCache(size_t cell_count)
{
    if (sc_cache.size() == 0) {
        sc_cache.resize(cell_count);
    } else {
        if (cell_count != sc_cache.size()) {
            cerr << "Cache size does not match the required size of " << cell_count << ".\n";
            exit(-1);
        }
    }
}

void CloneTreeNode::UpdateCache(size_t c, double val)
{
    sc_cache[c] += val;
}

double CloneTreeNode::GetScCache(size_t c)
{
    return sc_cache[c];
}

void CloneTreeNode::construct_datum2node(vector<CloneTreeNode *> &all_nodes,
                                     unordered_map<const BulkDatum *, CloneTreeNode *> &datum2node)
{
    for (auto it = all_nodes.begin(); it != all_nodes.end(); ++it) {
        for (const BulkDatum *datum : (*it)->get_data())
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
            v = v->get_parent_node();
            unordered_set<const BulkDatum *> d = v->get_data();
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

double ScLikelihood(size_t loci_idx,
                    const BulkDatum *bulk,
                    const SingleCellData *sc,
                    bool has_snv,
                    const ModelParams &model_params) {
    double log_lik = 0.0;

    // if sc has mutation s, then there are 3 cases
    // 1. non-bursty
    // 2. bursty for variant
    // 3. bursty for reference
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
        log_lik_biallelic += log(1-locus.get_dropout_prob());
        double log_lik_dropout = log_beta_binomial_pdf(var_reads,
                                                       total_reads,
                                                       model_params.GetScBurstyDistributionAlpha0(),
                                                       model_params.GetScBurstyDistributionBeta0());
        log_lik_dropout += log(locus.get_dropout_prob());
        log_lik = log_add(log_lik_biallelic, log_lik_dropout);
//        log_lik_biallelic += log(model_params.GetScBiallelicProportion());
//        double log_lik_bursty_variant = log_beta_binomial_pdf(var_reads,
//                                                              total_reads,
//                                                              model_params.GetScBurstyDistributionAlpha0(),
//                                                              model_params.GetScBurstyDistributionBeta0());
//        log_lik_bursty_variant += log(model_params.GetScBurstyVariantProportion());
//        double log_lik_dropout = log_beta_binomial_pdf(var_reads,
//                                                      total_reads,
//                                                      model_params.GetScDropoutDistributionAlpha0(),
//                                                      model_params.GetScDropoutDistributionBeta0());
//        log_lik_dropout += log(model_params.GetScDropoutProportion());
//        log_lik = log_add(log_lik_dropout, log_lik_biallelic);
//        log_lik = log_add(log_lik, log_lik_bursty_variant);
    } else {
        log_lik = log_beta_binomial_pdf(var_reads,
                                        total_reads,
                                        model_params.GetSequencingError(),
                                        1 - model_params.GetSequencingError());
    }

    return log_lik;
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
    if (node->get_parent_node() == 0) {
        // this node is the root, represents the healthy population
        return log(gsl_ran_binomial_pdf(var_reads,
                                        seq_err,
                                        total_reads));
    }
    
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
    if (node->get_parent_node() == 0) {
        // this node is the root, represents the healthy population
        return log(gsl_ran_binomial_pdf(var_reads, seq_err, total_reads));
    }
    
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
    if (node->get_parent_node() == 0) {
        // This node is the root, represents the healthy population.
        return log(gsl_ran_binomial_pdf(var_reads,
                                        seq_err,
                                        total_reads));
    }

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

