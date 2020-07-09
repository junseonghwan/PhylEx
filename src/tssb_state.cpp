#include "tssb_state.hpp"
#include "utils.hpp"

TSSBState::TSSBState(const gsl_rng *random,
                     CloneTreeNode *root,
                     const ModelParams &params,
                     double (*log_lik_datum)(size_t regions,
                                             const CloneTreeNode *v,
                                             const BulkDatum *s,
                                             const ModelParams &params),
                     double (*log_lik_sc_at_site)(size_t loci_idx,
                                                  const BulkDatum *s,
                                                  const SingleCellData *c,
                                                  bool has_snv,
                                                  const ModelParams &params),
                     vector<BulkDatum *> *bulk_data,
                     vector<SingleCellData*> *sc_data) :
root(root),
bulk_data_(bulk_data),
has_sc_coverage_(bulk_data_->size(), false),
sc_data(sc_data),
log_lik_datum(log_lik_datum),
log_lik_sc_at_site(log_lik_sc_at_site)
{
//    if (sc_data != 0) {
//        sc_cache.resize(sc_data->size());
//    }
    
    // Pre-compute single cell likelihood.
    ProcessSingleCellData(params);
    
    root->sample_node_parameters(random, params, 0);
    for (size_t idx = 0; idx < bulk_data->size(); idx++) {
        this->initialize_data_assignment(random, idx, params);
    }
    
    // initialize the log likelihoods
    log_lik_bulk = compute_log_likelihood_bulk(params);
    cout << "Initializing cache..." << endl;
    log_lik_sc = compute_log_likelihood_sc_cached(params);
}

void TSSBState::ProcessSingleCellData(const ModelParams &model_params) {
    size_t cell_count = sc_data->size();
    size_t mutation_count = bulk_data_->size();
    sc_presence_matrix_ = EigenMatrix::Zero(cell_count, mutation_count);
    sc_absence_matrix_ = EigenMatrix::Zero(cell_count, mutation_count);
    // Evaluate single cell log likelihoods.
    size_t snv_sc_coverage_count = 0;
    for (size_t n = 0; n < mutation_count; n++) {
        size_t sc_coverage_count = 0;
        for (size_t c = 0; c < cell_count; c++) {
            sc_presence_matrix_(c,n) = log_lik_sc_at_site(n, bulk_data_->at(n), sc_data->at(c), true, model_params);
            sc_absence_matrix_(c,n) = log_lik_sc_at_site(n, bulk_data_->at(n), sc_data->at(c), false, model_params);
            auto total_reads = sc_data->at(c)->GetTotalReads(n);
            sc_coverage_count += (total_reads > 0) ? 1 : 0;
        }
        if (sc_coverage_count > 0) {
            has_sc_coverage_.push_back(true);
            snv_sc_coverage_count++;
        }
    }
    cout << "SNV with single cell coverage: " << snv_sc_coverage_count << "\n";
}

/*******************
 private functions
 *******************/
////void TSSBState::initialize_data_assignment(const gsl_rng *random, BulkDatum *datum, const ModelParams &params)
//{
//    // sample an initial tree
////    double u = uniform(random, 0.0, 1.0);
////    CloneTreeNode *node = CloneTreeNode::find_node(random, u, root, params); // this call may generate new nodes
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
//        //CloneTreeNode::find_node(random, u, root, params); // this call may generate new nodes
//    }
//
//    // assign the datum to the first child of root
//    assert(root->get_num_children() > 0);
//    CloneTreeNode *first_child_node = root->get_child(0).second;
//    assign_data_point(0, first_child_node, datum, params);
//}

void TSSBState::initialize_data_assignment(const gsl_rng *random, size_t mut_id, const ModelParams &params)
{
    if (root->get_num_children() == 0) {
        // Spawn one child.
        root->InitializeChild(random, params);
    }
    
    // assign the datum to the first child of root
    assert(root->get_num_children() > 0);
    CloneTreeNode *first_child_node = root->get_child(0).second;
    assign_data_point(0, first_child_node, mut_id, params, false);
}

void TSSBState::InitializeCacheForNode(CloneTreeNode *v)
{
    if (!v->IsCacheAllocated(sc_data->size())) {
        v->AllocateCache(sc_data->size());
        double log_lik;
        for (size_t c = 0; c < sc_data->size(); c++) {
            log_lik = compute_loglik_sc(v, c);
            v->UpdateCache(c, log_lik);
        }
    }
}

void TSSBState::update_sc_cache(CloneTreeNode *curr_node, CloneTreeNode *new_node, size_t mut_id, const ModelParams &params)
{
    //cout << "===== Update sc cache =====" << endl;
    if (sc_data == 0) {
        return;
    }
    
    vector<CloneTreeNode *> subtree_curr;
    vector<CloneTreeNode *> subtree_new;
    get_all_nodes(false, curr_node, subtree_curr);
    get_all_nodes(false, new_node, subtree_new);
    
    // Ensure cache is allocated for each node.
    for (CloneTreeNode *v : subtree_curr) {
        InitializeCacheForNode(v);
    }
    for (CloneTreeNode *v : subtree_new) {
        InitializeCacheForNode(v);
    }
    
    bool exp_mut_status;
    for (size_t c = 0; c < sc_data->size(); c++) {
        for (CloneTreeNode *v : subtree_curr) {
            double x = sc_presence_matrix_(c, mut_id);
            v->UpdateCache(c, -x);
            exp_mut_status = v->IsDescendantOf(new_node) ? true : false;
            double y = exp_mut_status ? sc_presence_matrix_(c, mut_id) : sc_absence_matrix_(c, mut_id);
            v->UpdateCache(c, y);
        }
        
        for (CloneTreeNode *v : subtree_new) {
            exp_mut_status = v->IsDescendantOf(curr_node) ? true : false;
            double x = exp_mut_status ? sc_presence_matrix_(c, mut_id) : sc_absence_matrix_(c, mut_id);
            v->UpdateCache(c, -x);
            double y = sc_presence_matrix_(c, mut_id);
            v->UpdateCache(c, y);
        }
    }
    
}

double TSSBState::LogLikDatum(CloneTreeNode *node,
                              BulkDatum *datum,
                              const ModelParams &model_params) {
    size_t region_count = node->get_node_parameter().GetRegionCount();
    double log_lik = 0.0;
    for (size_t i = 0; i < region_count; i++) {
        log_lik += log_lik_datum(i, node, datum, model_params);
    }
    return log_lik;
}

void TSSBState::slice_sample_data_assignment(const gsl_rng *random,
                                             size_t mut_id,
                                             const ModelParams &model_params)
{
    auto datum = bulk_data_->at(mut_id);
    double u_min = 0.0, u_max = 1.0;
    double u;
    CloneTreeNode *curr_node = datum2node[datum];
    CloneTreeNode *new_node = 0;

    double curr_log_lik_bulk = LogLikDatum(curr_node, datum, model_params);
    double log_slice = log(uniform(random, 0.0, 1.0)); // sample the slice
    double new_log_lik_bulk = 0.0;
    size_t iter = 0;
    
    while ((abs(u_max - u_min) > 1e-3)) {
        
        iter++;
        u = uniform(random, u_min, u_max);
        new_node = CloneTreeNode::find_node(random, u, root, model_params); // this call may generate new nodes
        
        if (new_node == curr_node) {
            break; // no need to carry out further computation
        }
        
        // assign data point from curr_node to new_node
        assign_data_point(curr_node, new_node, mut_id, model_params);
        
        // compute the log likelihood
        new_log_lik_bulk = LogLikDatum(new_node, datum, model_params);
        if (new_log_lik_bulk > (log_slice + curr_log_lik_bulk)) {
            // update log_lik_bulk, log_lik_sc, log_lik
            log_lik_bulk -= curr_log_lik_bulk;
            log_lik_bulk += new_log_lik_bulk;
            log_lik = log_lik_bulk + log_lik_sc;
            break;
        } else {
            // revert the changes
            assign_data_point(new_node, curr_node, mut_id, model_params);
            if (CloneTreeNode::less(new_node, curr_node)) {
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
    
    if (abs(u_max - u_min) <= 1e-3) {
    }
}

void TSSBState::slice_sample_data_assignment_with_sc(const gsl_rng *random,
                                                     size_t mut_id,
                                                     const ModelParams &model_params)
{
    auto datum = bulk_data_->at(mut_id);
    double u_min = 0.0, u_max = 1.0;
    double u;
    CloneTreeNode *curr_node = datum2node[datum];
    CloneTreeNode *new_node = 0;
    
    double curr_log_lik_bulk = LogLikDatum(curr_node, datum, model_params);
    double curr_log_lik_sc = compute_log_likelihood_sc_cached(model_params);
    double curr_log_lik = curr_log_lik_bulk + curr_log_lik_sc;
    
    double log_slice = log(uniform(random, 0.0, 1.0)); // sample the slice
    double new_log_lik_bulk = 0.0, new_log_lik_sc = 0.0, new_log_lik = 0.0;
    size_t iter = 0;
    
    while ((abs(u_max - u_min) > 1e-3)) {
        
        iter++;
        u = uniform(random, u_min, u_max);
        new_node = CloneTreeNode::find_node(random, u, root, model_params); // this call may generate new nodes
        
        if (new_node == curr_node) {
            break; // no need to carry out further computation
        }
        
        // assign data point from curr_node to new_node
        assign_data_point(curr_node, new_node, mut_id, model_params);
        
        // compute the log likelihood
        new_log_lik_bulk = LogLikDatum(new_node, datum, model_params);
        new_log_lik_sc = compute_log_likelihood_sc_cached(model_params);
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
            if (CloneTreeNode::less(new_node, curr_node)) {
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

const vector<BulkDatum *> &TSSBState::get_data() const
{
    return *bulk_data_;
}

////void TSSBState::assign_data_point(CloneTreeNode *curr_node, CloneTreeNode *new_node, BulkDatum *datum, const ModelParams &model_params)
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

void TSSBState::assign_data_point(CloneTreeNode *curr_node, CloneTreeNode *new_node, size_t mut_id, const ModelParams &model_params, bool update_cache)
{
    auto datum = bulk_data_->at(mut_id);
    if (curr_node != 0) {
        curr_node->remove_datum(datum);
    }
    
    // update the map : datum -> node
    datum2node[datum] = new_node;
    // add the datum to the new node
    new_node->add_datum(datum);
    
    // update single cell cache.
    // TODO: do it only if mut_id has single cell coverage.
    if (update_cache && sc_data != 0 && sc_data->size() > 0)
        update_sc_cache(curr_node, new_node, mut_id, model_params);
}


/**********
 public functions
 **********/
double TSSBState::get_log_prior_assignment(CloneTreeNode *root)
{
    // log prior of the data assignment requires computing the probability of
    // assignment to a node and the number of data points at that node.
    vector<pair<double, CloneTreeNode *> > mixture;
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

double TSSBState::get_log_lik()
{
    double log_prior_assignment = get_log_prior_assignment(get_root());
    log_lik = log_lik_bulk + log_lik_sc + log_prior_assignment;
    return log_lik;
}

double TSSBState::compute_log_likelihood_sc(bool verbose)
{
    if (sc_data->size() == 0) {
        return 0.0;
    }

    double log_lik_sc = 0.0;

    // Marginalize over assignment to the nodes.
    vector<CloneTreeNode *> all_nodes;
    get_all_nodes(true, root, all_nodes);
    
    double log_prior_assignment = -log(all_nodes.size());
    for (size_t c = 0; c < sc_data->size(); c++) {
        double log_lik_cell = DOUBLE_NEG_INF;
        for (size_t i = 0; i < all_nodes.size(); i++) {
            CloneTreeNode *node = all_nodes[i];
            double log_val = compute_loglik_sc(node, c);
            if (verbose) {
                cout << "[" << node->get_name() << "]: " << log_val << endl;
            }
            log_lik_cell = log_add(log_val, log_lik_cell);
        }
        log_lik_cell += log_prior_assignment;
        log_lik_sc += log_lik_cell;
    }
    
    return log_lik_sc;
}

double TSSBState::compute_log_likelihood_sc_cached(const ModelParams &params, bool verbose)
{
    if (sc_data == 0 || sc_data->size() == 0) {
        return 0.0;
    }

    // get all non-empty nodes, these are possible candidates for assignment of single cells
    vector<CloneTreeNode *> all_nodes;
    get_all_nodes(true, all_nodes);

    for (size_t i = 0; i < all_nodes.size(); i++) {
        CloneTreeNode *v = all_nodes[i];
        InitializeCacheForNode(v);
    }
    
    double log_prior_assignment = -log(all_nodes.size());
    
    double log_lik_sc = 0.0;
    double log_val;
    // TODO: this can be parallelized over cells
    for (size_t c = 0; c < sc_data->size(); c++) {
        // Marginalize over the nodes.
        double log_lik_cell = DOUBLE_NEG_INF;
        for (size_t i = 0; i < all_nodes.size(); i++) {
            CloneTreeNode *v = all_nodes[i];
            log_val = v->GetScCache(c);
            if (verbose) {
                cout << "[" << v->get_name() << "]: " << log_val << endl;
            }
            log_lik_cell = log_add(log_val, log_lik_cell);
        }
        log_val = (log_lik_cell + log_prior_assignment);
        //cout << "Cell " << c << ": " << log_val << endl;
        log_lik_sc += log_val;
    }

    //print_cache();
    return log_lik_sc;
}

double TSSBState::compute_loglik_sc(CloneTreeNode *v, size_t cell_id)
{
    double log_lik = 0.0;
    bool has_snv;
    auto loci_idxs = sc_data->at(cell_id)->GetLoci();
    for (size_t loci_idx : loci_idxs) {
        const BulkDatum *snv = bulk_data_->at(loci_idx);
        auto u = datum2node[snv];
        has_snv = u->IsAncestorOf(v);
        double log_val = has_snv ? sc_presence_matrix_(cell_id, loci_idx) :
        sc_absence_matrix_(cell_id, loci_idx);
        log_lik += log_val;
        //cout << "Has SNV: " << has_snv << ", " << log_val << endl;
    }
    return log_lik;
}

void TSSBState::resample_data_assignment(const gsl_rng *random, const ModelParams &params)
{
    for (size_t mut_id = 0; mut_id < bulk_data_->size(); mut_id++)
    {
        if (has_sc_coverage_[mut_id]) {
            slice_sample_data_assignment_with_sc(random, mut_id, params);
        } else {
            slice_sample_data_assignment(random, mut_id, params);
        }
    }
}

void TSSBState::move_datum(CloneTreeNode *new_node, size_t mut_id, const ModelParams &model_params)
{
    // should check if this node is part of the tree
    CloneTreeNode *curr_node = new_node;
    // if curr_node has the root as ancestor, it is a valid node
    while (curr_node != root) {
        if (curr_node == 0) {
            cerr << "Error: Node is not part of the tree!" << endl;
        }
        curr_node = curr_node->get_parent_node();
    }
    
    curr_node = datum2node[bulk_data_->at(mut_id)];
    assign_data_point(curr_node, new_node, mut_id, model_params);
}

// helper functions for updating the sticks
void TSSBState::reorder_sticks(const gsl_rng *random, const ModelParams &model_params)
{
    queue<CloneTreeNode *> q;
    q.push(root);
    while (!q.empty()) {
        CloneTreeNode *node = q.front();
        q.pop();
        
        node->reorder_sticks(random, model_params);
        // enqueue children nodes to the queue
        for (auto it = node->get_idx2child().begin(); it != node->get_idx2child().end(); ++it) {
            q.push(it->second.second);
        }
    }
}

void descend_and_update_names(CloneTreeNode *node)
{
    node->reset_children_names();
    CloneTreeNode *child = 0;
    
    unordered_map<size_t, pair<double, CloneTreeNode *> > &children = node->get_idx2child();
    for (auto it = children.begin(); it != children.end(); it++) {
        child = it->second.second;
        descend_and_update_names(child);
    }
}

size_t descend_and_cull(CloneTreeNode *node)
{
    size_t n_data = node->get_num_data();
    size_t n_data_desc = 0;
    CloneTreeNode *child = 0;
    
    unordered_map<size_t, pair<double, CloneTreeNode *> > &children = node->get_idx2child();
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
pair<size_t, size_t> TSSBState::descend_and_sample_sticks(const gsl_rng *random, CloneTreeNode *node, const ModelParams &params)
{
    size_t n_data = node->get_num_data();
    size_t total_num_data_at_desc = 0; // number of data points below this node (over all descendant nodes)
    
    unordered_map<size_t, pair<double, CloneTreeNode *> > &children = node->get_idx2child();
    vector<pair<size_t, size_t> > n_data_desc; // n_data_desc[i] stores <n_data at i, total_num_desc of i> (note: the second count does not include n_data at i)
    for (size_t i = 0; i < children.size(); i++)
    {
        CloneTreeNode *child = children.at(i).second;
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
size_t TSSBState::get_num_nodes()
{
    vector<CloneTreeNode *> nodes;
    get_all_nodes(nodes);
    return nodes.size();
}

// get the node where this SomaticSNV (ssm) first appears in
CloneTreeNode *TSSBState::get_node(const BulkDatum *datum)
{
    return datum2node[datum];
}

void TSSBState::get_all_nodes(bool non_empty, vector<CloneTreeNode *> &ret) const
{
    get_all_nodes(non_empty, root, ret);
}

void TSSBState::get_all_nodes(bool non_empty, CloneTreeNode *root_node, vector<CloneTreeNode *> &ret)
{
    // get nodes in the subtree rooted at root
    if (ret.size() > 0) {
        cout << "Warning: get_all_nodes() is going to clear vector ret." << endl;
        ret.clear();
    }
    if (root_node == 0)
        return;
    
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

// get the nodes along with their mass
// the nodes are ordered by depth first traversal
void TSSBState::get_mixture(CloneTreeNode *node, double mass, vector<pair<double, CloneTreeNode *> > &ret)
{
    double node_mass = mass * node->get_nu_stick();
    ret.push_back(make_pair(node_mass, node));
    unordered_map<size_t, pair<double, CloneTreeNode *> > &children = node->get_idx2child();
    double cum_prod = 1.0;
    vector<double> weights;
    for (size_t i = 0; i < children.size(); i++)
    {
        CloneTreeNode *child = children[i].second;
        double val = cum_prod * (1 - children[i].first);
        weights.push_back(1 - val);
        cum_prod *= (1 - children[i].first);
        get_mixture(child, (1 - node->get_nu_stick()) * mass * weights[i], ret);
    }
}

//void TSSBState::insert_datum(const gsl_rng *random, BulkDatum * datum, const ModelParams &params)
//{
//    double u = gsl_ran_flat(random, 0, 1);
//    CloneTreeNode *node = CloneTreeNode::find_node(random, u, root, params);
//    bulk_data->push_back(datum);
//    assign_data_point(0, node, datum, params);
//    
//    // do not update sc_cache as it will initialized and called in compute_likelihood_sc_cached
//}

void TSSBState::initialize_sc_cache(const ModelParams &model_params)
{
    //sc_cache.resize(sc_data->size());
    vector<CloneTreeNode *> all_nodes;
    get_all_nodes(false, all_nodes);
    unordered_map<CloneTreeNode *, unordered_set<const BulkDatum *> > node2snvs;
    for (CloneTreeNode *v : all_nodes)
    {
        unordered_set<const BulkDatum *> snvs;
        CloneTreeNode::GetDataset(v, snvs);
        node2snvs[v] = snvs;
    }
    for (auto v : all_nodes) {
        // Check if cache is initialized.
        v->AllocateCache(sc_data->size());
        for (size_t c = 0; c < sc_data->size(); c++) {
            v->UpdateCache(c, compute_loglik_sc(v, c));
        }
    }
}

void TSSBState::set_sc_data(vector<SingleCellData *> *sc_data, const ModelParams &model_params)
{
    this->sc_data = sc_data;
    initialize_sc_cache(model_params);
}

// get all nodes in this tree with at least one datum
// via breadth-first traversal
// populate ret
void TSSBState::get_all_nodes(vector<CloneTreeNode *> &ret) const
{
    get_all_nodes(false, root, ret);
}

gsl_matrix *TSSBState::get_ancestral_matrix(TSSBState &state)
{
    size_t N = state.bulk_data_->size();
    vector<unordered_set<const BulkDatum *> > ancestors(N);
    for (size_t i = 0; i < N; i++) {
        const BulkDatum *datum = state.bulk_data_->at(i);
        CloneTreeNode *v = state.get_node(datum);
        // trace up to the root and set the row of A
        while (v != state.root) {
            v = v->get_parent_node();
            unordered_set<const BulkDatum *> data = v->get_data();
            ancestors[i].insert(data.begin(), data.end());
        }
    }

    // A_ij = 1 if i occurs before j
    gsl_matrix *A = gsl_matrix_alloc(N, N);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            if (ancestors[j].count(state.bulk_data_->at(i))) {
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

double TSSBState::compute_log_likelihood_bulk(size_t region,
                                              const ModelParams &model_params)
{
    // compute the log likelihood over data S
    double log_lik_bulk_region = 0.0;
    for (size_t i = 0; i < bulk_data_->size(); i++) {
        CloneTreeNode *node = datum2node[(*bulk_data_)[i]];
        double log_val = log_lik_datum(region, node, (*bulk_data_)[i], model_params);
        log_lik_bulk_region += log_val;
    }
    return log_lik_bulk_region;
}

// compute (or get) the log likelihood of the current state
double TSSBState::compute_log_likelihood_bulk(const ModelParams &model_params)
{
    size_t region_count = root->get_node_parameter().GetRegionCount();
    double log_lik_bulk = 0.0;
    for (size_t i = 0; i < region_count; i++) {
        log_lik_bulk += compute_log_likelihood_bulk(i, model_params);
    }
    return log_lik_bulk;
}

// update stick lengths, and exchange stick ordering
void TSSBState::update_sticks(const gsl_rng *random, const ModelParams &params)
{
    descend_and_sample_sticks(random, root, params);
    reorder_sticks(random, params);
    descend_and_sample_sticks(random, root, params);
}


double log_prod_beta(vector<CloneTreeNode *> &nodes, const ModelParams &params)
{
    double log_sum = 0.0;
    double alpha;
    for (CloneTreeNode *node : nodes) {
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

void TSSBState::sample_alpha0(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<CloneTreeNode *> &nodes)
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

void TSSBState::sample_lambda(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<CloneTreeNode *> &nodes)
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

void TSSBState::sample_gamma(const gsl_rng *random, size_t n_mh_iter, ModelParams &params, vector<double> psi_sticks)
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

void TSSBState::update_hyper_params(const gsl_rng *random,
                                           size_t n_mh_iter,
                                           TSSBState &state,
                                           ModelParams &params)
{
    vector<CloneTreeNode *> all_nodes;
    state.get_all_nodes(all_nodes);
    vector<double> psi_sticks;
    for (CloneTreeNode *node : all_nodes) {
        for (size_t i = 0; i < node->get_idx2child().size(); i++) {
            psi_sticks.push_back(node->get_idx2child()[i].first);
        }
    }
    
    sample_alpha0(random, n_mh_iter, params, all_nodes);
    sample_lambda(random, n_mh_iter, params, all_nodes);
    sample_gamma(random, n_mh_iter, params, psi_sticks);
}

string TSSBState::print()
{
    string ret = "";
    queue<CloneTreeNode *> q;
    q.push(root);
    while (!q.empty()) {
        CloneTreeNode *node = q.front();
        q.pop();
        
        ret += node->print() + "\n";
        // enqueue children nodes to the queue
        unordered_map<size_t, pair<double, CloneTreeNode *> > idx2child = node->get_idx2child();
        for (size_t i = 0; i < idx2child.size(); i++) {
            q.push(idx2child[i].second);
        }
    }
    
    ret += "\n";
    // print data assignments
    //double log_lik = 0.0, log_lik_datum = 0.0;
    //    for (size_t idx = 0; idx < bulk_data->size(); idx++) {
    //        CloneTreeNode *node = datum2node[(*bulk_data)[idx]];
    //        ret += "Data point " + to_string(idx) + " assigned to <" + (node->is_root() ? "root" : node->get_name()) + ">\n";
    //BulkDatum *datum = (*data)[idx];
    //        log_lik_datum = node->compute_log_likelihood_of_datum(*datum, params);
    //        ret += "Log likelihood: " + to_string(log_lik_datum) + "\n";
    //        log_lik += log_lik_datum;
    //    }
    //    ret += "Complete data log lik: " + to_string(log_lik) + "\n";
    return ret;
}

//void TSSBState::clear_cache()
//{
//    if (sc_data->size() == 0)
//        return;
//    
//    vector<CloneTreeNode *> active_nodes_vec;
//    get_all_nodes(false, root, active_nodes_vec);
//    unordered_set<CloneTreeNode *> active_nodes(active_nodes_vec.begin(), active_nodes_vec.end());
//    
//    unordered_set<CloneTreeNode *> purge_set;
//}

double update_cellular_prev_recursive(size_t region, CloneTreeNode * node)
{
    unordered_map<size_t, pair<double, CloneTreeNode *> > &children = node->get_idx2child();
    const CloneTreeNodeParam &param = node->get_node_parameter();
    double children_prev = 0.0;
    for (auto it = children.begin(); it != children.end(); ++it)
    {
        CloneTreeNode *child = it->second.second;
        children_prev += update_cellular_prev_recursive(region, child);
    }
    auto cell_prev = param.get_clone_freqs(region) + children_prev;
    node->set_cellular_prev(region, cell_prev);
    return cell_prev;
}

void update_params(size_t region,
                   vector<CloneTreeNode *> &nodes,
                   double *new_clone_freq)
{
    // update clone frequencies
    for (size_t i = 0; i < nodes.size(); i++) {
        nodes[i]->set_clone_freq(region, new_clone_freq[i]);
    }
    // update cellular prevalences
    update_cellular_prev_recursive(region, nodes[0]);
}

void get_clone_freqs(size_t region,
                     TSSBState &state,
                     double *clone_freqs)
{
    vector<CloneTreeNode *> nodes;
    state.get_all_nodes(false, nodes);
    CloneTreeNode *node = 0;
    for (size_t i = 0; i < nodes.size(); i++)
    {
        node = nodes[i];
        clone_freqs[i] = node->get_node_parameter().get_cellular_prevs(region);
        unordered_map<size_t, pair<double, CloneTreeNode *> > &children = node->get_idx2child();
        for (auto it = children.begin(); it != children.end(); ++it)
        {
            CloneTreeNode *child = it->second.second;
            clone_freqs[i] -= child->get_node_parameter().get_cellular_prevs(region);
        }
    }
}

bool check_clone_freq(size_t region, CloneTreeNode *root)
{
    double sum_of_freqs = 0.0;
    queue<CloneTreeNode *> q;
    q.push(root);
    CloneTreeNode *node = 0;
    while (!q.empty()) {
        node = q.front();
        q.pop();
        sum_of_freqs += node->GetCloneFreqs(region);
        unordered_map<size_t, pair<double, CloneTreeNode *> > &idx2child = node->get_idx2child();
        for (size_t i = 0; i < idx2child.size(); i++) {
            q.push(idx2child[i].second);
        }
    }
    return(abs(sum_of_freqs - 1.0) < 1e-3);
}

// This function is for sampling the parameters for one region.
double sample_params_dirichlet(size_t region,
                               const gsl_rng *random,
                               size_t n_mh_iter,
                               TSSBState &tree,
                               const ModelParams &params)
{
    double curr_log_lik = tree.compute_log_likelihood_bulk(region, params);
    double new_log_lik;
    
    // Perform MH to update the parameters
    // collect current parameters via tree traversal
    // sample from DP centered at current parameters multipled by constant factor
    double K = params.get_dir_conc_mult_factor();
    vector<CloneTreeNode *> nodes;
    tree.get_all_nodes(nodes);
    size_t dim = nodes.size();
    double *curr_clone_freq = new double[dim];
    double *new_clone_freq = new double[dim];
    double *dir_concentration_params_curr = new double[dim];
    double *dir_concentration_params_new = new double[dim];
    
    get_clone_freqs(region, tree, curr_clone_freq);
    
    size_t n_accepts = 0;
    double log_proposal_new = 0.0, log_proposal_curr = 0.0, log_accept, log_unif;
    for (size_t i = 0; i < n_mh_iter; i++) {
        for (size_t i = 0; i < nodes.size(); i++) {
            dir_concentration_params_curr[i] = K * curr_clone_freq[i] + 1;
        }
        gsl_ran_dirichlet(random, nodes.size(), dir_concentration_params_curr, new_clone_freq);
        for (size_t i = 0; i < nodes.size(); i++) {
            dir_concentration_params_new[i] = K * new_clone_freq[i] + 1;
        }
        update_params(region, nodes, new_clone_freq);
        
        new_log_lik = tree.compute_log_likelihood_bulk(region, params);
        log_proposal_new = log_dirichlet_pdf(dim, dir_concentration_params_curr, new_clone_freq);
        log_proposal_curr = log_dirichlet_pdf(dim, dir_concentration_params_new, curr_clone_freq);
        log_unif = log(uniform(random, 0, 1));
        log_accept = (new_log_lik - curr_log_lik) + (log_proposal_curr - log_proposal_new);
        if (log_unif < log_accept) {
            copy(new_clone_freq, new_clone_freq+dim, curr_clone_freq);
            copy(dir_concentration_params_new, dir_concentration_params_new + dim, dir_concentration_params_curr);
            curr_log_lik = new_log_lik;
            n_accepts++;
        } else {
            update_params(region, nodes, curr_clone_freq);
        }
    }
    cout << "n accepts: " << n_accepts << endl;
    double r = (double) n_accepts / n_mh_iter;
    
    delete [] curr_clone_freq;
    delete [] new_clone_freq;
    delete [] dir_concentration_params_curr;
    delete [] dir_concentration_params_new;
    
    return r;
}

double sample_params_dirichlet(const gsl_rng *random,
                               size_t n_mh_iter,
                               TSSBState &tree,
                               const ModelParams &params)
{
    size_t region_count = tree.get_root()->get_node_parameter().GetRegionCount();
    double ar = 0.0;
    for (size_t i = 0; i < region_count; i++) {
        ar += sample_params_dirichlet(i, random, n_mh_iter, tree, params);
    }
    return ar/region_count;
}

void cull(CloneTreeNode *root)
{
    descend_and_cull(root);
    descend_and_update_names(root);
    
    vector<CloneTreeNode *> nodes;
    CloneTreeNode::breadth_first_traversal(root, nodes);
    
    // update the clone frequencies
    size_t region_count = root->get_node_parameter().GetRegionCount();
    CloneTreeNode *curr_node;
    for (size_t region = 0; region < region_count; region++) {
        for (size_t i = 0; i < nodes.size(); i++) {
            curr_node = (CloneTreeNode *)nodes[i];
            unordered_map<size_t, pair<double, CloneTreeNode *> > &idx2child = curr_node->get_idx2child();
            double children_cellular_prev = 0.0;
            for (size_t idx = 0; idx < idx2child.size(); idx++) {
                CloneTreeNode *child_node = (CloneTreeNode *)idx2child[idx].second;
                children_cellular_prev += child_node->GetCellularPrevs(region);
            }
            // cull does not alter cellular prevalence but it does change the clone frequncy because some children gets removed
            curr_node->set_clone_freq(region, curr_node->GetCellularPrevs(region) - children_cellular_prev);
        }
    }
}
