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

void CloneTreeNodeParam::set_clone_freq(double val)
{
    this->clone_freq = val;
}
void CloneTreeNodeParam::set_cellular_prev(double val)
{
    this->cellular_prev = val;
}

bool CloneTreeNodeParam::is_consistent()
{
    return (clone_freq <= cellular_prev);
}

CloneTreeNode::CloneTreeNode(size_t child_idx,
                             Node<BulkDatum,CloneTreeNodeParam> *parent) :
Node<BulkDatum,CloneTreeNodeParam>(child_idx, parent)
{
}

//CloneTreeNode::CloneTreeNode(size_t child_idx,
//                             //Node<BulkDatum,CloneTreeNodeParam> *parent,
//                             CloneTreeNode *parent,
//                             TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> *tssb) :
//Node<BulkDatum,CloneTreeNodeParam>(child_idx, parent, tssb)
//{
////    this->ancestral_data.insert(parent->ancestral_data.begin(), parent->ancestral_data.end());
////    this->ancestral_data.insert(parent->data.begin(), parent->data.end());
//}

const CloneTreeNodeParam &CloneTreeNode::get_node_parameter() const
{
    return param;
}

//void CloneTreeNode::sample_node_parameters_prior(const gsl_rng *random, const ModelParams &params, Node<BulkDatum,CloneTreeNodeParam> *parent)
//{
//    if (parent_node == 0) {
//        this->param.set_cellular_prev(1.0);
//        this->param.set_clone_freq(1.0);
//    } else {
//        CloneTreeNode *parent_node = (CloneTreeNode *)parent;
//        double parent_clone_freq = parent_node->get_node_parameter().get_clone_freq();
//        double curr_cellular_prev = this->get_node_parameter().get_cellular_prev();
//
//        double min = 0.0;
//        double max = 1.0;
//        if (curr_cellular_prev == 0.0) { // new node
//            max = parent_clone_freq;
//            // sample from the prior distribution
//            curr_cellular_prev = uniform(random) * parent_clone_freq;
//            this->param.set_cellular_prev(curr_cellular_prev);
//            this->param.set_clone_freq(curr_cellular_prev);
//        } else {
//            // find minimum and maximum values
//            for (size_t i = 0; i < get_idx2child().size(); i++) {
//                min += get_idx2child()[i].second->get_node_parameter().get_cellular_prev();
//            }
//            max = curr_cellular_prev + parent_clone_freq;
//
//            double new_cellular_prev = gsl_ran_flat(random, min, max);
//            this->param.set_cellular_prev(new_cellular_prev);
//            this->param.set_clone_freq(new_cellular_prev - min);
//        }
//
//        parent_node->param.set_clone_freq(max - param.get_cellular_prev());
//    }
//}

void CloneTreeNode::sample_node_parameters(const gsl_rng *random, const ModelParams &params,  Node<BulkDatum,CloneTreeNodeParam> *parent)
{
    if (parent_node == 0) {
        this->param.set_cellular_prev(1.0);
        this->param.set_clone_freq(1.0);
    } else {
        CloneTreeNode *parent_node = (CloneTreeNode *)parent;
        double parent_clone_freq = parent_node->get_node_parameter().get_clone_freq();
        double curr_cellular_prev = uniform(random) * parent_clone_freq;
        this->param.set_cellular_prev(curr_cellular_prev);
        this->param.set_clone_freq(curr_cellular_prev);
        parent_node->param.set_clone_freq(parent_clone_freq - param.get_cellular_prev());
//        double curr_cellular_prev = this->get_node_parameter().get_cellular_prev();
//
//        double min = 0.0;
//        double max = 1.0;
//        if (curr_cellular_prev == 0.0) { // new node
//            max = parent_clone_freq;
//            // sample from the prior distribution
//            curr_cellular_prev = uniform(random) * parent_clone_freq;
//            this->param.set_cellular_prev(curr_cellular_prev);
//            this->param.set_clone_freq(curr_cellular_prev);
//        } else {
//            // find minimum and maximum values
//            for (size_t i = 0; i < get_idx2child().size(); i++) {
//                min += get_idx2child()[i].second->get_node_parameter().get_cellular_prev();
//            }
//            max = curr_cellular_prev + parent_clone_freq;
//
//            double new_cellular_prev = gsl_ran_flat(random, min, max);
//            this->param.set_cellular_prev(new_cellular_prev);
//            this->param.set_clone_freq(new_cellular_prev - min);
//        }
        
//        parent_node->param.set_clone_freq(max - param.get_cellular_prev());
    }

//    if (parent_node == 0) {
//        this->param.set_cellular_prev(1.0);
//        this->param.set_clone_freq(1.0);
//    } else {
//        CloneTreeNode *parent_node = (CloneTreeNode *)parent;
//        double parent_clone_freq = parent_node->get_node_parameter().get_clone_freq();
//        double curr_cellular_prev = this->get_node_parameter().get_cellular_prev();
//        double min = 0.0;
//        double max = 1.0;
//        if (curr_cellular_prev == 0.0) { // new node
//            max = parent_clone_freq;
//            // sample from the prior distribution
//            curr_cellular_prev = uniform(random) * parent_clone_freq;
//            this->param.set_cellular_prev(curr_cellular_prev);
//            this->param.set_clone_freq(curr_cellular_prev);
//        } else {
//            // find minimum and maximum values
//            for (size_t i = 0; i < get_idx2child().size(); i++) {
//                min += get_idx2child()[i].second->get_node_parameter().get_cellular_prev();
//            }
//            max = curr_cellular_prev + parent_clone_freq;
//        }
//
//        double curr_log_lik = compute_log_likelihood(params);
//        // independent MH
//        for (size_t i = 0; i < 50; i++) {
//            double new_cellular_prev = gsl_ran_flat(random, min, max);
//            this->param.set_cellular_prev(new_cellular_prev);
//            this->param.set_clone_freq(new_cellular_prev - min);
//            // compute log likelihood of data assigned to this node
//            double new_log_lik = compute_log_likelihood(params);
//            double log_unif = log(gsl_ran_flat(random, 0, 1));
//            if (log_unif < (new_log_lik - curr_log_lik)) {
//                curr_log_lik = new_log_lik;
//                curr_cellular_prev = new_cellular_prev;
//            } else {
//                this->param.set_clone_freq(curr_cellular_prev - min);
//                this->param.set_cellular_prev(curr_cellular_prev);
//            }
//            parent_node->param.set_clone_freq(max - param.get_cellular_prev());
//        }
//    }
}

void CloneTreeNode::get_snvs(Node<BulkDatum,CloneTreeNodeParam> *node, unordered_set<Locus> &ret)
{
    // traverse to the root to retrieve all SNVs
    Node<BulkDatum,CloneTreeNodeParam> *v = node;
    while (v != 0) {
        unordered_set<BulkDatum *> data = v->get_data();
        for (BulkDatum *datum : data) {
            ret.insert(datum->GetLocus());
        }
        v = v->get_parent_node();
    }
}

double get_branch_length(Node<BulkDatum,CloneTreeNodeParam> *child, Node<BulkDatum,CloneTreeNodeParam> *parent)
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


//double CloneTreeNode::compute_expectation_of_xi(const ModelParams &model_params)
//{
//    double xi = 0.0;
//
//    // perform dynamic programming on the tree
//    vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
//    tssb->get_all_nodes(false, this, nodes); // breadth first retrieval
//
//    // DP table:
//    // for each node, store 2D-array of dimension M x M
//    size_t M = model_params.get_max_cn();
//    unordered_map<Node<BulkDatum,CloneTreeNodeParam> *, vector<vector<double> > > dp_table;
//    for (int i = nodes.size() - 1; i >= 0; i--) {
//        auto node = nodes[i];
//        if (node->is_leaf()) {
//            vector<vector<double> > empty(M);
//            for (size_t chi_r = 0; chi_r < M; chi_r++) {
//                for (size_t chi_v = 0; chi_v < M; chi_v++) {
//                    empty[chi_r].push_back(0);
//                }
//            }
//            dp_table[node] = empty;
//        } else {
//            // get children
//            unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &idx2child = node->get_idx2child();
//            Node<BulkDatum,CloneTreeNodeParam> *child;
//
//            // fill the DP table
//            vector<vector<double> > table(M);
//            for (size_t chi_r = 0; chi_r < M; chi_r++) {
//                for (size_t chi_v = 0; chi_v < M; chi_v++) {
//                    double sum = 0.0;
//                    for (size_t j = 0; j < node->get_num_children(); j++) {
//                        child = idx2child.at(j).second;
//                        if (dp_table.count(child) == 0) // if child does not have any data, nothing to compute so continue
//                            continue;
//
//                        double eta = child->get_node_parameter().get_clone_freq();
//                        // retrieve the DP table for the child node
//                        vector<vector<double> > &child_table = dp_table[child];
//
//                        for (size_t chi_v_child = 0; chi_v_child < M; chi_v_child++) {
//                            for (size_t chi_r_child = 0; chi_r_child < M; chi_r_child++) {
//                                size_t chi_t_child = chi_v_child + chi_r_child;
//                                double transition_prob_r = model_params.compute_transition_prob(chi_r, chi_r_child);
//                                double transition_prob_v = model_params.compute_transition_prob(chi_v, chi_v_child);
//                                double transition_prob = transition_prob_r * transition_prob_v;
//                                double g = 0.0;
//                                if (chi_v_child == 0) {
//                                    g = model_params.get_seq_error();
//                                } else {
//                                    g = (double)chi_v_child/chi_t_child;
//                                }
//                                double val = transition_prob * (eta * g + child_table[chi_r_child][chi_v_child]);
//                                
//                                //printf("chi_r=%lu, chi_v=%lu, g=%f, p_v=%f, p_r=%f, val=%f\n", chi_r_child, chi_v_child, g, transition_prob_v, transition_prob_r, val);
//                                sum += val;
//                            }
//                        }
//                    }
//                    table[chi_r].push_back(sum);
//                }
//            }
//            dp_table[node] = table;
//        }
//    }
//
//    // sum over the copy number at the root
//    // compute the initial distribution
//    double branch_length = get_branch_length(this, tssb->get_root());
//    gsl_matrix *trans_mat = gsl_matrix_alloc(model_params.get_max_cn(), model_params.get_max_cn());
//    model_params.compute_transition_matrix(branch_length - 1, trans_mat);
//    vector<vector<double> > &table = dp_table[this];
//    double eta = this->get_node_parameter().get_clone_freq();
//    for (size_t chi_v = 0; chi_v < M; chi_v++) {
//        for (size_t chi_r = 0; chi_r < M; chi_r++) {
//            size_t chi_t = chi_r + chi_v;
//            double initial_prob = model_params.compute_initial_prob(trans_mat, chi_r, chi_v, 2);
//            double g = 0.0;
//            if (chi_v == 0) {
//                g = model_params.get_seq_error();
//            } else {
//                g = (double)chi_v/chi_t;
//            }
//            double val = initial_prob * (eta * g + table[chi_r][chi_v]);
//            //printf("chi_r=%lu, chi_v=%lu, g=%f, p=%f, val=%f\n", chi_r, chi_v, g, initial_prob, val);
//            xi += val;
//        }
//    }
//
//    xi += (1 - get_node_parameter().get_cellular_prev()) * model_params.get_seq_error();
//
//    gsl_matrix_free(trans_mat);
//    return xi;
//}

//double CloneTreeNode::compute_log_likelihood_of_datum(BulkDatum &datum, const ModelParams &model_params)
//{
//    if (datum.get_n_reads() == 0) {
//        if (datum.get_n_variant_reads() > 0) {
//            cerr << "Error in the data. Num reads: " << datum.get_n_reads() << ". Num variants: " << datum.get_n_variant_reads() << endl;
//            exit(-1);
//        } else {
//            return 0.0;
//        }
//    }
//
//    double seq_err = model_params.get_seq_error();
//    if (get_parent_node() == 0) {
//        // this node is the root, represents the healthy population
//        return log(gsl_ran_binomial_pdf(datum.get_n_variant_reads(), seq_err, datum.get_n_reads()));
//    }
//    double phi = param.get_cellular_prev();
//    double xi = 0.0;
//
//    if (datum.is_var_cn_observed()) {
//        xi += (1 - phi) * seq_err;
//        xi += phi * datum.get_var_cn()/datum.get_total_cn();
//        if (xi > 1) {
//            cerr << "Error: probability of success: " << xi <<  " > 1." << endl;
//            xi = 1.0;
//        } else if (xi < 0) {
//            cerr << "Error: negative probability of success." << endl;
//            exit(-1);
//        }
//    } else {
//        xi = compute_expectation_of_xi(model_params);
//        assert(xi < 1+1e-6);
//        assert(xi > 0);
//    }
//    double log_val = log(gsl_ran_binomial_pdf(datum.get_n_variant_reads(), xi, datum.get_n_reads()));
//    return log_val;
//}

// compute log p(a_{c,s} | \zeta_c = this)
//double CloneTreeNode::compute_log_likelihood_of_sc_at_site(const BulkDatum *s, const SingleCellData *sc, size_t exp_mut_status, const ModelParams &hyper_params)
//{
//    const unordered_map<Locus, size_t> &sc_mut_map = sc->get_mutation_map();
//    size_t obs_mut_status = sc_mut_map.at(s->get_loci());
//    double log_lik = 0.0;
//    if (obs_mut_status == 0) {
//        // mutation is not observed -- either a false negative or a true negative
//        double beta = hyper_params.get_sc_fn_rate();
//        log_lik += exp_mut_status * log(beta) + (1 - exp_mut_status) * log(1 - beta);
//    } else if (obs_mut_status == 1) {
//        // mutation is harboured -- either a false positive or a true positive
//        double alpha = hyper_params.get_sc_fp_rate();
//        log_lik += exp_mut_status * log(1-alpha) + (1 - exp_mut_status) * log(alpha);
//    }
//    return log_lik;
//}
//
//double CloneTreeNode::compute_log_likelihood_of_sc(const SingleCellData *sc_datum, const unordered_set<BulkDatum *> &snas, const ModelParams &hyper_params)
//{
//    double log_lik = 0.0;
//
//    // extract Loci from snas
//    unordered_set<Locus> muts_harboured;
//    for (BulkDatum *datum : snas) {
//        muts_harboured.insert(datum->get_loci());
//    }
//
//    const unordered_map<Locus, size_t> &sc_mut_map = sc_datum->get_mutation_map();
//    for (auto it = sc_mut_map.begin(); it != sc_mut_map.end(); ++it) {
//        const Locus &loci = it->first;
//        size_t observed_mut_status = it->second;
//
//        // check if mutation is harboured at the node
//        size_t expected_mut_status = (muts_harboured.count(loci) > 0) ? 1 : 0;
//
//        if (observed_mut_status == 0) {
//            // mutation is not harboured
//            double beta = hyper_params.get_sc_fn_rate();
//            log_lik += expected_mut_status * log(beta) + (1 - expected_mut_status) * log(1 - beta);
//        } else if (observed_mut_status == 1) {
//            // mutation is harboured
//            double alpha = hyper_params.get_sc_fp_rate();
//            log_lik += expected_mut_status * log(1-alpha) + (1 - expected_mut_status) * log(alpha);
//        } else {
//            // N/A: no reads
//            // does not contribute to the likelihood
//        }
//    }
//    
//    if (log_lik > 0) {
//        cout << "here" << endl;
//    }
//
//    return log_lik;
//}

//double CloneTreeNode::compute_log_likelihood(const ModelParams &hyper_params)
//{
//    double log_lik = 0.0;
//    for (auto it = get_data().begin(); it != get_data().end(); ++it) {
//        log_lik += compute_log_likelihood_of_datum(*(*it), hyper_params);
//    }
//    return log_lik;
//}

Node<BulkDatum,CloneTreeNodeParam> *CloneTreeNode::spawn_child(double psi)
{
    size_t j = get_idx2child().size();
    //Node<BulkDatum,CloneTreeNodeParam> *child = new CloneTreeNode(j, this, tssb);
    CloneTreeNode *child = new CloneTreeNode(j, this);
    idx2child[j] = make_pair(psi, child);
    return child;
}

string CloneTreeNode::print()
{
    string psi_stick_str = "(";
    unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &children = get_idx2child();
    for (size_t i = 0; i < get_num_children(); i++) {
        psi_stick_str += to_string(children[i].first);
        if (i < get_num_children() - 1) {
            psi_stick_str += ", ";
        }
    }
    psi_stick_str += ")";
    // print the Node<S,P> name, nu-stick, Node<S,P> params
    string ret = "[" + (parent_node == 0 ? "root" : this->get_name()) + ", ";
    ret += "phi=" + to_string(get_node_parameter().get_cellular_prev()) + ", ";
    ret += "eta=" + to_string(get_node_parameter().get_clone_freq()) + ", ";
    ret += "nu=" + to_string(nu) + ", ";
    ret += "psi=" + psi_stick_str + ", ";
    
    ret += "data=( ";
    for (BulkDatum *datum : get_data()) {
        ret += datum->GetId() + " ";
    }
    ret += ")";
    ret += "]";

    return ret;
}

CloneTreeNode *CloneTreeNode::create_root_node()
{
    auto node = new CloneTreeNode(0, 0);
    return node;
}

void CloneTreeNode::change_clone_freq(double change)
{
    double new_val = param.get_clone_freq() + change;
    param.set_clone_freq(new_val);
    
    new_val = param.get_cellular_prev() + change;
    param.set_cellular_prev(new_val);
}

void CloneTreeNode::set_clone_freq(double new_val)
{
    param.set_clone_freq(new_val);
}

void CloneTreeNode::set_cellular_prev(double new_val)
{
    param.set_cellular_prev(new_val);
}

double ScLikelihood(const BulkDatum *s, const SingleCellData *sc, bool has_snv, const ModelParams &model_params) {
    double log_lik = 0.0;

    cout << s->GetId() << endl;
    // if sc has mutation s, then there are 3 cases
    // 1. non-bursty
    // 2. bursty for variant
    // 3. bursty for reference
    const Locus &somatic_locus = s->GetLocus();
    const LocusDatum *locus_datum = sc->get_locus_datum(somatic_locus);
    if (locus_datum == 0 || locus_datum->get_n_total_reads() == 0) {
        return 0.0;
    }
    size_t total_reads = locus_datum->get_n_total_reads();
    if (has_snv) {
        double alpha = somatic_locus.get_alpha();
        double beta = somatic_locus.get_beta();
        double log_lik_non_bursty = log_beta_binomial_pdf(locus_datum->get_n_var_reads(),
                                                          total_reads,
                                                          alpha,
                                                          beta);
        log_lik_non_bursty += log(1 - somatic_locus.get_dropout_prob());
        double log_lik_bursty = log_beta_binomial_pdf(locus_datum->get_n_var_reads(),
                                                      total_reads,
                                                      model_params.GetScDropoutDistributionAlpha0(),
                                                      model_params.GetScDropoutDistributionBeta0());
        log_lik_bursty += log(somatic_locus.get_dropout_prob());
        log_lik = log_add(log_lik_bursty, log_lik_non_bursty);
    } else {
        log_lik = log_beta_binomial_pdf(locus_datum->get_n_var_reads(),
                                        total_reads,
                                        model_params.get_seq_error(),
                                        1 - model_params.get_seq_error());
    }

    return log_lik;
}

double BulkLogLikWithTotalCopyNumber(const Node<BulkDatum, CloneTreeNodeParam> *node,
                                     const BulkDatum *datum,
                                     const ModelParams &model_params)
{
    if (datum->GetReadCount() == 0) {
        if (datum->GetVariantReadCount() > 0) {
            cerr << "Error in the data.\n";
            cerr << "Num reads: " << datum->GetReadCount() << ".\n";
            cerr << "Num variants: " << datum->GetVariantReadCount() << ".\n";
            exit(-1);
        } else {
            return 0.0;
        }
    }

    double seq_err = model_params.get_seq_error();
    if (node->get_parent_node() == 0) {
        // this node is the root, represents the healthy population
        return log(gsl_ran_binomial_pdf(datum->GetVariantReadCount(),
                                        seq_err,
                                        datum->GetReadCount()));
    }
    
    const CloneTreeNodeParam &param = node->get_node_parameter();
    double phi = param.get_cellular_prev();
    size_t total_cn = datum->GetTotalCN();
    double xi = 0.0;

    if (total_cn == 0) {
        xi = seq_err;
        return log(gsl_ran_binomial_pdf(datum->GetVariantReadCount(),
                                        xi,
                                        datum->GetReadCount()));
    }

    double log_val;
    double log_prior_genotype;
    double log_lik = DOUBLE_NEG_INF;
    double log_prior_norm = log(1 - gsl_ran_binomial_pdf(0, model_params.get_var_cp_prob(), total_cn));
    // We marginalize over the number of variant copies
    for (size_t var_cn = 1; var_cn <= total_cn; var_cn++) {
        xi = (1 - phi) * seq_err;
        if (var_cn == total_cn) {
            xi += phi * (1 - seq_err);
        } else {
            xi += phi * var_cn / total_cn;
        }
        log_val = log_binomial_pdf(datum->GetVariantReadCount(), xi, datum->GetReadCount());
        log_prior_genotype = log_binomial_pdf(var_cn, model_params.get_var_cp_prob(), total_cn);
        log_prior_genotype -= log_prior_norm;
        log_lik = log_add(log_lik, log_val + log_prior_genotype);
    }
    return log_lik;
}

double BulkLogLikWithCopyNumberProfile(const Node<BulkDatum, CloneTreeNodeParam> *node,
                                       const BulkDatum *datum,
                                       const ModelParams &model_params)
{
    if (datum->GetReadCount() == 0) {
        if (datum->GetVariantReadCount() > 0) {
            cerr << "Error in the data. Num reads: " << datum->GetReadCount() << ". Num variants: " << datum->GetVariantReadCount() << endl;
            exit(-1);
        } else {
            return 0.0;
        }
    }

    double seq_err = model_params.get_seq_error();
    if (node->get_parent_node() == 0) {
        // this node is the root, represents the healthy population
        return log(gsl_ran_binomial_pdf(datum->GetVariantReadCount(), seq_err, datum->GetReadCount()));
    }
    
    const CloneTreeNodeParam &param = node->get_node_parameter();
    auto cn_profile = datum->GetCopyNumberProfile();
    double phi = param.get_cellular_prev();
    double xi = 0.0;
    double log_val;
    double log_prior_genotype;

    // Initialize log_lik for total_cn = 0.
    // When total_cn = 0, variant is only seen (1 - phi) of the cells as a result
    // of sequencing error.
    xi = (1 - phi) * seq_err;
    double log_lik = log_binomial_pdf(datum->GetVariantReadCount(),
                                      xi,
                                      datum->GetReadCount());

    for (size_t total_cn = 1; total_cn < cn_profile.size(); total_cn++) {
        double log_prior_norm = log(1 - gsl_ran_binomial_pdf(0, model_params.get_var_cp_prob(), total_cn));
        double log_lik_cn = DOUBLE_NEG_INF;
        for (size_t var_cn = 1; var_cn <= total_cn; var_cn++) {
            xi = (1 - phi) * seq_err;
            if (var_cn == total_cn) {
                xi += phi * (1 - seq_err);
            } else if (var_cn == 0) {
                xi += phi * seq_err;
            } else {
                xi += phi * var_cn / total_cn;
            }
            log_val = log_binomial_pdf(datum->GetVariantReadCount(),
                                       xi,
                                       datum->GetReadCount());
            log_prior_genotype = log(gsl_ran_binomial_pdf(var_cn,
                  model_params.get_var_cp_prob(), total_cn)) - log_prior_norm;
            log_lik_cn = log_add(log_lik_cn, log_val + log_prior_genotype);
        }
        log_lik = log_add(log_lik, log_lik_cn + log(cn_profile[total_cn]));
    }
    return log_lik;
}

double BulkLogLikWithGenotype(const Node<BulkDatum, CloneTreeNodeParam> *node,
                              const BulkDatum *datum,
                              const ModelParams &model_params)
{
    if (datum->GetMajorCN() < datum->GetMinorCN()) {
        cerr << "Error: major copy number < minor copy number.\n";
        exit(-1);
    }
    if (datum->GetReadCount() == 0) {
        if (datum->GetVariantReadCount() > 0) {
            cerr << "Error in the data.\n";
            cerr << "Num reads: " << datum->GetReadCount() << ".\n";
            cerr << "Num variants: " << datum->GetVariantReadCount() << ".\n";
            exit(-1);
        } else {
            return 0.0;
        }
    }

    double seq_err = model_params.get_seq_error();
    if (node->get_parent_node() == 0) {
        // this node is the root, represents the healthy population
        return log(gsl_ran_binomial_pdf(datum->GetVariantReadCount(),
                                        seq_err,
                                        datum->GetReadCount()));
    }

    const CloneTreeNodeParam &param = node->get_node_parameter();
    double phi = param.get_cellular_prev();
    size_t major_cn = datum->GetMajorCN();
    size_t minor_cn = datum->GetMinorCN();
    size_t total_cn = major_cn + minor_cn;
    double xi = 0.0;
    
    if (total_cn == 0) {
        xi = seq_err;
        return log(gsl_ran_binomial_pdf(datum->GetVariantReadCount(),
                                        xi,
                                        datum->GetReadCount()));
    }
    
    double log_val;
    double log_lik = DOUBLE_NEG_INF;
    
    // Marginalize over the var_cn over minor_cn.
    for (size_t var_cn = 1; var_cn <= minor_cn; var_cn++) {
        xi = (1 - phi) * seq_err;
        if (var_cn == total_cn) {
            xi += phi * (1 - seq_err);
        } else {
            xi += phi * var_cn / total_cn;
        }
        log_val = log_binomial_pdf(datum->GetVariantReadCount(), xi, datum->GetReadCount());
        log_lik = log_add(log_lik, log_val);
    }
    // Marginalize over the var_cn over major_cn.
    for (size_t var_cn = 1; var_cn <= major_cn; var_cn++) {
        xi = (1 - phi) * seq_err;
        if (var_cn == total_cn) {
            xi += phi * (1 - seq_err);
        } else {
            xi += phi * var_cn / total_cn;
        }
        log_val = log_binomial_pdf(datum->GetVariantReadCount(), xi, datum->GetReadCount());
        log_lik = log_add(log_lik, log_val);
    }
    log_lik -= log(total_cn);
    return log_lik;
}

double update_cellular_prev_recursive(Node<BulkDatum,CloneTreeNodeParam> * node)
{
    unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &children = node->get_idx2child();
    double children_prev = 0.0;
    for (auto it = children.begin(); it != children.end(); ++it)
    {
        Node<BulkDatum,CloneTreeNodeParam> *child = it->second.second;
        children_prev += update_cellular_prev_recursive(child);
    }
    const CloneTreeNodeParam &param = node->get_node_parameter();
    double prevalence = param.get_clone_freq() + children_prev;
    ((CloneTreeNode *)node)->set_cellular_prev(prevalence);
    return prevalence;
}

double get_children_cellular_prev(Node<BulkDatum,CloneTreeNodeParam> *node)
{
    double children_prev = 0.0;
    unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &children = node->get_idx2child();
    for (auto it = children.begin(); it != children.end(); ++it)
    {
        Node<BulkDatum,CloneTreeNodeParam> *child = it->second.second;
        children_prev += child->get_node_parameter().get_cellular_prev();
    }
    return children_prev;
}

void update_params(Node<BulkDatum,CloneTreeNodeParam> *node, double curr_cellular_prev, double new_cellular_prev)
{
    auto clone_node = (CloneTreeNode *)node;
    double children_prev = get_children_cellular_prev(node);
    clone_node->set_cellular_prev(new_cellular_prev);
    clone_node->set_clone_freq(new_cellular_prev - children_prev);
    auto parent_node = (CloneTreeNode *)(node->get_parent_node());
    double parent_clone_freq = parent_node->get_node_parameter().get_clone_freq();
    double diff = new_cellular_prev - curr_cellular_prev;
    parent_node->set_clone_freq(parent_clone_freq - diff);
}

void update_params(vector<Node<BulkDatum,CloneTreeNodeParam> *> &nodes, double *new_clone_freq)
{
    // update clone frequencies
    for (size_t i = 0; i < nodes.size(); i++) {
        ((CloneTreeNode *)nodes[i])->set_clone_freq(new_clone_freq[i]);
    }
    // update cellular prevalences
    update_cellular_prev_recursive(nodes[0]);
}

void get_clone_freqs(TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &state,
                     double *clone_freqs)
{
    vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
    state.get_all_nodes(false, nodes);
    Node<BulkDatum,CloneTreeNodeParam> *node = 0;
    for (size_t i = 0; i < nodes.size(); i++)
    {
        node = nodes[i];
        clone_freqs[i] = node->get_node_parameter().get_cellular_prev();
        unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &children = node->get_idx2child();
        for (auto it = children.begin(); it != children.end(); ++it)
        {
            Node<BulkDatum,CloneTreeNodeParam> *child = it->second.second;
            clone_freqs[i] -= child->get_node_parameter().get_cellular_prev();
        }
    }
}

void sample_params_bottom_up(const gsl_rng *random,
                             size_t n_mh_iter,
                             TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &tree,
                             const ModelParams &params)
{
    // start from the leaf on up
    // update the clone freq using independent MH
    // note that likelihood over all mutations must be re-computed due to DP
    //double curr_log_lik = tree.compute_log_likelihood_bulk(params);
    double new_log_lik;
    vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
    tree.get_all_nodes(nodes);
    // Note: nodes[0] is the root node, we do not update its cellulare prevalence, as it is always set to 1
    for (int i = nodes.size() - 1; i > 0; i--)
    {
        //cout << nodes[i]->get_name() << endl;
        double curr_phi = nodes[i]->get_node_parameter().get_cellular_prev();
        double new_phi;
        double children_cellular_prev = get_children_cellular_prev(nodes[i]);
        double parent_clone_req = nodes[i]->get_parent_node()->get_node_parameter().get_clone_freq();
        double minimum = children_cellular_prev;
        double maximum = parent_clone_req + curr_phi;
        double range = maximum - minimum;
        double interval = range / n_mh_iter;
        vector<double> log_liks;
        vector<double> new_phis;
        double log_norm = DOUBLE_NEG_INF;
        for (size_t j = 0; j < n_mh_iter; j++) {
            new_phi = gsl_ran_flat(random, minimum + j*interval, minimum + (j+1)*interval);
            update_params(nodes[i], curr_phi, new_phi);
            new_log_lik = tree.compute_log_likelihood_bulk(params);
            log_liks.push_back(new_log_lik);
            new_phis.push_back(new_phi);
            log_norm = log_add(log_norm, new_log_lik);
            update_params(nodes[i], new_phi, curr_phi);
            //cout << new_phi << ", " << new_log_lik << endl;
        }
        if (log_norm != DOUBLE_NEG_INF) {
            normalize_destructively(log_liks);
            size_t idx = multinomial(random, log_liks);
            update_params(nodes[i], curr_phi, new_phis[idx]);
        }
    }
    
}

double sample_params_dirichlet(const gsl_rng *random,
                             size_t n_mh_iter,
                             TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &tree,
                             const ModelParams &params)
{
    double curr_log_lik = tree.compute_log_likelihood_bulk(params);
    double new_log_lik;
    
    // Perform MH to update the parameters
    // collect current parameters via tree traversal
    // sample from DP centered at current parameters multipled by constant factor
    double K = params.get_dir_conc_mult_factor();
    vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
    tree.get_all_nodes(nodes);
    size_t dim = nodes.size();
    double *curr_clone_freq = new double[dim];
    double *new_clone_freq = new double[dim];
    double *dir_concentration_params_curr = new double[dim];
    double *dir_concentration_params_new = new double[dim];
    
    get_clone_freqs(tree, curr_clone_freq);
    
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
        update_params(nodes, new_clone_freq);
        new_log_lik = tree.compute_log_likelihood_bulk(params);
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
            update_params(nodes, curr_clone_freq);
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

void cull(Node<BulkDatum,CloneTreeNodeParam> *root)
{
    descend_and_cull(root);
    descend_and_update_names(root);
    
    vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
    Node<BulkDatum,CloneTreeNodeParam>::breadth_first_traversal(root, nodes);
    
    // update the clone frequencies
    CloneTreeNode *curr_node;
    for (size_t i = 0; i < nodes.size(); i++) {
        curr_node = (CloneTreeNode *)nodes[i];
        unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &idx2child = curr_node->get_idx2child();
        double children_cellular_prev = 0.0;
        for (size_t idx = 0; idx < idx2child.size(); idx++) {
            CloneTreeNode *child_node = (CloneTreeNode *)idx2child[idx].second;
            children_cellular_prev += child_node->get_node_parameter().get_cellular_prev();
        }
        // cull does not alter cellular prevalence but it does change the clone frequncy because some children gets removed
        curr_node->set_clone_freq(curr_node->get_node_parameter().get_cellular_prev() - children_cellular_prev);
    }
}

//void cull(TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> *state)
//{
//    //state->cull();
//    descend_and_cull(state->get_root());
//    descend_and_update_names(state->get_root());
//    
//    vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
//    state->get_all_nodes(nodes);
//    unordered_set<Node<BulkDatum,CloneTreeNodeParam> *> node_set;
//    node_set.insert(nodes.begin(), nodes.end());
//    state->clear_cache(node_set);
//    
//    CloneTreeNode *curr_node;
//    for (size_t i = 0; i < nodes.size(); i++) {
//        curr_node = (CloneTreeNode *)nodes[i];
//        unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &idx2child = curr_node->get_idx2child();
//        double children_cellular_prev = 0.0;
//        for (size_t idx = 0; idx < idx2child.size(); idx++) {
//            CloneTreeNode *child_node = (CloneTreeNode *)idx2child[idx].second;
//            children_cellular_prev += child_node->get_node_parameter().get_cellular_prev();
//        }
//        // cull does not alter cellular prevalence but it does change the clone frequncy because some children gets removed
//        curr_node->set_clone_freq(curr_node->get_node_parameter().get_cellular_prev() - children_cellular_prev);
//    }
//}

