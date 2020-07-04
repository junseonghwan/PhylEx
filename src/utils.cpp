//
//  utils.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-07-30.
//

#include <gsl/gsl_sf.h>

#include "utils.hpp"

//void check_parameter_feasibility_helper(Node<BulkDatum,CloneTreeNodeParam> *node)
//{
//    double cellular_prev = node->get_node_parameter().get_cellular_prev();
//    double clone_freq = node->get_node_parameter().get_clone_freq();
//    assert( clone_freq >= 0.0 );
//    unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &idx2child = node->get_idx2child();
//    double children_cellular_prev = 0.0;
//    for (size_t i = 0; i < idx2child.size(); i++) {
//        Node<BulkDatum,CloneTreeNodeParam> *child_node = idx2child[i].second;
//        children_cellular_prev += child_node->get_node_parameter().get_cellular_prev();
//    }
//    assert( cellular_prev >= children_cellular_prev );
//    double diff = cellular_prev - children_cellular_prev;
//    assert( abs(diff - clone_freq) < 1e-3 );
//}

//void check_parameter_feasibility(Node<BulkDatum,CloneTreeNodeParam> *root, bool print = false)
//{
//    queue<Node<BulkDatum,CloneTreeNodeParam> *> q;
//    q.push(root);
//    Node<BulkDatum,CloneTreeNodeParam> *node = 0;
//    while (!q.empty()) {
//        node = q.front();
//        q.pop();
//        if (print)
//            cout << node->get_name() << ": " << node->get_node_parameter().get_cellular_prev() << ", " << node->get_node_parameter().get_clone_freq() << endl;
//        check_parameter_feasibility_helper(node);
//        unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &idx2child = node->get_idx2child();
//        for (size_t i = 0; i < idx2child.size(); i++) {
//            q.push(idx2child[i].second);
//        }
//    }
//}

bool check_clone_freq(Node<BulkDatum,CloneTreeNodeParam> *root)
{
    double sum_of_freqs = 0.0;
    queue<Node<BulkDatum,CloneTreeNodeParam> *> q;
    q.push(root);
    Node<BulkDatum,CloneTreeNodeParam> *node = 0;
    while (!q.empty()) {
        node = q.front();
        q.pop();
        sum_of_freqs += node->get_node_parameter().get_clone_freq();
        unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &idx2child = node->get_idx2child();
        for (size_t i = 0; i < idx2child.size(); i++) {
            q.push(idx2child[i].second);
        }
    }
    return(abs(sum_of_freqs - 1.0) < 1e-3);
}

double v_measure(double beta,
                 vector<size_t> true_labels,
                 vector<size_t> predicted_labels)
{
    if (true_labels.size() != predicted_labels.size()) {
        cerr << "True labels and predicted labels do not have the same length!" << endl;
        exit(-1);
    }
    
    size_t C0 = *min_element(true_labels.begin(), true_labels.end());
    size_t K0 = *min_element(predicted_labels.begin(), predicted_labels.end());
    
    size_t C = *max_element(true_labels.begin(), true_labels.end());
    size_t K = *max_element(predicted_labels.begin(), predicted_labels.end());
    
    size_t N = true_labels.size();
    
    // re-label from 0, ..., C-C0 (resp. K-K0)
    for (size_t n = 0; n < N; n++) {
        true_labels[n] -= C0;
        predicted_labels[n] -= K0;
    }
    C = C - C0 + 1;
    K = K - K0 + 1;
    
    // allocate and initialize the contigency table
    gsl_matrix *A = gsl_matrix_alloc(C, K);
    for (size_t c = 0; c < C; c++) {
        for (size_t k = 0; k < K; k++) {
            gsl_matrix_set(A, c, k, 0.0);
        }
    }
    
    // formulate the contingency table
    for (size_t n = 0; n < N; n++) {
        size_t c = true_labels[n];
        size_t k = predicted_labels[n];
        double val = gsl_matrix_get(A, c, k);
        gsl_matrix_set(A, c, k, val + 1);
    }
    
    vector<double> sum_over_c(K); // sum_over_c[k] = sum_{c} a_{ck}
    for (size_t k = 0; k < K; k++) {
        sum_over_c[k] = 0;
        for (size_t c = 0; c < C; c++) {
            sum_over_c[k] += gsl_matrix_get(A, c, k);
        }
    }
    vector<double> sum_over_k(C); // sum_over_k[c] = sum_{k} a_{ck}
    for (size_t c = 0; c < C; c++) {
        sum_over_k[c] = 0;
        for (size_t k = 0; k < K; k++) {
            sum_over_k[c] += gsl_matrix_get(A, c, k);
        }
    }
    
    double h_ck = 0.0;
    double a_ck = 0.0;
    for (size_t k = 0; k < K; k++) {
        for (size_t c = 0; c < C; c++) {
            a_ck = gsl_matrix_get(A, c, k);
            if (a_ck == 0)
                continue;
            h_ck += (a_ck / N) * log(a_ck / sum_over_c[k]);
        }
    }
    double h_c = 0.0;
    for (size_t c = 0; c < C; c++) {
        h_c += (sum_over_k[c]/N) * log(sum_over_k[c]/N);
    }
    
    double h = 1 - (-h_ck/-h_c);
    
    double h_kc = 0.0;
    for (size_t c = 0; c < C; c++) {
        for (size_t k = 0; k < K; k++) {
            a_ck = gsl_matrix_get(A, c, k);
            if (a_ck == 0)
                continue;
            h_kc += (a_ck / N) * log(a_ck / sum_over_k[c]);
        }
    }
    double h_k = 0.0;
    for (size_t k = 0; k < K; k++) {
        h_k += (sum_over_c[k]/N) * log(sum_over_c[k]/N);
    }
    
    double c = 1 - (-h_kc/-h_k);
    double v = (1 + beta) * h * c / (beta * h + c);
    cout << "h: " << h << endl;
    cout << "c: " << c << endl;
    cout << "V-measure: " << v << endl;
    
    gsl_matrix_free(A);
    
    return v;
}

double v_measure(double beta,
                 const vector<BulkDatum *> &data,
                 TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &true_state,
                 TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &sampled_state)
{
    // determine the classes and clusters
    vector<Node<BulkDatum,CloneTreeNodeParam> *> classes;
    true_state.get_all_nodes(true, classes);
    unordered_map<Node<BulkDatum,CloneTreeNodeParam> *, size_t> node2class;
    for (size_t c = 0; c < classes.size(); c++) {
        node2class[classes[c]] = c;
    }
    
    vector<Node<BulkDatum,CloneTreeNodeParam> *> clusters;
    sampled_state.get_all_nodes(true, clusters);
    unordered_map<Node<BulkDatum,CloneTreeNodeParam> *, size_t> node2cluster;
    for (size_t k = 0; k < clusters.size(); k++) {
        node2cluster[clusters[k]] = k;
    }
    
    size_t N = data.size();
    
    // formulate true and predicted labels
    vector<size_t> true_labels;
    vector<size_t> predicted_labels;
    Node<BulkDatum,CloneTreeNodeParam> *node = 0;
    for (size_t i = 0; i < N; i++) {
        node = true_state.get_node(data[i]);
        size_t c = node2class[node];
        
        node = sampled_state.get_node(data[i]);
        size_t k = node2cluster[node];
        
        true_labels.push_back(c);
        predicted_labels.push_back(k);
    }
    
    return v_measure(beta, true_labels, predicted_labels);
}

void compute_metrics(const vector<BulkDatum *> &data,
                     TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &true_state,
                     TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &sampled_state)
{
    // get ancestral matrices
    gsl_matrix *true_A = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::get_ancestral_matrix(true_state);
    gsl_matrix *sampled_A = TSSBState<BulkDatum, SingleCellData,CloneTreeNodeParam>::get_ancestral_matrix(sampled_state);
    // compute the differences between the ancestral matrices
    double err = 0.0;
    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data.size(); j++) {
            if (gsl_matrix_get(true_A, i, j) != gsl_matrix_get(sampled_A, i, j)) {
                err += 1;
            }
        }
    }
    cout << "Err: " << err << "/" << (data.size() * (data.size() - 1) / 2) << endl;
    
    // compute L1-loss on cellular prevalences for each of SNA
    Node<BulkDatum,CloneTreeNodeParam> *node = 0;
    double l1_loss = 0.0;
    for (size_t i = 0; i < data.size(); i++) {
        node = true_state.get_node(data[i]);
        double true_cellular_prev = node->get_node_parameter().get_cellular_prev();
        
        node = sampled_state.get_node(data[i]);
        double sampled_cellular_prev = node->get_node_parameter().get_cellular_prev();
        
        l1_loss += abs(sampled_cellular_prev - true_cellular_prev);
    }
    cout << "L1 loss: " << l1_loss/data.size() << endl;
    
    gsl_matrix_free(true_A);
    gsl_matrix_free(sampled_A);
    
    v_measure(1, data, true_state, sampled_state);
}

unordered_map<Locus *, size_t> *get_mutations()
{
    unordered_map<Locus *, size_t> *mut_map = new unordered_map<Locus *, size_t>();

    return mut_map;
}

void call_mutation(const unordered_map<Locus, LocusDatum*> &loci_data, const ModelParams &model_params, unordered_map<Locus, size_t> &mut_map, size_t max_cn)
{
    // log_priors stores prior belief on the number of variant copies
    // before assigning any mutation to a node, we assume that a node is likely to have
    // 0 variant copies (does not harbour the mutation), which is reflected by
    // assignment of largest possible prior prob to value 0: log_priors[0] = log(max_cn)
    // the rest are assigned a value in decreasing order
    vector<double> log_priors(max_cn);
    log_priors[0] = log(max_cn);
    for (size_t i = 1; i < max_cn; i++) {
        log_priors[i] = log(max_cn - i);
    }
    // normalize
    double log_norm = log_add(log_priors);
    for (size_t i = 0; i < max_cn; i++) {
        log_priors[i] -= log_norm;
    }

    for (auto it = loci_data.begin(); it != loci_data.end(); ++it) {
        const Locus &locus = it->first;
        const LocusDatum *loci_datum = it->second;
        size_t n_reads = loci_datum->get_n_total_reads();
        size_t n_var_reads = loci_datum->get_n_var_reads();

        if (n_reads == 0) {
            mut_map[locus] = 3; // 3: N/A
            continue;
        }
        
        double alpha = 0.0, beta = 0.0, log_sum, log_prior;
        vector<double> log_w(max_cn);
        for (size_t a = 0; a < max_cn; a++) {
            alpha = (a == 0) ? model_params.get_amp_error() : a;
            log_sum = DOUBLE_NEG_INF;
            log_prior = -log(max_cn - a + 1);
            for (size_t b = 0; b < max_cn - a; b++) {
                beta = (b == 0) ? model_params.get_amp_error() : b;
                log_sum = log_add(log_sum, log_beta_binomial_pdf(n_var_reads, n_reads, alpha, beta) + log_prior);
            }
            log_w[a] = log_sum + log_priors[a];
        }

        // convert to probabilities
        log_norm = log_add(log_w);
        double prob_no_mut = exp(log_w[0] - log_norm);
        if (prob_no_mut > 0.5) {
            mut_map[locus] = 0;
        } else {
            mut_map[locus] = 1;
        }
    }
}
