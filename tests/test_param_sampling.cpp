//
//  test_param_sampling.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-07-30.
//

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <unordered_map>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

#include "clone_node.hpp"
#include "model_params.hpp"
#include "node.hpp"
#include "numerical_utils.hpp"
#include "sampling_utils.hpp"
#include "single_cell.hpp"
#include "utils.hpp"

#include "simul_data.hpp"

void check_parameter_feasibility_helper(Node<BulkDatum,CloneTreeNodeParam> *node)
{
    double cellular_prev = node->get_node_parameter().get_cellular_prev();
    double clone_freq = node->get_node_parameter().get_clone_freq();
    BOOST_TEST( clone_freq >= 0.0 );
    unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &idx2child = node->get_idx2child();
    double children_cellular_prev = 0.0;
    for (size_t i = 0; i < idx2child.size(); i++) {
        Node<BulkDatum,CloneTreeNodeParam> *child_node = idx2child[i].second;
        children_cellular_prev += child_node->get_node_parameter().get_cellular_prev();
    }
    BOOST_TEST( cellular_prev >= children_cellular_prev );
    double diff = cellular_prev - children_cellular_prev;
    BOOST_TEST( abs(diff - clone_freq) < 1e-3 );
}

void check_parameter_feasibility(Node<BulkDatum,CloneTreeNodeParam> *root, bool print = false)
{
    queue<Node<BulkDatum,CloneTreeNodeParam> *> q;
    q.push(root);
    Node<BulkDatum,CloneTreeNodeParam> *node = 0;
    while (!q.empty()) {
        node = q.front();
        q.pop();
        if (print)
            cout << node->get_name() << ": " << node->get_node_parameter().get_cellular_prev() << ", " << node->get_node_parameter().get_clone_freq() << endl;
        check_parameter_feasibility_helper(node);
        unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &idx2child = node->get_idx2child();
        for (size_t i = 0; i < idx2child.size(); i++) {
            q.push(idx2child[i].second);
        }
    }
}

BOOST_AUTO_TEST_CASE( test_sampling_alpha0 )
{
    gsl_rng *random = generate_random_object(3);

    double ll = 0;
    double uu = 50;

    size_t n_reps = 10;
    size_t n_mcmc_iter = 5000;
    size_t n_mh_iter = 10;

    cout << "Testing sampling for alpha0." << endl;
    double chi2_alpha0 = 0.0;
    for (size_t n = 0; n < n_reps; n++)
    {
        cout << "Rep: " << (n + 1) << endl;
        double alpha0 = gsl_ran_flat(random, ll, uu);

        ModelParams model_params(alpha0, 2, 0.8, 0.01, 1);
        model_params.set_alpha0_bound(false, ll);
        model_params.set_alpha0_bound(true, uu);

        // sample tree
        CloneTreeNode *root = CloneTreeNode::create_root_node();
        root->set_nu_stick(bounded_beta(random, 1, model_params.alpha(root->get_name_vec())));
        root->sample_node_parameters(random, model_params, 0);
        auto state = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::construct_trivial_state(root, bulk_likelihood, sc_likelihood);

        // insert data, will generate new nodes
        size_t n_data = 10;
        for (size_t i = 0; i < n_data; i++) {
            Locus *loci = new Locus("1", 1);
            BulkDatum *datum = new BulkDatum("s" + to_string(i+1), *loci);
            state->insert_datum(random, datum, model_params);
        }

        // get all nodes
        vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
        state->get_all_nodes(nodes);

        // compute the likelihood at truth
        double log_lik_truth = log_prod_beta(nodes, model_params);

        // perform posterior sampling -- compute quantile for the true parameter
        double q_alpha0 = 0;

        // randomly initialize alpha0, lambda
        model_params.set_alpha0(gsl_ran_flat(random, ll, uu));
        size_t count = 0;
        for (size_t n = 0; n < n_mcmc_iter; n++) {
//            if (n % 100 == 0) {
//                cout << "Iter: " << n << endl;
//            }
            TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::sample_alpha0(random, n_mh_iter, model_params, nodes);
            if (n % 20 == 0) {
                if (model_params.get_alpha0() < alpha0) {
                    q_alpha0++;
                }
                count++;
            }
        }

        q_alpha0 /= count;

        double x_alpha0 = gsl_cdf_gaussian_Pinv(q_alpha0, 1);
        chi2_alpha0 += pow(x_alpha0, 2);

        cout << "Truth: " << alpha0 << ", " << log_lik_truth << endl;
        cout << "Quantile: " << q_alpha0 << endl;
        cout << "Std Normal: " << x_alpha0 << endl;
    }

    double pval = gsl_cdf_chisq_P(chi2_alpha0, n_reps);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );
}

BOOST_AUTO_TEST_CASE( test_sampling_lambda )
{
    cout << "Testing sampling for lambda." << endl;
    gsl_rng *random = generate_random_object(32);

    double ll = 0.2;
    double uu = 1;

    size_t n_reps = 10;
    size_t n_mcmc_iter = 5000;
    size_t n_mh_iter = 10;

    double chi2 = 0.0;
    for (size_t n = 0; n < n_reps; n++)
    {
        cout << "Rep: " << (n + 1) << endl;
        double lambda = gsl_ran_flat(random, ll, uu);

        ModelParams model_params(5, 2, lambda, 0.01, 1);
        CloneTreeNode *root = CloneTreeNode::create_root_node();
        root->set_nu_stick(0);
        root->sample_node_parameters(random, model_params, 0);
        auto state = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::construct_trivial_state(root, bulk_likelihood, sc_likelihood);

        // sample a tree
        // insert some data
        size_t n_data = 10;
        for (size_t i = 0; i < n_data; i++) {
            Locus *loci = new Locus("1", 1);
            BulkDatum *datum = new BulkDatum("s" + to_string(i+1), *loci);
            state->insert_datum(random, datum, model_params);
        }

        vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
        state->get_all_nodes(nodes);

        // compute the likelihood at truth
        double log_lik_truth = log_prod_beta(nodes, model_params);

        // perform posterior sampling -- compute quantile for the true parameter
        double q = 0;

        // randomly initialize alpha0, lambda
        model_params.set_lambda(gsl_ran_flat(random, ll, uu));
        size_t count = 0;
        for (size_t n = 0; n < n_mcmc_iter; n++) {
//            if (n % 100 == 0) {
//                cout << "Iter: " << n << endl;
//            }

            TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::sample_lambda(random, n_mh_iter, model_params, nodes);
            if (n % 20 == 0) {
                if (model_params.get_lambda() < lambda) {
                    q++;
                }
                count++;
            }
        }

        q /= count;

        double x = gsl_cdf_gaussian_Pinv(q, 1);
        chi2 += pow(x, 2);

        cout << "Truth: " << lambda << ", " << log_lik_truth << endl;
        cout << "Quantile: " << q << endl;
        cout << "Std Normal: " << x << endl;
    }

    double pval = gsl_cdf_chisq_P(chi2, n_reps);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );
}

BOOST_AUTO_TEST_CASE( test_sampling_gamma )
{
    cout << "Testing sampling for gamma." << endl;
    gsl_rng *random = generate_random_object(132);

    double ll = 0.5;
    double uu = 10;

    size_t n_reps = 10;
    size_t n_mcmc_iter = 5000;
    size_t n_mh_iter = 10;

    double chi2 = 0.0;
    for (size_t n = 0; n < n_reps; n++)
    {
        cout << "Rep: " << (n + 1) << endl;
        double gamma = gsl_ran_flat(random, ll, uu);

        ModelParams model_params(5, gamma, 0.8, 0.01, 1);
        CloneTreeNode *root = CloneTreeNode::create_root_node();
        root->set_nu_stick(0);
        root->sample_node_parameters(random, model_params, 0);
        auto state = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::construct_trivial_state(root, bulk_likelihood, sc_likelihood);

        // sample a tree
        // insert some data
        size_t n_data = 10;
        for (size_t i = 0; i < n_data; i++) {
            Locus *loci = new Locus("1", 1);
            BulkDatum *datum = new BulkDatum("s" + to_string(i+1), *loci);
            state->insert_datum(random, datum, model_params);
        }

        vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
        state->get_all_nodes(nodes);
        vector<double> psi_sticks;
        for (Node<BulkDatum,CloneTreeNodeParam> *node : nodes) {
            for (size_t i = 0; i < node->get_idx2child().size(); i++) {
                psi_sticks.push_back(node->get_idx2child()[i].first);
            }
        }


        // compute the likelihood at truth
        double log_lik_truth = log_prod_beta(nodes, model_params);

        // perform posterior sampling -- compute quantile for the true parameter
        double q = 0;

        // randomly initialize alpha0, lambda
        model_params.set_gamma(gsl_ran_flat(random, ll, uu));
        size_t count = 0;
        for (size_t n = 0; n < n_mcmc_iter; n++) {
//            if (n % 100 == 0) {
//                cout << "Iter: " << n << endl;
//            }

            TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::sample_gamma(random, n_mh_iter, model_params, psi_sticks);
            if (n % 20 == 0) {
                if (model_params.get_gamma() < gamma) {
                    q++;
                }
                count++;
            }
        }

        q /= count;

        double x = gsl_cdf_gaussian_Pinv(q, 1);
        chi2 += pow(x, 2);

        cout << "Truth: " << gamma << ", " << log_lik_truth << endl;
        cout << "Quantile: " << q << endl;
        cout << "Std Normal: " << x << endl;
    }

    double pval = gsl_cdf_chisq_P(chi2, n_reps);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );
}

BOOST_AUTO_TEST_CASE( test_sampling_parameters )
{
    cout << "Testing parameter sampling." << endl;
    gsl_rng *random = generate_random_object(9132);
    ModelParams model_params(1, 1, 0.3, 0.01, 100);

    SimulationConfig simul_config;
    size_t n_reps = 10;

    size_t n_mcmc_iter = 100000;

    simul_config.bulk_mean_depth = 1000;
    simul_config.max_cn = 6;
    simul_config.birth_rate = 0.1;
    simul_config.death_rate = 0.1;

    double chi2 = 0.0;
    
    CloneTreeNode *root = CloneTreeNode::create_root_node();
    //root->set_nu_stick(bounded_beta(random, 1, model_params.alpha("")));
    root->set_nu_stick(0.0);
    root->sample_node_parameters(random, model_params, 0);
    auto state = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::construct_trivial_state(root, bulk_likelihood, sc_likelihood);

    // sample a tree
    // insert some data
    size_t n_data = 10;
    vector<BulkDatum *> data;
    for (size_t i = 0; i < n_data; i++) {
        Locus *loci = new Locus("1", 1);
        BulkDatum *datum = new BulkDatum("s" + to_string(i+1), *loci);
        state->insert_datum(random, datum, model_params);
        data.push_back(datum);
    }
    cull(state->get_root());

    Node<BulkDatum,CloneTreeNodeParam> *node = state->get_node(data[0]);
    for (size_t n = 0; n < n_reps; n++)
    {
        cout << "Rep: " << (n + 1) << endl;

        // sample parameter for node from prior
        //node->sample_node_parameters(random, model_params, node->get_parent_node());
        check_parameter_feasibility(root);
        node->sample_node_parameters(random, model_params, node->get_parent_node());
        check_parameter_feasibility(root);
        // generate read counts using the truth parameters
        unordered_map<Node<BulkDatum,CloneTreeNodeParam>*, vector<pair<size_t, size_t> > > cn_profile;
        generate_bulk_data(random, simul_config, data, model_params, root, cn_profile);
        
        // get true clone frequency for the node
        double true_eta = node->get_node_parameter().get_clone_freq();
        double true_log_lik = state->compute_log_likelihood_bulk(model_params);

        // run MCMC to sample parameters
        // first, randomly initialize the parameter for the node
        node->sample_node_parameters(random, model_params, node->get_parent_node());
        check_parameter_feasibility(root);
        
        double curr_eta = node->get_node_parameter().get_clone_freq();
        double curr_log_lik = state->compute_log_likelihood_bulk(model_params);

        double q = 0.0;
        size_t count = 0;
        for (size_t i = 0; i < n_mcmc_iter; i++) {

            //sample_params(random, n_mh_iter, *state, model_params);
            node->sample_node_parameters(random, model_params, node->get_parent_node());
            curr_log_lik = state->compute_log_likelihood_bulk(model_params);

            // store parameter for one of the datum and compute quantile
            curr_eta = node->get_node_parameter().get_clone_freq();

            if ((i % 20) == 0) {
                // test feasibility of cellular prevalences and clone frequencies
                check_parameter_feasibility(root);

                if (true_eta > curr_eta) {
                    q++;
                }
                count++;
            }
        }

        q /= count;
        double x = gsl_cdf_gaussian_Pinv(q, 1);
        chi2 += pow(x, 2);

        cout << "Truth: " << true_eta << ", " << true_log_lik << endl;
        cout << "Quantile: " << q << endl;
        cout << "Std Normal: " << x << endl;
    }

    double pval = gsl_cdf_chisq_P(chi2, n_reps);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );
}

