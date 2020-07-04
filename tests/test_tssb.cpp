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
#include "utils.hpp"

using namespace std;

bool check_node_name_consistency(Node<BulkDatum,CloneTreeNodeParam> *node)
{
    const vector<size_t> &name_vec = node->get_name_vec();
    unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &idx2child = node->get_idx2child();
    for (size_t i = 0; i < idx2child.size(); i++) {
        const vector<size_t> &child_name_vec = idx2child[i].second->get_name_vec();
        size_t n = child_name_vec.size();
        for (size_t j = 0; j < n - 1; j++) {
            if (name_vec[j] != child_name_vec[j]) {
                return false;
            }
        }
        if (child_name_vec[n-1] != i) {
            return false;
        }
    }
    return true;
}

vector<double> collect_psi_sticks(TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &state)
{
    // collect psi-sticks
    vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
    state.get_all_nodes(false, nodes);
    vector<double> psi_sticks;
    for (auto it = nodes.begin(); it != nodes.end(); ++it)
    {
        unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &idx2child = (*it)->get_idx2child();
        for (auto it2 = idx2child.begin(); it2 != idx2child.end(); it2++)
        {
            psi_sticks.push_back(it2->second.first);
        }
    }
    return psi_sticks;
}

// returns true if any subtree (including the leaf is empty -- no datum)
bool has_empty_subtree(Node<BulkDatum,CloneTreeNodeParam> *node)
{
    Node<BulkDatum,CloneTreeNodeParam> *child = 0;
    unordered_map<size_t, pair<double, Node<BulkDatum,CloneTreeNodeParam> *> > &children = node->get_idx2child();
    if (children.size() == 0) {
        if (node->get_num_data() == 0) {
            return true;
        } else {
            return false;
        }
    }

    for (int i = 0; i < children.size(); i++) {

        child = children[i].second;

        bool ret = has_empty_subtree(child);
        // if any of the child has an empy subtree, return true immediately
        if (ret) {
            return true;
        }
    }
    
    // return false if non of the children are empty
    return false;
}

// set psi stick for the first child of node to val
// the rest of the children share the breaks the stick randomly
// the last child gets psi stick of 1 so that the all of the stick is broken up
void impute_psi_sticks(const gsl_rng *random, Node<BulkDatum,CloneTreeNodeParam> *node, double val)
{
    size_t n_children = node->get_idx2child().size();
    // set psi sticks to 0.5 for the first child
    node->get_idx2child()[0].first = 0.5;
    
    // set other psi-sticks to be equal but add up to use the remaining sticks (sum to 1)
    double stick_remaining = val;
    double beta;
    for (size_t i = 1; i < n_children; i++) {
        if (i == (n_children - 1)) {
            beta = 1;
        } else {
            beta = bounded_beta(random, 1, 1);
        }
        node->get_idx2child()[i].first = beta;
        stick_remaining -= beta*stick_remaining;
    }
}

BOOST_AUTO_TEST_CASE( test_tssb_state )
{
    gsl_rng *random = generate_random_object(3);
    ModelParams model_params(5, 2, 0.8, 0.01, 1);

    CloneTreeNode *root = CloneTreeNode::create_root_node();
    root->set_nu_stick(0.2);
    root->sample_node_parameters(random, model_params, 0);
    auto state = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::construct_trivial_state(root, 0, 0);

    Locus *loci = new Locus("1",1);
    BulkDatum *datum = new BulkDatum("s0", *loci);
    state->insert_datum(random, datum, model_params);

    BOOST_TEST( state->compute_log_likelihood_bulk(model_params) == 0.0);

    datum->set_n_variant_reads(20);
    datum->set_n_reads(100);
    datum->set_var_cn(1);
    datum->set_total_cn(2);
    datum->set_var_cn_observed(true);
    double phi = state->get_node(datum)->get_node_parameter().get_cellular_prev();
    double xi = (1 - phi) * model_params.get_seq_error() + phi * datum->get_var_cn()/datum->get_total_cn();
    double expected = log(gsl_ran_binomial_pdf(datum->get_n_variant_reads(), xi, datum->get_n_reads()));

    BOOST_TEST( state->compute_log_likelihood_bulk(model_params) == expected );

    // insert large number of data
    size_t n_data = 5000;
    for (size_t i = 0; i < n_data; i++) {
        Locus *loci = new Locus("1", 1);
        datum = new BulkDatum("s" + to_string(i+1), *loci);
        state->insert_datum(random, datum, model_params);
    }

    // perform distribution test of psi-sticks against Beta(1, gamma)
    // use CLT on the sample mean of the psi sticks and compare the mean against mean of Beta(1, gamma)
    vector<double> psi_sticks = collect_psi_sticks(*state);
    double mean = gsl_stats_mean(psi_sticks.data(), 1, psi_sticks.size());
    // use sd under null: (1 + gamma)/( (1+gamma)^2 * (1 + gamma + 1) )
    double gamma = model_params.get_gamma();
    double sd_under_null = sqrt((1 + gamma) / ( pow(1 + gamma, 2) * (1 + gamma + 1) ));
    double expected_mean = 1/(1 + gamma);
    double x = sqrt(psi_sticks.size()) * (mean - expected_mean)/sd_under_null;
    double p = gsl_cdf_gaussian_P(x, 1.0);
    BOOST_TEST((p > 0.025 && p < 0.975));

    // re-assign data and perform cull many times -- ensure correctness
    for (size_t i = 0; i < 100; i++) {
        state->resample_data_assignment(random, model_params);
        cull(state->get_root());
        BOOST_TEST( has_empty_subtree(root) == false );

        // get all nodes with data
        vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes_pre_cull;
        state->get_all_nodes(true, nodes_pre_cull);

        // cull the tree
        cull(state->get_root());
        // check that there are no empty subtree
        BOOST_TEST( has_empty_subtree(root) == false );

        vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes_post_cull;
        state->get_all_nodes(true, nodes_post_cull);

//        cout << "pre-cull nodes with data: " << nodes_pre_cull.size() << endl;
//        cout << "post-cull nodes with data: " << nodes_post_cull.size() << endl;
        BOOST_TEST( nodes_pre_cull.size() == nodes_post_cull.size() );

        // order of the nodes should be the same for two trees and hence, we can compare Node pointers
        size_t total1 = 0, total2 = 0;
        for (size_t i = 0; i < nodes_pre_cull.size(); i++) {
            total1 += nodes_pre_cull[i]->get_num_data();
            total2 += nodes_pre_cull[i]->get_num_data();
            BOOST_TEST( nodes_pre_cull[i] == nodes_post_cull[i] );
        }
        BOOST_TEST( total1 == total2 );
    }
}

BOOST_AUTO_TEST_CASE( test_resample_stick_orders )
{
    gsl_rng *random = generate_random_object(12);
    ModelParams model_params(5, 2, 0.8, 0.01, 1);

    CloneTreeNode *root = CloneTreeNode::create_root_node();
    root->set_nu_stick(0);
    root->sample_node_parameters(random, model_params, 0);
    auto state = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::construct_trivial_state(root, bulk_likelihood, sc_likelihood);

    // insert some data
    size_t n_data = 10;
    for (size_t i = 0; i < n_data; i++) {
        Locus *loci = new Locus("1", 1);
        BulkDatum *datum = new BulkDatum("s" + to_string(i+1), *loci);
        state->insert_datum(random, datum, model_params);
    }

    BOOST_TEST( check_node_name_consistency(root) );
    cull(state->get_root());
    BOOST_TEST( check_node_name_consistency(root) ); // test that cull maintains the consistency of the node names
    vector<double> psi_sticks_pre = collect_psi_sticks(*state);
    sort(psi_sticks_pre.begin(), psi_sticks_pre.end());

    // reorder the sticks
    state->reorder_sticks(random);
    BOOST_TEST( check_node_name_consistency(root) ); // test that cull maintains the consistency of the node names

    vector<double> psi_sticks_post = collect_psi_sticks(*state);
    sort(psi_sticks_post.begin(), psi_sticks_post.end());

    // the number of psi-sticks should not have changed
    BOOST_TEST(psi_sticks_pre.size() == psi_sticks_post.size());
    // check that psi stick lengths did not change
    for (size_t i = 0; i < psi_sticks_pre.size(); i++) {
        BOOST_TEST( psi_sticks_pre[i] == psi_sticks_post[i] );
    }

    // set psi-sticks manually to test for any boundary conditions
    size_t n_children = root->get_idx2child().size();
    while (n_children <= 2) {
        // assign more data points
        Locus *loci = new Locus("1", 1);
        BulkDatum *datum = new BulkDatum("s" + to_string(n_data), *loci);
        state->insert_datum(random, datum, model_params);
        n_data++;
        n_children = root->get_idx2child().size();
    }

    // reorder the sticks many times, we expect child0 to appear first about half of the times
    // define a Bernoulli RV, that takes a value of 1 if child0 is the first child of the root and 0 otherwise
    // we can compute the Monte Carlo mean and sd and use CLT to compare against expected mean of 0.5
    size_t n_reps = 10000;
    double mean = 0.0;
    for (size_t i = 0; i < n_reps; i++) {

        impute_psi_sticks(random, root, 0.5);
        Node<BulkDatum,CloneTreeNodeParam> *child0 = root->get_idx2child()[0].second;
        state->reorder_sticks(random);
        BOOST_TEST( check_node_name_consistency(root) );

        Node<BulkDatum,CloneTreeNodeParam> *node = root->get_idx2child()[0].second;
        if (node == child0) {
            mean += 1;
        }
    }

    mean /= n_reps;
    //double sd = sqrt(mean * (1 - mean) / n_reps);
    double expected_mean = 0.5;
    double sd_under_null = sqrt(expected_mean * (1 - expected_mean)) / sqrt(n_reps);
    double x = (mean - expected_mean)/sd_under_null;
    double pval = gsl_cdf_gaussian_P(x, 1);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );
}

BOOST_AUTO_TEST_CASE( test_resample_sticks )
{
    gsl_rng *random = generate_random_object(4321);
    ModelParams model_params(5, 2, 0.8, 0.01, 1);

    CloneTreeNode *root = CloneTreeNode::create_root_node();
    root->set_nu_stick(0);
    auto state = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::construct_trivial_state(root, bulk_likelihood, sc_likelihood);

    // insert some data to random nodes in the tree
    vector<BulkDatum *> data;
    size_t n_data = 10;
    for (size_t i = 0; i < n_data; i++) {
        Locus *loci = new Locus("1", 1);
        BulkDatum *datum = new BulkDatum("s" + to_string(i+1), *loci);
        state->insert_datum(random, datum, model_params);
        data.push_back(datum);
    }

    // move all data to first child of the root
    // at least one child is expected to exist because the nu stick for the root is 0
    // therefore, each datum must be assigned to non-root node
    if (root->get_idx2child().size() == 0) {
        BOOST_ERROR("Number of children equal to 0.");
    }
    Node<BulkDatum,CloneTreeNodeParam> *child0 = root->get_idx2child()[0].second;
    for (size_t i = 0; i < n_data; i++) {
        state->move_datum(child0, data[i], model_params);
    }

    // repeatedly sample sticks,
    size_t n_reps = 10000;
    vector<double> nu_sticks;
    vector<double> psi_sticks;
    for (size_t i = 0; i < n_reps; i++) {
        state->descend_and_sample_sticks(random, root, model_params);
        nu_sticks.push_back(child0->get_nu_stick());
        psi_sticks.push_back(root->get_idx2child()[0].first);
    }

    // check that the distribution of the nu stick should be Beta(n_data + 1, alpha(child0->get_name()))
    double alpha = n_data + 1;
    double beta = model_params.alpha(child0->get_name_vec());
    double expected_mean = alpha / (alpha + beta);
    double mean = gsl_stats_mean(nu_sticks.data(), 1, n_reps);
    double sd_under_null = sqrt((alpha * beta) / ( pow(alpha + beta, 2) * (alpha + beta + 1) ));
    double x = sqrt(n_reps) * (mean - expected_mean) / sd_under_null;
    double pval = gsl_cdf_gaussian_P(x, 1);
//    printf("%f vs %f\n", expected_mean, mean);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );

    // check the distribution of the psi stick is Beta(n_data + 1, gamma)
    beta = model_params.get_gamma();
    expected_mean = alpha / (alpha + beta);
    mean = gsl_stats_mean(psi_sticks.data(), 1, n_reps);
    sd_under_null = sqrt((alpha * beta) / ( pow(alpha + beta, 2) * (alpha + beta + 1) ));
    x = sqrt(n_reps) * (mean - expected_mean) / sd_under_null;
    pval = gsl_cdf_gaussian_P(x, 1);
//    printf("%f vs %f\n", expected_mean, mean);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );

    // assign some data to child of child0, check distribution of nu-sticks
    size_t n_grand_children = child0->get_idx2child().size();

    // set psi-stick for child0 to be 0.999 and nu stick to be 0.01
    // this ensures that when we re-assign data, a child of child0 will be selected with high probability
    root->get_idx2child()[0].first = 0.999;
    child0->set_nu_stick(0.01);
    double u;
    while (n_grand_children <= 0) {
        // call re-assign data
        for (size_t i = 0; i < n_data; i++) {
            u = gsl_ran_flat(random, 0, 1);
            Node<BulkDatum,CloneTreeNodeParam> *node = Node<BulkDatum,CloneTreeNodeParam>::find_node(random, u, root, model_params);
            state->move_datum(node, data[i], model_params);
        }
        n_grand_children = child0->get_idx2child().size();
    }

    // keep some data at the child0, move rest to the first child of child0
    Node<BulkDatum,CloneTreeNodeParam> *child0_0 = child0->get_idx2child()[0].second;
    size_t n_data_to_child0 = 3;
    for (size_t i = 0; i < n_data; i++) {
        if (i < n_data_to_child0) {
            state->move_datum(child0, data[i], model_params);
        } else {
            state->move_datum(child0_0, data[i], model_params);
        }
    }

    //state->cull();
    //BOOST_TEST( state->get_num_nodes() == 3 );

    // resample sticks
    for (size_t i = 0; i < n_reps; i++)
    {
        state->descend_and_sample_sticks(random, root, model_params);
        nu_sticks[i] = child0->get_nu_stick();
    }
    alpha = n_data_to_child0 + 1;
    beta = (n_data - n_data_to_child0) + model_params.alpha(child0->get_name_vec());
    expected_mean = alpha / (alpha + beta);
    mean = gsl_stats_mean(nu_sticks.data(), 1, n_reps);
    sd_under_null = sqrt((alpha * beta) / ( pow(alpha + beta, 2) * (alpha + beta + 1) ));
    x = sqrt(n_reps) * (mean - expected_mean) / sd_under_null;
    pval = gsl_cdf_gaussian_P(x, 1);
//    printf("%f vs %f\n", expected_mean, mean);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );

    // assign some data to sibling of child0, check distribution of psi-sticks
    size_t n_children = root->get_idx2child().size();
    root->set_psi_stick(0, 0.0001); // set psi stick to a very small number so that probability of being assigned to child0 is small
    // get at least 3 children
    while (n_children <= 2) {
        for (size_t i = 0; i < n_data; i++) {
            u = gsl_ran_flat(random, 0, 1);
            Node<BulkDatum,CloneTreeNodeParam> *node = Node<BulkDatum,CloneTreeNodeParam>::find_node(random, u, root, model_params);
            state->move_datum(node, data[i], model_params);
        }
        n_children = root->get_idx2child().size();
    }

    // keep some data in child0, move the rest of the data to child1 and child2
    Node<BulkDatum,CloneTreeNodeParam> *child1 = root->get_idx2child()[1].second;
    Node<BulkDatum,CloneTreeNodeParam> *child2 = root->get_idx2child()[2].second;
    // partition the data (specify end points of the intervals)
    n_data_to_child0 = 3;
    size_t n_data_to_child0_0 = 5;
    size_t n_data_to_child1 = 8;
    for (size_t i = 0; i < n_data; i++) {
        if (i < n_data_to_child0) {
            state->move_datum(child0, data[i], model_params);
        } else if (i < n_data_to_child0_0) {
            state->move_datum(child0_0, data[i], model_params);
        } else if (i <n_data_to_child1) {
            state->move_datum(child1, data[i], model_params);
        } else {
            // remaining goes to child2
            state->move_datum(child2, data[i], model_params);
        }
    }

    cull(state->get_root());
    BOOST_TEST( state->get_num_nodes() == 5 );

    // resample sticks
    // store psi-sticks for child0 into psi_sticks, child1 into psi_sticks1, child2 into psi_sticks2
    vector<double> psi_sticks1, psi_sticks2;
    for (size_t i = 0; i < n_reps; i++)
    {
        state->descend_and_sample_sticks(random, root, model_params);
        psi_sticks[i] = root->get_idx2child()[0].first;
        psi_sticks1.push_back(root->get_idx2child()[1].first);
        psi_sticks2.push_back(root->get_idx2child()[2].first);
    }

    // the psi-stick for child0 should depend on the number of data points that went down that branch (n_data to 0 and 0_0)
    alpha = n_data_to_child0_0 + 1;
    beta = (n_data - n_data_to_child0_0) + model_params.get_gamma();
    expected_mean = alpha / (alpha + beta);
    mean = gsl_stats_mean(psi_sticks.data(), 1, n_reps);
    sd_under_null = sqrt((alpha * beta) / ( pow(alpha + beta, 2) * (alpha + beta + 1) ));
    x = sqrt(n_reps) * (mean - expected_mean) / sd_under_null;
    pval = gsl_cdf_gaussian_P(x, 1);
//    printf("%f vs %f\n", expected_mean, mean);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );

    alpha = (n_data_to_child1 - n_data_to_child0_0) + 1;
    beta = (n_data - n_data_to_child1) + model_params.get_gamma();
    expected_mean = alpha / (alpha + beta);
    mean = gsl_stats_mean(psi_sticks1.data(), 1, n_reps);
    sd_under_null = sqrt((alpha * beta) / ( pow(alpha + beta, 2) * (alpha + beta + 1) ));
    x = sqrt(n_reps) * (mean - expected_mean) / sd_under_null;
    pval = gsl_cdf_gaussian_P(x, 1);
//    printf("%f vs %f\n", expected_mean, mean);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );

    alpha = (n_data - n_data_to_child1) + 1;
    beta = model_params.get_gamma();
    expected_mean = alpha / (alpha + beta);
    mean = gsl_stats_mean(psi_sticks2.data(), 1, n_reps);
    sd_under_null = sqrt((alpha * beta) / ( pow(alpha + beta, 2) * (alpha + beta + 1) ));
    x = sqrt(n_reps) * (mean - expected_mean) / sd_under_null;
    pval = gsl_cdf_gaussian_P(x, 1);
//    printf("%f vs %f\n", expected_mean, mean);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );
}

// test find_node and locate_child
BOOST_AUTO_TEST_CASE( test_find_node )
{
    // sample a tree, set nu and psi sticks
    // check that the expected node is returned

    gsl_rng *random = generate_random_object(1);
    ModelParams model_params(1, 0.5, 0.8, 0.01, 1);

    CloneTreeNode *root = CloneTreeNode::create_root_node();
    root->set_nu_stick(0.2);
    auto state = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::construct_trivial_state(root, bulk_likelihood, sc_likelihood);
    double u = 0.1;
    BOOST_TEST( CloneTreeNode::find_node(random, u, root, model_params) == root );

    // insert some data to random node in the tree
    vector<BulkDatum *> data;
    size_t n_data = 10;
    for (size_t i = 0; i < n_data; i++) {
        Locus *loci = new Locus("1", 1);
        BulkDatum *datum = new BulkDatum("s" + to_string(i+1), *loci);
        state->insert_datum(random, datum, model_params);
        data.push_back(datum);
    }

    root->get_idx2child()[0].first = 0.0001;
    size_t n_children = root->get_idx2child().size();
    while (n_children <= 1) {
        for (size_t i = 0; i < n_data; i++) {
            u = gsl_ran_flat(random, 0, 1);
            Node<BulkDatum,CloneTreeNodeParam> *node = Node<BulkDatum,CloneTreeNodeParam>::find_node(random, u, root, model_params);
            state->move_datum(node, data[i], model_params);
        }
        n_children = root->get_idx2child().size();
    }

    root->get_idx2child()[0].first = 0.8;
    root->get_idx2child()[1].first = 0.2;

    Node<BulkDatum,CloneTreeNodeParam> *child0 = root->get_idx2child()[0].second;
    Node<BulkDatum,CloneTreeNodeParam> *child1 = root->get_idx2child()[1].second;

    root->set_nu_stick(0.2);
    child0->set_nu_stick(0.5);
    child1->set_nu_stick(0.6);

    // first stick length is 0.8 of (1-0.2) = 0.64 so [0.2, 0.84)
    // of that, the first child is responsible for interval [0.2, 0.2 + (1-0.2)*0.8*0.5) = [0.2, 0.52)
    u = 0.51999;
    BOOST_TEST( CloneTreeNode::find_node(random, u, root, model_params) == child0 );
    u = 0.2;
    BOOST_TEST( CloneTreeNode::find_node(random, u, root, model_params) == child0 );
    u = 0.521;
    BOOST_TEST( CloneTreeNode::find_node(random, u, root, model_params) != child0 );

    // second stick length is 0.2 of (1-0.2)*(1-0.8) = 0.032 so [0.84, 0.872)
    // second child is responsible for interval [0.84, 0.84 + (1-0.2)*(1-0.8)*0.2*0.6] = [0.84, 0.8592)
    u = 0.8591999;
    BOOST_TEST( CloneTreeNode::find_node(random, u, root, model_params) == child1 );
    u = 0.8401;
    BOOST_TEST( CloneTreeNode::find_node(random, u, root, model_params) == child1 );
    u = 0.85921;
    BOOST_TEST( CloneTreeNode::find_node(random, u, root, model_params) != child1 );
}

BOOST_AUTO_TEST_CASE( test_model_params )
{
    double alpha0 = 25;
    double gamma = 2;
    double lambda = 1;
    double seq_err = 0.01;
    ModelParams model_params(alpha0, gamma, lambda, seq_err, 1);
    BOOST_TEST( (model_params.get_alpha0() == alpha0) );
    BOOST_TEST( (model_params.get_lambda() == lambda) );
    BOOST_TEST( (model_params.get_gamma() == gamma) );
    BOOST_TEST( model_params.check_bounds() );

    model_params.set_alpha0_bound(true, 24.9);
    BOOST_TEST( (model_params.check_bounds() == false) );
    alpha0 = 20;
    model_params.set_alpha0(alpha0);
    BOOST_TEST( (model_params.check_bounds() == true) );

    model_params.set_lambda(1.2);
    BOOST_TEST( (model_params.check_bounds() == false) );
    model_params.set_lambda(lambda);
    BOOST_TEST( model_params.get_lambda() == 1.0 );

    vector<size_t> name_vec;
    name_vec.push_back(0);
    name_vec.push_back(1);
    name_vec.push_back(2);
    double val = model_params.alpha(name_vec);
    BOOST_TEST( (val == alpha0) );

    lambda = 0.8;
    model_params.set_lambda(lambda);
    val = model_params.alpha(name_vec);
    BOOST_TEST( (val == alpha0 * pow(lambda, 2)) );

}


