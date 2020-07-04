//
//  test_node.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-07-29.
//

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <unordered_map>
#include <unordered_set>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

#include "clone_node.hpp"
#include "loci.hpp"
#include "model_params.hpp"
#include "node.hpp"
#include "numerical_utils.hpp"
#include "sampling_utils.hpp"

using namespace std;

BOOST_AUTO_TEST_CASE( test_node )
{
    CloneTreeNode *node = CloneTreeNode::create_root_node();
    BOOST_TEST( node->get_name() == "0");
    
    size_t n_data = 5;
    BulkDatum *datum;
    Locus *loci;
    vector<BulkDatum *> data;
    for (size_t i = 0; i < n_data; i++) {
        loci = new Locus("1", 1);
        datum = new BulkDatum("s" + to_string(i), *loci);
        node->add_datum(datum);
        data.push_back(datum);
    }
    
    BOOST_TEST( node->get_num_data() == n_data);
    
    double nu = 0.5;
    node->set_nu_stick(nu);
    BOOST_TEST( node->get_nu_stick() == nu);
    
    BOOST_TEST( (node->get_parent_node() == 0) );

    BOOST_TEST( node->get_num_data() == node->get_data().size() );
    for (size_t i = 0; i < n_data; i++) {
        BOOST_TEST( node->get_data().count(data[i]) > 0 );
    }

    // remove datum
    node->remove_datum(data[0]);
    BOOST_TEST( (node->get_data().count(data[0]) == 0) );
    BOOST_TEST( (node->contains_datum(data[0]) == false) );

    BOOST_TEST( (node->get_num_children() == 0) );
    BOOST_TEST( node->is_root() );

//    BOOST_TEST( (node->get_parent_string(node->get_name()) == "0_2") );
//    BOOST_TEST( (node->form_node_string(node->get_name(), 1) == "0_2_3_1") );
}

BOOST_AUTO_TEST_CASE( test_clone_node )
{
    CloneTreeNode *root = CloneTreeNode::create_root_node();
    Node<BulkDatum,CloneTreeNodeParam> *child0 = root->spawn_child(0.5);
    Node<BulkDatum,CloneTreeNodeParam> *child1 = root->spawn_child(0.5);
    
    BOOST_TEST( (child0->get_parent_node() == root) );
    BOOST_TEST( (child0->get_name() == "0_0") );
    BOOST_TEST( (child1->get_name() == "0_1") );
    BOOST_TEST( (root->get_num_children() == 2) );
    BOOST_TEST( (root->get_child(0).second == child0) );
    BOOST_TEST( (root->get_child(1).second == child1) );
    BOOST_TEST ( (child0->is_root() == false) );

    gsl_rng *random = generate_random_object(1);
    ModelParams model_params(5, 2, 0.5, 0.01, 1);

    vector<double> params;
    size_t n_reps = 5000;
    root->sample_node_parameters(random, model_params, 0);
    root->set_cellular_prev(0.5);
    root->set_clone_freq(0.5);
    for (size_t i = 0; i < n_reps; i++) {
        child0->sample_node_parameters(random, model_params, root);
        params.push_back(child0->get_node_parameter().get_cellular_prev());
    }

    double mean = gsl_stats_mean(params.data(), 1, n_reps);
    double expected_mean = 0.25;
    double sd_under_null = sqrt((1./12) * pow(0.5, 2));
    double x = sqrt(n_reps) * (mean - expected_mean)/sd_under_null;
    double pval = gsl_cdf_gaussian_P(x, 1);
    BOOST_TEST( (pval > 0.025 && pval < 0.975) );
}

BOOST_AUTO_TEST_CASE( test_clone_node_param )
{
    CloneTreeNodeParam p;
    double eta = 0.2;
    double phi = 0.7;
    p.set_clone_freq(eta);
    p.set_cellular_prev(phi);
    BOOST_TEST( (p.get_clone_freq() == eta) );
    BOOST_TEST( (p.get_cellular_prev() == phi) );
    BOOST_TEST( p.is_consistent() );    
}

