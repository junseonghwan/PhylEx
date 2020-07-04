#include <boost/test/unit_test.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "clone_node.hpp"
#include "numerical_utils.hpp"
#include "model_params.hpp"
#include "tssb_state.hpp"

// 1. generate rate matrix, check that the components are correctly filled
// 2. compute transition matrix, ensure that the rows sum to 1
// 3. compute initial probabilities for copy numbers
// 4. compute expectation for various simple cases
// 4a. assignment to leaf node
// 4b. assignment to a node with one child
// 4c. assignment to a node with two children
// 4d. chain of length 3 (3 nodes)
// 4e. tree with height = 2, and there are 2 leaf nodes
// 4f. tree with height = 2, and there are 2 internal nodes, each of them has one child

//BOOST_AUTO_TEST_CASE( test_rate_matrix_construction )
//{
//    double b_rate = 0.1;
//    double d_rate = 0.2;
//    size_t N = 5;
//    gsl_matrix *Q = gsl_matrix_alloc(N, N);
//    construct_rate_matrix(b_rate, d_rate, Q);
//    BOOST_TEST( gsl_matrix_get(Q, 0, 0) == 0.0 );
//    BOOST_TEST( gsl_matrix_get(Q, 0, 1) == 0.0 );
//    
//    BOOST_TEST( gsl_matrix_get(Q, N-1, N-1) == -d_rate );
//    BOOST_TEST( gsl_matrix_get(Q, N-1, N-2) == d_rate );
//
//    BOOST_TEST( gsl_matrix_get(Q, 3, 2) == d_rate );
//    BOOST_TEST( gsl_matrix_get(Q, 3, 3) == -(b_rate + d_rate) );
//    BOOST_TEST( gsl_matrix_get(Q, 3, 4) == b_rate );
//
//}

//BOOST_AUTO_TEST_CASE( test_transition_matrix_construction )
//{
//    double b_rate = 0.1;
//    double d_rate = 0.2;
//    size_t N = 5;
//
//    gsl_matrix *Q = gsl_matrix_alloc(N, N);
//    gsl_matrix *Q_test = gsl_matrix_alloc(N, N);
//
//    gsl_matrix *U = gsl_matrix_alloc(N, N);
//    gsl_matrix *U_inv = gsl_matrix_alloc(N, N);
//
//    gsl_vector *d = gsl_vector_alloc(N);
//    gsl_matrix *D = gsl_matrix_alloc(N, N);
//
//    gsl_matrix *P = gsl_matrix_alloc(N, N);
//
//    construct_rate_matrix(b_rate, d_rate, Q);
//    eigen(Q, U, U_inv, d);
//
//    // check that U*U_inv = I
//    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, U, U_inv, 0.0, D);
//    for (size_t i = 0; i < N; i++) {
//        for (size_t j = 0; j < N; j++) {
//            if (i == j) {
//                BOOST_TEST( abs(gsl_matrix_get(D, i, j) - 1) < 1e-4 );
//            } else {
//                BOOST_TEST( gsl_matrix_get(D, i, j) < 1e-4 );
//            }
//        }
//    }
//
//    for (size_t i = 0; i < N; i++) {
//        for (size_t j = 0; j < N; j++) {
//            gsl_matrix_set(D, i, j, 0.0);
//        }
//        gsl_matrix_set(D, i, i, gsl_vector_get(d, i));
//    }
//
//    // check that UDU^{-1} = Q
//    gsl_matrix *temp = gsl_matrix_alloc(N, N);
//    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, D, 0, temp);
//    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp, U_inv, 0, Q_test);
//    for (size_t i = 0; i < N; i++) {
//        for (size_t j = 0; j < N; j++) {
//            double val1 = gsl_matrix_get(Q, i, j);
//            double val2 = gsl_matrix_get(Q_test, i, j);
//            //cout << val1 << ", " << val2 << endl;
//            BOOST_TEST( abs(val1 - val2) < 1e-4 );
//        }
//    }
//
//    construct_transition_matrix(1, U, U_inv, d, P);
//
//    // check that the rows sum to 1
//    for (size_t i = 0; i < N; i++) {
//        double sum = 0.0;
//        for (size_t j = 0; j < N; j++) {
//            sum += gsl_matrix_get(P, i, j);
//        }
//        BOOST_TEST( abs(sum - 1) < 1e-4 );
//    }
//
//    // test the elements of the transition matrix against hand computed values
//    ModelParams params(1, 1, 1, 0.01, 0.1);
//    params.set_birth_rate(0.1);
//    params.set_death_rate(0.2);
//    params.recompute_transition_matrix();
//    size_t parent_cn = 3;
//    size_t child_cn = 1;
//    double prob = params.compute_transition_prob(parent_cn, child_cn);
//    double exp_prob = 0.01491688;
//    BOOST_TEST( abs(prob - exp_prob) < 1e-4 );
//
//    cout << "Transition matrix test done!" << endl;
//}
//
//BOOST_AUTO_TEST_CASE( test_initial_transition_prob )
//{
//    ModelParams params(1, 1, 1, 0.01, 0.1);
//    params.set_birth_rate(0.1);
//    params.set_death_rate(0.1);
//    params.recompute_transition_matrix();
//
//    double b = 5;
//    size_t chi_r = 2, chi_v = 2, parent_cn = 2;
//    double prob = params.compute_initial_prob(b, chi_r, chi_v, parent_cn);
//    double exp_prob = 0.03034709;
//    BOOST_TEST( abs(prob - exp_prob) < 1e-4 );
//
//    // test that initial prob sum to 1
//    double sum = 0.0;
//    for (size_t i = 0; i < params.get_max_cn(); i++) {
//        for (size_t j = 0; j < params.get_max_cn(); j++) {
//            sum += params.compute_initial_prob(b, i, j, parent_cn);
//        }
//    }
//    BOOST_TEST( abs(1 - sum) < 1e-4 );
//}

//BOOST_AUTO_TEST_CASE( test_dynamic_program_a )
//{
//    // generate TSSB object, create a small tree, insert data
//    ModelParams model_params(1, 0.5, 0.8, 0.05, 1);
//    model_params.set_birth_rate(0.1);
//    model_params.set_death_rate(0.1);
//    model_params.recompute_transition_matrix();
//
//    CloneTreeNode *root = CloneTreeNode::create_root_node();
//    auto state = TSSBState<BulkDatum,SingleCellDNA,CloneTreeNodeParam>::construct_trivial_state(root);
//
//    CloneTreeNode *child0 = (CloneTreeNode *)root->spawn_child(0.5);
//
//    // assign datum to each of the 3 nodes
//    Locus *s1 = new Locus(1, 1);
//    BulkDatum d1("s1", s1, 10, 100, 2.0);
//
//    state->move_datum(child0, &d1, model_params);
//
//    double freq = 0.2;
//    child0->set_cellular_prev(freq);
//    child0->set_clone_freq(freq);
//
//    root->set_cellular_prev(1);
//    root->set_clone_freq(1 - freq);
//
//    double xi = child0->compute_expectation_of_xi(model_params);
//    double expected_val = 0.1731112;
//    BOOST_TEST( abs(xi - expected_val) < 1e-3 );
//}
//
//BOOST_AUTO_TEST_CASE( test_dynamic_program_b )
//{
//    // generate TSSB object, create a small tree, insert data
//    double seq_err = 0.05;
//    double b_rate = 0.1, d_rate = 0.1;
//    ModelParams model_params(1, 0.5, 0.8, 0.01, seq_err);
//    model_params.set_birth_rate(b_rate);
//    model_params.set_death_rate(d_rate);
//    model_params.recompute_transition_matrix();
//
//    CloneTreeNode *root = CloneTreeNode::create_root_node();
//    auto state = TSSBState<BulkDatum,SingleCellDNA,CloneTreeNodeParam>::construct_trivial_state(root);
//
//    CloneTreeNode *child0 = (CloneTreeNode *)root->spawn_child(0.5);
//    CloneTreeNode *child1 = (CloneTreeNode *)child0->spawn_child(0.5);
//
//    // assign datum to each of the 3 nodes
//    Locus *s1 = new Locus(1, 1);
//    Locus *s2 = new Locus(1, 2);
//    BulkDatum d1("s1", s1, 10, 100, 2.0);
//    BulkDatum d2("s2", s2, 5, 100, 2.0);
//
//    state->move_datum(child0, &d1, model_params);
//    state->move_datum(child1, &d2, model_params);
//
//    double freq0 = 0.2;
//    double freq1 = 0.1;
//    child0->set_cellular_prev(freq0 + freq1);
//    child0->set_clone_freq(freq0);
//
//    child1->set_cellular_prev(freq1);
//    child1->set_clone_freq(freq1);
//
//    root->set_cellular_prev(1);
//    root->set_clone_freq(1 - freq0 - freq1);
//
//    double xi = child0->compute_expectation_of_xi(model_params);
//    double expected_val = 0.2339236;
//    BOOST_TEST( abs(xi - expected_val) < 1e-3 );
//    cout << "done." << endl;
//}
//
//BOOST_AUTO_TEST_CASE( test_dynamic_program_c )
//{
//    // generate TSSB object, create a small tree, insert data
//    double seq_err = 0.05;
//    double b_rate = 0, d_rate = 0;
//    ModelParams model_params(1, 0.5, 0.8, 0.01, seq_err);
//    model_params.set_birth_rate(b_rate);
//    model_params.set_death_rate(d_rate);
//    model_params.recompute_transition_matrix();
//
//    CloneTreeNode *root = CloneTreeNode::create_root_node();
//    auto state = TSSBState<BulkDatum,SingleCellDNA,CloneTreeNodeParam>::construct_trivial_state(root);
//
//    CloneTreeNode *child0 = (CloneTreeNode *)root->spawn_child(0.5);
//    CloneTreeNode *child1 = (CloneTreeNode *)child0->spawn_child(0.5);
//    CloneTreeNode *child2 = (CloneTreeNode *)child0->spawn_child(0.5);
//
//    // assign datum to each of the 3 nodes
//    Locus *s1 = new Locus(1, 1);
//    Locus *s2 = new Locus(1, 2);
//    Locus *s3 = new Locus(1, 3);
//    BulkDatum d1("s1", s1, 30, 100, 2.0);
//    BulkDatum d2("s2", s2, 5, 100, 2.0);
//    BulkDatum d3("s3", s3, 15, 100, 2.0);
//
//    state->move_datum(child0, &d1, model_params);
//    state->move_datum(child1, &d2, model_params);
//    state->move_datum(child2, &d3, model_params);
//
//    double freq0 = 0.2;
//    double freq1 = 0.1;
//    double freq2 = 0.3;
//    child0->set_cellular_prev(freq0 + freq1 + freq2);
//    child0->set_clone_freq(freq0);
//
//    child1->set_cellular_prev(freq1);
//    child1->set_clone_freq(freq1);
//
//    child2->set_cellular_prev(freq2);
//    child2->set_clone_freq(freq2);
//
//    root->set_cellular_prev(1);
//    root->set_clone_freq(1 - child0->get_node_parameter().get_cellular_prev());
//
//    double xi = child0->compute_expectation_of_xi(model_params);
//    double expected_val = 0.42;
//    BOOST_TEST( abs(xi - expected_val) < 1e-3 );
//    
//    b_rate = 1, d_rate = 1;
//    model_params.set_birth_rate(b_rate);
//    model_params.set_death_rate(d_rate);
//    model_params.recompute_transition_matrix();
//
//    xi = child0->compute_expectation_of_xi(model_params);
//    expected_val = 0.314858;
//    BOOST_TEST( abs(xi - expected_val) < 1e-3 );
//
//    cout << "done." << endl;
//}
//
//BOOST_AUTO_TEST_CASE( test_dynamic_program_d )
//{
//    // generate TSSB object, create a small tree, insert data
//    double seq_err = 0.05;
//    double b_rate = 0, d_rate = 0;
//    ModelParams model_params(1, 0.5, 0.8, 0.01, seq_err);
//    model_params.set_birth_rate(b_rate);
//    model_params.set_death_rate(d_rate);
//    model_params.recompute_transition_matrix();
//
//    CloneTreeNode *root = CloneTreeNode::create_root_node();
//    auto state = TSSBState<BulkDatum,SingleCellDNA,CloneTreeNodeParam>::construct_trivial_state(root);
//
//    CloneTreeNode *child0 = (CloneTreeNode *)root->spawn_child(0.5);
//    CloneTreeNode *child1 = (CloneTreeNode *)child0->spawn_child(0.5);
//    CloneTreeNode *child2 = (CloneTreeNode *)child1->spawn_child(0.5);
//
//    double freq0 = 0.2;
//    double freq1 = 0.1;
//    double freq2 = 0.3;
//    child0->set_cellular_prev(freq0 + freq1 + freq2);
//    child0->set_clone_freq(freq0);
//
//    child1->set_cellular_prev(freq1 + freq2);
//    child1->set_clone_freq(freq1);
//
//    child2->set_cellular_prev(freq2);
//    child2->set_clone_freq(freq2);
//
//    root->set_cellular_prev(1);
//    root->set_clone_freq(1 - child0->get_node_parameter().get_cellular_prev());
//
//    double xi = child0->compute_expectation_of_xi(model_params);
//    double expected_val = 0.42;
//    BOOST_TEST( abs(xi - expected_val) < 1e-3 );
//
//    b_rate = 1, d_rate = 1;
//    model_params.set_birth_rate(b_rate);
//    model_params.set_death_rate(d_rate);
//    model_params.recompute_transition_matrix();
//
//    xi = child0->compute_expectation_of_xi(model_params);
//    expected_val = 0.2971249;
//    BOOST_TEST( abs(xi - expected_val) < 1e-3 );
//
//    cout << "done." << endl;
//}
//
//BOOST_AUTO_TEST_CASE( test_dynamic_program_e )
//{
//    // generate TSSB object, create a small tree, insert data
//    double seq_err = 0.05;
//    double b_rate = 0, d_rate = 0;
//    ModelParams model_params(1, 0.5, 0.8, 0.01, seq_err);
//    model_params.set_birth_rate(b_rate);
//    model_params.set_death_rate(d_rate);
//    model_params.recompute_transition_matrix();
//    
//    CloneTreeNode *root = CloneTreeNode::create_root_node();
//    auto state = TSSBState<BulkDatum,SingleCellDNA,CloneTreeNodeParam>::construct_trivial_state(root);
//    
//    CloneTreeNode *child0 = (CloneTreeNode *)root->spawn_child(0.5);
//    CloneTreeNode *child1 = (CloneTreeNode *)child0->spawn_child(0.5);
//    CloneTreeNode *child2 = (CloneTreeNode *)child1->spawn_child(0.5);
//    CloneTreeNode *child3 = (CloneTreeNode *)child1->spawn_child(0.5);
//    
//    double freq0 = 0.2;
//    double freq1 = 0.1;
//    double freq2 = 0.3;
//    double freq3 = 0.1;
//    child0->set_cellular_prev(freq0 + freq1 + freq2 + freq3);
//    child0->set_clone_freq(freq0);
//
//    child1->set_cellular_prev(freq1 + freq2 + freq3);
//    child1->set_clone_freq(freq1);
//
//    child2->set_cellular_prev(freq2);
//    child2->set_clone_freq(freq2);
//
//    child3->set_cellular_prev(freq3);
//    child3->set_clone_freq(freq3);
//
//    root->set_cellular_prev(1);
//    root->set_clone_freq(1 - child0->get_node_parameter().get_cellular_prev());
//
//    double xi = child0->compute_expectation_of_xi(model_params);
//    double expected_val = 0.2333333 + 0.2333333 + 0.015;
//    BOOST_TEST( abs(xi - expected_val) < 1e-3 );
//
//    b_rate = 1, d_rate = 1;
//    model_params.set_birth_rate(b_rate);
//    model_params.set_death_rate(d_rate);
//    model_params.recompute_transition_matrix();
//
//    xi = child0->compute_expectation_of_xi(model_params);
//    expected_val = 0.3322373;
//    BOOST_TEST( abs(xi - expected_val) < 1e-3 );
//
//    cout << "done." << endl;
//}
//
//BOOST_AUTO_TEST_CASE( test_dynamic_program_f )
//{
//    // generate TSSB object, create a small tree, insert data
//    double seq_err = 0.05;
//    double b_rate = 0, d_rate = 0;
//    ModelParams model_params(1, 0.5, 0.8, 0.01, seq_err);
//    model_params.set_birth_rate(b_rate);
//    model_params.set_death_rate(d_rate);
//    model_params.recompute_transition_matrix();
//    
//    CloneTreeNode *root = CloneTreeNode::create_root_node();
//    auto state = TSSBState<BulkDatum,SingleCellDNA,CloneTreeNodeParam>::construct_trivial_state(root);
//    
//    CloneTreeNode *u = (CloneTreeNode *)root->spawn_child(0.5);
//    CloneTreeNode *v = (CloneTreeNode *)u->spawn_child(0.5);
//    CloneTreeNode *w = (CloneTreeNode *)u->spawn_child(0.5);
//    CloneTreeNode *y = (CloneTreeNode *)v->spawn_child(0.5);
//    CloneTreeNode *z = (CloneTreeNode *)w->spawn_child(0.5);
//    
//    double freq_u = 0.2;
//    double freq_v = 0.1;
//    double freq_w = 0.3;
//    double freq_y = 0.05;
//    double freq_z = 0.1;
//    u->set_cellular_prev(freq_u + freq_v + freq_w + freq_y + freq_z);
//    u->set_clone_freq(freq_u);
//
//    v->set_cellular_prev(freq_v + freq_y);
//    v->set_clone_freq(freq_v);
//
//    w->set_cellular_prev(freq_w + freq_z);
//    w->set_clone_freq(freq_w);
//
//    y->set_cellular_prev(freq_y);
//    y->set_clone_freq(freq_y);
//
//    z->set_cellular_prev(freq_z);
//    z->set_clone_freq(freq_z);
//
//    root->set_cellular_prev(1);
//    root->set_clone_freq(1 - u->get_node_parameter().get_cellular_prev());
//
//    double xi = u->compute_expectation_of_xi(model_params);
//    double expected_val = 0.5125;
//    BOOST_TEST( abs(xi - expected_val) < 1e-3 );
//
//    b_rate = 1, d_rate = 1;
//    model_params.set_birth_rate(b_rate);
//    model_params.set_death_rate(d_rate);
//    model_params.recompute_transition_matrix();
//    
//    xi = u->compute_expectation_of_xi(model_params);
//    expected_val = 0.3675266;
//    BOOST_TEST( abs(xi - expected_val) < 1e-3 );
//    
//    cout << "done." << endl;
//}

