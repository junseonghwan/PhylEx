#define BOOST_TEST_MODULE tssb test
#include <boost/test/unit_test.hpp>

#include <math.h>
#include <iostream>

#include "clone_node.hpp"
#include "single_cell.hpp"

const double ALPHA0 = 1.5;
const double GAMMA = 0.5;
const double LAMBDA = 1.1;
const double SEQ_ERR = 0.001;
const double VAR_CP_PROB = 0.5;

ModelParams model_params(ALPHA0, GAMMA, LAMBDA, SEQ_ERR);

void InitializeTestSetup() {
    model_params.set_var_cp_prob(VAR_CP_PROB);
    model_params.set_sc_dropout_alpha0(0.01);
    model_params.set_sc_dropout_beta0(0.01);
}

BOOST_AUTO_TEST_CASE( TestBulkLogLikWithTotalCopyNumber )
{
    size_t total_cn = 2;
    size_t var_read_count = 10;
    size_t total_read_count = 20;
    Locus locus("", "", 0, "");
    
    auto root = CloneTreeNode::create_root_node();
    auto child = (CloneTreeNode *)root->spawn_child(0.5);
    child->set_cellular_prev(0.5);
    child->set_clone_freq(0.5);

    BulkDatum datum1("s1", locus, 0, 0, total_cn);
    double realized = BulkLogLikWithTotalCopyNumber(child, &datum1, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(realized == 0.0);

    BulkDatum datum2("s2", locus, var_read_count, total_read_count, total_cn);
    realized = BulkLogLikWithTotalCopyNumber(child, &datum2, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -2.726685148) < 1e-6);
    
    // Change total_cn to 8.
    total_cn = 8;
    datum2.SetTotalCopyNumber(total_cn);
    realized = BulkLogLikWithTotalCopyNumber(child, &datum2, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -3.602852362) < 1e-6);
    
    total_cn = 0;
    BulkDatum datum3("s1", locus, var_read_count, total_read_count, total_cn);
    realized = BulkLogLikWithTotalCopyNumber(child, &datum3, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -56.96076648) < 1e-6);
}

BOOST_AUTO_TEST_CASE( TestBulkLogLikWithGenotype )
{
    size_t major_cn = 3;
    size_t minor_cn = 1;
    size_t var_read_count = 10;
    size_t total_read_count = 20;
    Locus locus("", "", 0, "");

    auto root = CloneTreeNode::create_root_node();
    auto child = (CloneTreeNode *)root->spawn_child(0.5);
    child->set_cellular_prev(0.5);
    child->set_clone_freq(0.5);
    
    BulkDatum datum1("s1", locus, var_read_count, total_read_count,
                     major_cn, minor_cn);
    double realized = BulkLogLikWithGenotype(child, &datum1, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -3.658831964) < 1e-6);
    
    major_cn = 2;
    minor_cn = 0;
    datum1.SetGenotype(make_pair(major_cn, minor_cn));
    realized = BulkLogLikWithGenotype(child, &datum1, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -2.373800323) < 1e-6);
    
    BulkDatum datum2("s2", locus, 0, 0,
                     major_cn, minor_cn);
    realized = BulkLogLikWithGenotype(child, &datum2, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - 0) < 1e-6);
}

BOOST_AUTO_TEST_CASE( TestScLikelihood )
{
    double alpha = 3.0;
    double beta = 12.0;
    double dropout_prob = 0.8;
    Locus locus("s1", "", 0, "");
    locus.set_alpha(alpha);
    locus.set_beta(beta);
    locus.set_dropout_prob(dropout_prob);
    
    LocusDatum locus_datum(20, 10);
    unordered_map<Locus, LocusDatum*> reads;
    reads[locus] = &locus_datum;
    
    SingleCellData c0("", reads);
    BulkDatum bulk("", locus);
    double realized = ScLikelihood(&bulk, &c0, true, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -5.488633383) < 1e-6);
    
    realized = ScLikelihood(&bulk, &c0, false, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -9.210441917) < 1e-6);
}
