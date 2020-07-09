#define BOOST_TEST_MODULE tssb test
#include <boost/test/unit_test.hpp>

#include <math.h>
#include <iostream>

#include "clone_node.hpp"
#include "single_cell.hpp"
#include "tssb_state.hpp"

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
    vector<size_t> total_cn(1, 2);
    vector<size_t> var_read_count(1, 10);
    vector<size_t> total_read_count(1, 20);
    
    auto root = CloneTreeNode::create_root_node(1);
    auto child = (CloneTreeNode *)root->spawn_child(0.5);
    child->set_cellular_prev(0, 0.5);
    child->set_clone_freq(0, 0.5);

    BulkDatum datum1("s1", "chr", 0);
    datum1.AddRegionData(0, 0, 1, 1);
    double realized = BulkLogLikWithTotalCopyNumber(0, child, &datum1, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(realized == 0.0);

    BulkDatum datum2("s2", "chr", 0, var_read_count, total_read_count, total_cn);
    realized = BulkLogLikWithTotalCopyNumber(0, child, &datum2, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -2.726685148) < 1e-6);
    
    // Change total_cn to 8.
    total_cn[0] = 8;
    BulkDatum datum3("s3", "chr", 0, var_read_count, total_read_count, total_cn);
    realized = BulkLogLikWithTotalCopyNumber(0, child, &datum3, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -3.602852362) < 1e-6);
    
    total_cn[0] = 0;
    BulkDatum datum4("s4", "chr", 0, var_read_count, total_read_count, total_cn);
    realized = BulkLogLikWithTotalCopyNumber(0, child, &datum4, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -56.96076648) < 1e-6);
}

BOOST_AUTO_TEST_CASE( TestBulkLogLikWithGenotype )
{
    vector<size_t> major_cn(1, 3);
    vector<size_t> minor_cn(1, 1);
    vector<size_t> var_read_count(1, 10);
    vector<size_t> total_read_count(1, 20);

    auto root = CloneTreeNode::create_root_node(1);
    auto child = (CloneTreeNode *)root->spawn_child(0.5);
    child->set_cellular_prev(0, 0.5);
    child->set_clone_freq(0, 0.5);
    
    BulkDatum datum1("s1", "chr", 0, var_read_count, total_read_count,
                     major_cn, minor_cn);
    double realized = BulkLogLikWithGenotype(0, child, &datum1, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -3.658831964) < 1e-6);

    major_cn[0] = 2;
    minor_cn[0] = 0;
    BulkDatum datum2("s2", "chr", 0, var_read_count, total_read_count,
                     major_cn, minor_cn);
    realized = BulkLogLikWithGenotype(0, child, &datum2, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -2.373800323) < 1e-6);

    var_read_count[0] = 0;
    total_read_count[0] = 0;
    BulkDatum datum3("s3", "chr", 0, var_read_count, total_read_count,
                     major_cn, minor_cn);
    realized = BulkLogLikWithGenotype(0, child, &datum3, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - 0) < 1e-6);
}

BOOST_AUTO_TEST_CASE( TestScLikelihood )
{
    double alpha = 3.0;
    double beta = 12.0;
    double dropout_prob = 0.8;
    
    BulkDatum bulk("", "chr", 0);
    bulk.SetLocuHyperParameters(alpha, beta, dropout_prob);
    
    SingleCellData c0("", 1);
    c0.InsertDatum(0, 10, 20);
    
    double realized = ScLikelihood(0, &bulk, &c0, true, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -5.488633383) < 1e-6);
    
    realized = ScLikelihood(0, &bulk, &c0, false, model_params);
    std::cout << realized << std::endl;
    BOOST_TEST(fabs(realized - -9.210441917) < 1e-6);
}

BOOST_AUTO_TEST_CASE( TestScCache )
{
    gsl_rng *random = generate_random_object(3);
    ModelParams model_params(3, 1, 0.8, 0.01);
    size_t region_count = 1;
    size_t single_cell_count = 10;
    size_t bulk_data_count = 20;
    size_t max_read_count = 100;

    // Generate some bulk data.
    vector<BulkDatum *> bulk_data;
    string mutation_id;
    for (size_t i = 0; i < bulk_data_count; i++)
    {
        vector<size_t> total_reads(1, gsl_rng_uniform_int(random, 1000));
        vector<size_t> var_reads(1, gsl_rng_uniform_int(random, total_reads[0]));
        vector<size_t> major_cns(1, 1);
        vector<size_t> minor_cns(1, 1);

        BulkDatum *datum = new BulkDatum(mutation_id, "chr", i,
                                         var_reads, total_reads,
                                         major_cns, minor_cns);
        bulk_data.push_back(datum);
    }
    
    // Generate single cell data.
    vector<SingleCellData *> sc_data;
    for (size_t c = 0; c < single_cell_count; c++) {
        vector<size_t> var_reads(bulk_data_count), total_reads(bulk_data_count);
        for (size_t loci_idx = 0; loci_idx < bulk_data_count; loci_idx++) {
            size_t total_read_count = gsl_rng_uniform_int(random, max_read_count);
            size_t var_read_count = 0;
            if (total_read_count > 0) {
                var_read_count = gsl_rng_uniform_int(random, total_read_count);
            }
            var_reads[loci_idx] = var_read_count;
            total_reads[loci_idx] = total_read_count;
        }
        SingleCellData *sc = new SingleCellData("c" + to_string(c), var_reads, total_reads);
        sc_data.push_back(sc);
    }

    CloneTreeNode *root = CloneTreeNode::create_root_node(region_count);
    TSSBState tree(random, root, model_params,
                   BulkLogLikWithGenotype, ScLikelihood,
                   &bulk_data, &sc_data);
    cout << tree.print() << endl;

    double sc_log_lik = TSSBState::compute_log_likelihood_sc(root,
                                                             bulk_data,
                                                             sc_data,
                                                             ScLikelihood,
                                                             model_params);
    double sc_log_lik_cache = tree.compute_log_likelihood_sc_cached(model_params);
    cout << "At init: " << sc_log_lik << ", " << sc_log_lik_cache << endl;
    BOOST_TEST( abs(sc_log_lik - sc_log_lik_cache) < 1e-3 );
    
    // Resample assignment of bulk data.
    tree.resample_data_assignment(random, model_params);
    cout << tree.print() << endl;
    sc_log_lik = TSSBState::compute_log_likelihood_sc(root,
                                                      bulk_data,
                                                      sc_data,
                                                      ScLikelihood,
                                                      model_params);
    sc_log_lik_cache = tree.compute_log_likelihood_sc_cached(model_params);
    cout << "After move: " << sc_log_lik << ", " << sc_log_lik_cache << endl;
    BOOST_TEST( abs(sc_log_lik - sc_log_lik_cache) < 1e-3 );
}
