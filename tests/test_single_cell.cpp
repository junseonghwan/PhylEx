//
//  test_single_cell.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-08-08.
//
#include <boost/test/unit_test.hpp>

#include "clone_node.hpp"
#include "data_util.hpp"
#include "model_params.hpp"
#include "node.hpp"
#include "numerical_utils.hpp"
#include "sampling_utils.hpp"
#include "utils.hpp"

#include "simul_data.hpp"

BOOST_AUTO_TEST_CASE( test_sc_likelihood_cache )
{
    SimulationConfig simul_config;

    gsl_rng *random = generate_random_object(3);
    ModelParams model_params(3, 1, 0.8, 0.01, 0.01);

    simul_config.bulk_mean_depth = 1000;
    simul_config.max_cn = 6;
    simul_config.sc_mean_depth = 20;
    simul_config.n_cells = 10;
    simul_config.prob_expression = 0.5;
    simul_config.zero_inflation_prob = 0.5;
    simul_config.birth_rate = 0;
    simul_config.death_rate = 0;
    simul_config.balance_ref_var_alleles = true;
    simul_config.read_length = 150;

    CloneTreeNode *root = CloneTreeNode::create_root_node();
    root->set_nu_stick(0);
    root->sample_node_parameters(random, model_params, 0);
    auto true_tree = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::construct_trivial_state(root, bulk_likelihood, sc_likelihood);

    size_t n_data = 20;
    vector<BulkDatum *> data;
    unordered_set<Locus> somatic_loci;
    for (size_t i = 0; i < n_data; i++)
    {
        string chr = convert_chr_to_string(gsl_rng_uniform_int(random, 23) + 1);
        size_t pos = gsl_rng_uniform_int(random, 1000000000);
        Locus *loci = new Locus(chr, pos);
        BulkDatum *datum = new BulkDatum("s" + to_string(i+1), *loci);
        true_tree->insert_datum(random, datum, model_params);
        data.push_back(datum);
        somatic_loci.insert(*loci);
    }
    
    unordered_map<Node<BulkDatum,CloneTreeNodeParam> *, vector<pair<size_t, size_t> > > cn_profile;
    generate_bulk_data(random, simul_config, data, model_params, root, cn_profile);

    cout << true_tree->print() << endl;

    vector<SingleCellData *> *sc_data = new vector<SingleCellData *>();
    //generate_scDNA_data(random, root, data, cn_profile, model_params, simul_config, *sc_data);
    generate_scRNA_data(random, root, somatic_loci, cn_profile, model_params, simul_config, *sc_data);
    true_tree->set_sc_data(sc_data, model_params);

    double sc_log_lik = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::compute_log_likelihood_sc(root, data, *sc_data, sc_likelihood, model_params);
    double sc_log_lik_cache = true_tree->compute_log_likelihood_sc_cached(model_params);
    cout << "At init: " << sc_log_lik << ", " << sc_log_lik_cache << endl;
    BOOST_TEST( abs(sc_log_lik - sc_log_lik_cache) < 1e-3 );

    BulkDatum *datum0 = data[0];
    Node<BulkDatum,CloneTreeNodeParam> *node0 = true_tree->get_node(datum0);
    BulkDatum *datum1 = data[2];
    true_tree->move_datum(node0, datum1, model_params);
    sc_log_lik = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::compute_log_likelihood_sc(root, data, *sc_data, sc_likelihood, model_params);
    sc_log_lik_cache = true_tree->compute_log_likelihood_sc_cached(model_params);
    cout << "After move: " << sc_log_lik << ", " << sc_log_lik_cache << endl;
    BOOST_TEST( abs(sc_log_lik - sc_log_lik_cache) < 1e-3 );

    for (size_t i = 0; i < 20; i++) {
        true_tree->resample_data_assignment(random, model_params);
        sc_log_lik = TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::compute_log_likelihood_sc(root, data, *sc_data, sc_likelihood, model_params);
        sc_log_lik_cache = true_tree->compute_log_likelihood_sc_cached(model_params);
        cull(true_tree->get_root()); // make sure that culling does not introduce any problem to the cache contents
        cout << "After resampling assignment: " << sc_log_lik << ", " << sc_log_lik_cache << endl;
        BOOST_TEST( abs(sc_log_lik - sc_log_lik_cache) < 1e-3 );
    }
}
