//
//  interface.cpp
//  run
//
//  Created by Seong-Hwan Jun on 2020-06-02.
//

#include "interface.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iterator>
#include <iostream>
#include <stdio.h>
#include <string>
#include <unordered_map>

#include <gsl/gsl_statistics.h>

#include "numerical_utils.hpp"
#include "utils.hpp"

void write_best_trees(string output_path,
                      const vector<BulkDatum *> &bulk,
                      vector<pair<double, shared_ptr<CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> > > > &best_states)
{
    size_t n_trees = best_states.size();
    for (size_t i = 0; i < best_states.size(); i++)
    {
        string path = output_path + "/tree" + to_string(i);
        size_t idx = n_trees-i-1;
        write_tree(path,
                   bulk,
                   *best_states[idx].second.get());
        WriteLogLikToFile(path + "/log_lik.txt", best_states[idx].first);
    }
}

bool comp_states(const pair<double, shared_ptr<CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> > > &s1,
                 const pair<double, shared_ptr<CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> > > &s2)
{
    return (s1.first < s2.first);
}

TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> *RunSliceSampler(
    const gsl_rng *random,
    ModelParams &params,
    Config &config,
    vector<BulkDatum *> *bulk_data,
    vector<SingleCellData *> *sc_data,
    CopyNumberInputType cn_input_type)
{
    CloneTreeNode *root = CloneTreeNode::create_root_node();
    root->sample_node_parameters(random, params, 0);

    TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> *tree;
    switch (cn_input_type) {
        case GENOTYPE:
            tree = new TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>(random,
                                                                              root,
                                                                              params,
                                                                              BulkLogLikWithGenotype,
                                                                              ScLikelihood,
                                                                              bulk_data,
                                                                              sc_data);
            break;
        case TOTAL_CN:
            tree = new TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>(random,
                                                                              root,
                                                                              params,
                                                                              BulkLogLikWithTotalCopyNumber,
                                                                              ScLikelihood,
                                                                              bulk_data,
                                                                              sc_data);
            break;
        case TOTAL_CN_PROFILE:
            tree = new TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>(random,
                                                                              root,
                                                                              params,
                                                                              BulkLogLikWithCopyNumberProfile,
                                                                              ScLikelihood,
                                                                              bulk_data,
                                                                              sc_data);
            break;
        case UNDETERMINED:
            exit(-1);
    }

    cull(tree->get_root());
    cout << "Current tree:" << endl;
    cout << tree->print() << endl;
    
    vector<double> alpha0, lambda, gamma;
    
    // Keep track of top 5 trees.
    size_t state_count = 5;
    vector<pair<double, shared_ptr<CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> > > > joint_best;
    vector<pair<double, shared_ptr<CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> > > > states;
    
    size_t n_trees = 0;
    size_t n_iter = config.burn_in + config.n_mcmc_iter;
    auto start = chrono::steady_clock::now();
    for (size_t iter = 0; iter < n_iter; iter++) {

        cout << "Iter: " << iter << endl;

        // 1. resample assignment
        tree->resample_data_assignment(random, params);

        cull(tree->get_root());
        tree->clear_cache();
        
        // 2. update params
        double ar = sample_params_dirichlet(random, config.n_mh_iter, *tree, params);
        if (iter < config.burn_in) {
            if (ar < 0.08 && params.get_dir_conc_mult_factor() < 10000) {
                params.set_dir_conc_mult_factor(params.get_dir_conc_mult_factor() * 2);
            } else if (ar > 0.5 && ar < 0.99) {
                params.set_dir_conc_mult_factor(params.get_dir_conc_mult_factor() / 2);
            }
        }
        
        double log_lik_post = tree->get_log_lik_bulk();
        double log_lik_prior_assignment = tree->get_log_prior_assignment(tree->get_root());
        double log_lik_sc = tree->get_log_lik_sc();
        double log_lik = log_lik_post + log_lik_prior_assignment + log_lik_sc;
        cout << "joint log lik: " << log_lik << endl;
        cout << "log lik_bulk: " << log_lik_post << endl;
        cout << "log lik_prior: " << log_lik_prior_assignment << endl;
        cout << "log lik_sc: " << log_lik_sc << endl;
        cout << "alpha0: " << params.get_alpha0() << endl;
        cout << "lambda: " << params.get_lambda() << endl;
        cout << "gamma: " << params.get_gamma() << endl;
        
        // 3. update sticks
        tree->update_sticks(random, params);
        
        // 4. resample hyper parameters
        TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::update_hyper_params(random, config.n_mh_iter, *tree, params);
        alpha0.push_back(params.get_alpha0());
        lambda.push_back(params.get_lambda());
        gamma.push_back(params.get_gamma());
        
        if (joint_best.size() < state_count) {
            auto compact_state = shared_ptr<CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam > >(new CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam >(*tree));
            joint_best.push_back(make_pair(log_lik, compact_state));
        } else {
            auto compact_state = shared_ptr<CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam > >(new CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam >(*tree));
            if (log_lik > joint_best[0].first) {
                joint_best[0] = make_pair(log_lik, compact_state);
                sort(joint_best.begin(), joint_best.end(), comp_states);
            }
        }
        
        if ((iter % config.thinning) == 0 || (iter + 1) == n_iter) {
            // store current tree
            auto compact_state = shared_ptr<CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> >(new  CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>(*tree));
            states.push_back(make_pair(log_lik, compact_state));
        }

        if ((iter % config.output_interval) == 0) {
            write_best_trees(config.output_path + "/joint", *bulk_data, joint_best);
            for (size_t i = n_trees; i < states.size(); i++) {
                write_tree(config.output_path + "/states/tree" + to_string(n_trees) + "/", *bulk_data, *states[i].second.get());
                WriteLogLikToFile(config.output_path + "/states/tree" + to_string(n_trees) + "/log_lik.txt", states[i].first);
                n_trees++;
            }
        }
    }
    auto end = chrono::steady_clock::now();
    cout << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms." << endl;
    
    // flush the states
    write_best_trees(config.output_path + "/joint", *bulk_data, joint_best);
    for (size_t i = n_trees; i < states.size(); i++) {
        write_tree(config.output_path + "/states/tree" + to_string(n_trees) + "/", *bulk_data, *states[i].second.get());
        WriteLogLikToFile(config.output_path + "/states/tree" + to_string(n_trees) + "/log_lik.txt", states[i].first);
        n_trees++;
    }
    
    return tree;
}

Config parse_config_file(string config_file_path)
{
    string line;
    ifstream config_file(config_file_path);
    if (!config_file.is_open())
    {
        cerr << "Could not open the file: " << config_file_path << endl;
        exit(-1);
    }
    
    Config config;
    vector<string> results;
    vector<double> dat;
    
    while ( getline (config_file, line) )
    {
        boost::split(results, line, boost::is_any_of(":"));
        boost::algorithm::trim(results[1]);
        if (results[0] == "seed") {
            config.seed = stoul(results[1]);
        } else if (results[0] == "bulk_data_path") {
            config.bulk_file = results[1];
        } else if (results[0] == "sc_rna_data_path") {
            config.scRNA_file = results[1];
        } else if (results[0] == "hyperparams_path") {
            config.sc_hyperparam_file = results[1];
        } else if (results[0] == "cn_prior_path") {
            config.cn_prior_path = results[1];
        } else if (results[0] == "output_path") {
            config.output_path = results[1];
        } else if (results[0] == "n_mcmc_iter") {
            config.n_mcmc_iter = stoul(results[1]);
        } else if (results[0] == "n_mh_iter") {
            config.n_mh_iter = stoul(results[1]);
        } else if (results[0] == "thinning") {
            config.thinning = stoul(results[1]);
        } else if (results[0] == "burn_in") {
            config.burn_in = stoul(results[1]);
        } else if (results[0] == "alpha0_max") {
            config.alpha0_max = stod(results[1]);
        } else if (results[0] == "lambda_max") {
            config.lambda_max = stod(results[1]);
        } else if (results[0] == "gamma_max") {
            config.gamma_max = stod(results[1]);
        } else if (results[0] == "seq_err") {
            config.seq_err = stod(results[1]);
        } else if (results[0] == "var_cp_prob") {
            config.var_cp_prob_ = stod(results[1]);
        } else if (results[0] == "sc_dropout_alpha0") {
            config.sc_dropout_alpha0 = stod(results[1]);
        } else if (results[0] == "sc_dropout_beta0") {
            config.sc_dropout_beta0 = stod(results[1]);
        } else {
            cerr << "Unknown option: " << results[0] << endl;
        }
    }
    config_file.close();
    return config;
}

void Run(string config_file)
{
    // Parse the configuration file
    Config config = parse_config_file(config_file);
    cout << "Finished reading config file:" << config_file << endl;
    
    gsl_rng *random = generate_random_object(config.seed);
    ModelParams params = ModelParams::RandomInit(random,
                                                 config.alpha0_max,
                                                 config.lambda_max,
                                                 config.gamma_max,
                                                 config.seq_err);
    params.set_sc_dropout_alpha0(config.sc_dropout_alpha0);
    params.set_sc_dropout_beta0(config.sc_dropout_beta0);
    params.set_var_cp_prob(config.var_cp_prob_);

    vector<BulkDatum *> bulk_data;
    unordered_map<string, Locus *> somatic_loci;
    CopyNumberInputType cn_input_type = ReadBulkData(config.bulk_file, bulk_data, somatic_loci);
    if (config.cn_prior_path != "") {
        ReadCnPrior(config.cn_prior_path, bulk_data);
        cn_input_type = CopyNumberInputType::TOTAL_CN_PROFILE;
    }
    if (cn_input_type == CopyNumberInputType::UNDETERMINED) {
        cerr << "Error: Copy number input type is undetermined.\n";
        exit(-1);
    }
    vector<SingleCellData *> sc_data;
    if (config.scRNA_file != "") {
        ReadScRnaData(config.scRNA_file, somatic_loci, sc_data);
        ReadScRnaHyperparams(config.sc_hyperparam_file, somatic_loci);
    }

    RunSliceSampler(random, params, config, &bulk_data, &sc_data, cn_input_type);
}
