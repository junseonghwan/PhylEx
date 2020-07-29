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

#include "data_util.hpp"
#include "numerical_utils.hpp"
#include "utils.hpp"

void write_best_trees(string output_path,
                      const vector<BulkDatum *> &bulk,
                      vector<pair<double, shared_ptr<CompactTSSBState > > > &best_states)
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

bool comp_states(const pair<double, shared_ptr<CompactTSSBState > > &s1,
                 const pair<double, shared_ptr<CompactTSSBState > > &s2)
{
    return (s1.first < s2.first);
}

TSSBState *Interface::RunSliceSampler(const gsl_rng *random,
                                      ModelParams &params)
{
    size_t region_count = bulk_data_[0]->GetRegionCount();
    CloneTreeNode *root = CloneTreeNode::create_root_node(region_count);
    root->sample_node_parameters(random, params, 0);

    TSSBState *tree;
    switch (cn_input_type_) {
        case CopyNumberInputType::GENOTYPE:
            tree = new TSSBState(random,
                                 root,
                                 params,
                                 BulkLogLikWithGenotype,
                                 ScLikelihood,
                                 &bulk_data_,
                                 &sc_data_);
            break;
        case CopyNumberInputType::TOTAL_CN:
            tree = new TSSBState(random,
                                 root,
                                 params,
                                 BulkLogLikWithTotalCopyNumber,
                                 ScLikelihood,
                                 &bulk_data_,
                                 &sc_data_);
            break;
        case CopyNumberInputType::TOTAL_CN_PROFILE:
            tree = new TSSBState(random,
                                 root,
                                 params,
                                 BulkLogLikWithCopyNumberProfile,
                                 ScLikelihood,
                                 &bulk_data_,
                                 &sc_data_);
            break;
        case CopyNumberInputType::UNDETERMINED:
            exit(-1);
    }

    cull(tree->get_root());
    cout << "Current tree:" << endl;
    cout << tree->print() << endl;
    
    vector<double> alpha0, lambda, gamma;
    
    // Keep track of top 5 trees.
    size_t state_count = 5;
    vector<pair<double, shared_ptr<CompactTSSBState > > > joint_best;
    vector<pair<double, shared_ptr<CompactTSSBState > > > states;
    
    size_t n_trees = 0;
    size_t n_iter = config_.burn_in + config_.n_mcmc_iter;
    auto start = chrono::steady_clock::now();
    for (size_t iter = 0; iter < n_iter; iter++) {

        cout << "Iter: " << iter << endl;

        // 1. resample assignment
        tree->resample_data_assignment(random, params);
        
        cull(tree->get_root());
        cout << tree->print() << "\n";
        //tree->clear_cache();
        
        // 2. update params
        double ar = sample_params_dirichlet(random, config_.n_mh_iter, *tree, params);
        if (iter < config_.burn_in) {
            if (ar < 0.08 && params.get_dir_conc_mult_factor() < 10000) {
                params.set_dir_conc_mult_factor(params.get_dir_conc_mult_factor() * 2);
            } else if (ar > 0.5 && ar < 0.99) {
                params.set_dir_conc_mult_factor(params.get_dir_conc_mult_factor() / 2);
            }
        }
        
        double log_lik_post = tree->compute_log_likelihood_bulk(params);
        double log_lik_prior_assignment = tree->get_log_prior_assignment(tree->get_root());
        double log_lik_sc = tree->compute_log_likelihood_sc();
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
        TSSBState::update_hyper_params(random, config_.n_mh_iter, *tree, params);
        alpha0.push_back(params.get_alpha0());
        lambda.push_back(params.get_lambda());
        gamma.push_back(params.get_gamma());
        
        if (joint_best.size() < state_count) {
            auto compact_state = shared_ptr<CompactTSSBState >(new CompactTSSBState(*tree));
            joint_best.push_back(make_pair(log_lik, compact_state));
        } else {
            auto compact_state = shared_ptr<CompactTSSBState >(new CompactTSSBState(*tree));
            if (log_lik > joint_best[0].first) {
                joint_best[0] = make_pair(log_lik, compact_state);
                sort(joint_best.begin(), joint_best.end(), comp_states);
            }
        }
        
        if ((iter % config_.thinning) == 0 || (iter + 1) == n_iter) {
            // store current tree
            auto compact_state = shared_ptr<CompactTSSBState >(new  CompactTSSBState(*tree));
            states.push_back(make_pair(log_lik, compact_state));
        }

        if ((iter % config_.output_interval) == 0) {
            write_best_trees(config_.output_path + "/joint", bulk_data_, joint_best);
            for (size_t i = n_trees; i < states.size(); i++) {
                write_tree(config_.output_path + "/states/tree" + to_string(n_trees) + "/", bulk_data_, *states[i].second.get());
                WriteLogLikToFile(config_.output_path + "/states/tree" + to_string(n_trees) + "/log_lik.txt", states[i].first);
                n_trees++;
            }
        }
    }
    auto end = chrono::steady_clock::now();
    cout << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms." << endl;
    
    // flush the states
    write_best_trees(config_.output_path + "/joint", bulk_data_, joint_best);
    for (size_t i = n_trees; i < states.size(); i++) {
        write_tree(config_.output_path + "/states/tree" + to_string(n_trees) + "/", bulk_data_, *states[i].second.get());
        WriteLogLikToFile(config_.output_path + "/states/tree" + to_string(n_trees) + "/log_lik.txt", states[i].first);
        n_trees++;
    }
    
    return tree;
}

void ProcessConfigFile(string config_file_path, Config &config)
{
    string line;
    ifstream config_file(config_file_path);
    if (!config_file.is_open())
    {
        cerr << "Could not open the file: " << config_file_path << endl;
        exit(-1);
    }
    
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
}

vector<size_t> ParseRegionalData(string line) {
    vector<string> result;
    boost::split(result, line, boost::is_any_of(","));
    vector<size_t> data(result.size());
    for (size_t i = 0; i < result.size(); i++) {
        data[i] = stoul(result[i]);
    }
    return data;
}

vector<double> ParseRegionalCNData(string line) {
    vector<string> result;
    vector<double> data;
    boost::split(result, line, boost::is_any_of(","));
    for (size_t i = 0; i < result.size(); i++) {
        data[i] = stod(result[i]);
    }
    return data;
}

void ProcessBulkWithTotalCopyNumberProfile(ifstream &dat_file,
                                           vector<BulkDatum *> &bulk_data,
                                           unordered_map<string, size_t> &mut_id2idx) {
    vector<string> results;
    string line;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of("\t"));
        string mut_id = results[0];
        string chr = results[1];
        size_t pos = stol(results[2]);
        auto n_vars = ParseRegionalData(results[6]);
        auto n_reads = ParseRegionalData(results[7]);
        
        if (mut_id2idx.count(mut_id) > 0) {
            cerr << "Error: " << mut_id << " already exists!" << endl;
            exit(-1);
        }
        
        BulkDatum *datum = new BulkDatum(mut_id, chr, pos, n_vars, n_reads);
        mut_id2idx[mut_id] = bulk_data.size();
        bulk_data.push_back(datum);
    }
}

void ProcessBulkWithTotalCopyNumber(ifstream &dat_file,
                                    vector<BulkDatum *> &bulk_data,
                                    unordered_map<string, size_t> &somatic_loci) {
    vector<string> results;
    string line;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of("\t"), boost::token_compress_on);
        string mut_id = results[0];
        string chr = results[1];
        size_t pos = stol(results[2]);
        auto n_vars = ParseRegionalData(results[5]);
        auto n_reads = ParseRegionalData(results[6]);
        auto total_cn = ParseRegionalData(results[7]);
        
        if (somatic_loci.count(mut_id) > 0) {
            cerr << "Error: " << mut_id << " already exists!" << endl;
            exit(-1);
        }
        
        BulkDatum *datum = new BulkDatum(mut_id, chr, pos,
                                         n_vars, n_reads, total_cn);
        somatic_loci[mut_id] = bulk_data.size();
        bulk_data.push_back(datum);
    }
}

void ProcessBulkWithGenotype(ifstream &dat_file,
                             vector<BulkDatum *> &bulk_data,
                             unordered_map<string, size_t> &somatic_loci)
{
    vector<string> results;
    string line;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of("\t"), boost::token_compress_on);
        string mut_id = results[0];
        string chr = results[1];
        size_t pos = stol(results[2]);
        
        vector<size_t> var_reads = ParseRegionalData(results[5]);
        vector<size_t> total_reads = ParseRegionalData(results[6]);
        
        vector<size_t> major_cns = ParseRegionalData(results[7]);
        vector<size_t> minor_cns = ParseRegionalData(results[8]);
        
        if (somatic_loci.count(mut_id) > 0) {
            cerr << "Error: " << mut_id << " already exists!" << endl;
            exit(-1);
        }
        
        BulkDatum *datum = new BulkDatum(mut_id, chr, pos,
                                         var_reads, total_reads,
                                         major_cns, minor_cns);
        somatic_loci[mut_id] = bulk_data.size();
        bulk_data.push_back(datum);
    }
}

void Interface::ReadBulkData()
{
    string line;
    ifstream dat_file (config_.bulk_file);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << config_.bulk_file << endl;
        exit(-1);
    }
    
    vector<string> results;
    
    // Retrieve the first line to determine the input format.
    getline(dat_file, line);
    boost::split(results, line, boost::is_any_of("\t"));
    
    if (results.size() == BULK_WITH_GENOTYPE_COLUMN_COUNT) {
        ProcessBulkWithGenotype(dat_file, bulk_data_, mut_id2bulk_idx_);
        cn_input_type_ = CopyNumberInputType::GENOTYPE;
    } else if (results.size() == BULK_WITH_TOTAL_CN_COLUMN_COUNT) {
        ProcessBulkWithTotalCopyNumber(dat_file, bulk_data_, mut_id2bulk_idx_);
        cn_input_type_ = CopyNumberInputType::TOTAL_CN;
    } else if (results.size() == BULK_WITH_TOTAL_CN_PRIOR_COLUMN_COUNT) {
        ProcessBulkWithTotalCopyNumberProfile(dat_file, bulk_data_, mut_id2bulk_idx_);
        cn_input_type_ = CopyNumberInputType::UNDETERMINED;
    } else {
        cerr << "Error: invalid bulk input format.\n";
        exit(-1);
    }
    
    dat_file.close();
}

// cn_prior_path points to a file with the following format:
// First column is the mutation ID used to look up BulkDatum from bulk_data.
// Each of the following column is comma separated, one for each region.
// The second column denotes the probability of copy number being 0.
// The third column denotes the probability of copy number being 1 and so on.
void Interface::ReadCnPrior()
{
    ifstream dat_file (config_.cn_prior_path);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << config_.cn_prior_path << endl;
        exit(-1);
    }
    
    unordered_map<string, BulkDatum *> id2bulk;
    for (auto bulk_datum : bulk_data_) {
        id2bulk[bulk_datum->GetId()] = bulk_datum;
    }
    
    vector<string> results;
    string line;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of(" "));
        if (!id2bulk.count(results[0])) {
            continue;
        }
        auto bulk_datum = id2bulk[results[0]];
        
        vector<vector<double> > cn_profile;
        for (size_t i = 1; i < results.size(); i++) {
            auto cn_probs = ParseRegionalCNData(results[i]);
            cn_profile.push_back(cn_probs);
        }
        bulk_datum->SetCopyNumberProbs(cn_profile);
    }
    
    dat_file.close();
}

void Interface::ReadScRnaData()
{
    // first column is header
    // file contains 4 columns:
    // mutation id, cell name, a, d
    string line;
    ifstream dat_file (config_.scRNA_file);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << config_.scRNA_file << endl;
        exit(-1);
    }
    
    vector<string> cell_name_order;
    unordered_map<string, SingleCellData *> dat;
    
    vector<string> results;
    size_t line_idx = 0;
    SingleCellData *sc = 0;
    string mutation_id, cell_name;
    size_t ref_reads, total_reads, var_reads;
    while ( getline (dat_file, line) )
    {
        if (line_idx == 0) {
            line_idx++;
            continue;
        }
        boost::split(results, line, boost::is_any_of("\t"));
        for (size_t i = 0; i < results.size(); i++) {
            boost::algorithm::trim(results[i]);
        }
        mutation_id = results[0];
        cell_name = results[1];
        ref_reads = stol(results[2]);
        total_reads = stol(results[3]);
        var_reads = total_reads - ref_reads;
        
        if (dat.count(cell_name) > 0) {
            sc = dat[cell_name];
        } else {
            sc = new SingleCellData(cell_name, bulk_data_.size());
            dat[cell_name] = sc;
            cell_name_order.push_back(cell_name);
        }
        
        // Find locus
        if (mut_id2bulk_idx_.count(mutation_id) == 0) {
            cerr << "Single cell RNA data contains loci " << mutation_id << " not found in the bulk." << endl;
            exit(-1);
        }
        size_t bulk_idx = mut_id2bulk_idx_[mutation_id];
        sc->InsertDatum(bulk_idx, var_reads, total_reads);
        
        line_idx++;
    }
    dat_file.close();
    
    for (string cell_name : cell_name_order) {
        sc_data_.push_back(dat[cell_name]);
    }
}

void Interface::ReadScRnaHyperparams()
{
    // first column is header
    // file contains 3 columns:
    // mutation id, alpha, beta
    string line;
    ifstream dat_file (config_.sc_hyperparam_file);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << config_.sc_hyperparam_file << endl;
        exit(-1);
    }
    
    vector<string> results;
    size_t line_idx = 0;
    string mutation_id;
    double alpha, beta, delta0;
    while ( getline (dat_file, line) )
    {
        if (line_idx == 0) {
            line_idx++;
            continue;
        }
        boost::split(results, line, boost::is_any_of("\t"));
        for (size_t i = 0; i < results.size(); i++) {
            boost::algorithm::trim(results[i]);
        }
        mutation_id = results[0];
        alpha = stod(results[1]);
        beta = stod(results[2]);
        delta0 = stod(results[3]);
        if (fabs(delta0 - 1.0) < 1e-9) {
            delta0 = 0.999;
        } else if (fabs(delta0) < 1e-9) {
            delta0 = 0.001;
        }
        // Find locus
        if (mut_id2bulk_idx_.count(mutation_id) == 0) {
            cerr << "Hyper parameter contains loci that is not found in the bulk." << endl;
            exit(-1);
        }
        size_t bulk_idx = mut_id2bulk_idx_[mutation_id];
        bulk_data_[bulk_idx]->SetLocuHyperParameters(alpha, beta, delta0);

        line_idx++;
    }
    dat_file.close();
    
}


Interface::Interface(string config_file) {
    // Parse the configuration file
    cout << "Reading config file:" << config_file << endl;
    ProcessConfigFile(config_file, config_);
    
    ReadBulkData();
    if (config_.cn_prior_path != "") {
        ReadCnPrior();
        cn_input_type_ = CopyNumberInputType::TOTAL_CN_PROFILE;
    }
    if (cn_input_type_ == CopyNumberInputType::UNDETERMINED) {
        cerr << "Error: Copy number input type is undetermined.\n";
        exit(-1);
    }
    vector<SingleCellData *> sc_data;
    if (config_.scRNA_file != "") {
        ReadScRnaData();
        ReadScRnaHyperparams();
    }
}

void Interface::Run()
{
    gsl_rng *random = generate_random_object(config_.seed);
    ModelParams params = ModelParams::RandomInit(random,
                                                 config_.alpha0_max,
                                                 config_.lambda_max,
                                                 config_.gamma_max,
                                                 config_.seq_err);
    params.set_sc_dropout_alpha0(config_.sc_dropout_alpha0);
    params.set_sc_dropout_beta0(config_.sc_dropout_beta0);
    params.set_var_cp_prob(config_.var_cp_prob_);

    RunSliceSampler(random, params);
}
