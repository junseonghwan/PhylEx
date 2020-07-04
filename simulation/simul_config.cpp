//
//  simul_config.cpp
//  lib_simul
//
//  Created by Seong-Hwan Jun on 2019-10-23.
//

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "simul_config.hpp"

vector<double> parse_cn_probs(string val)
{
    vector<double> cn_probs;
    vector<string> results;
    boost::split(results, val, boost::is_any_of(","));
    double sum = 0.0;
    for (size_t i = 0; i < results.size(); i++) {
        cn_probs.push_back(stod(results[i]));
        sum += cn_probs[i];
    }
    assert(abs(sum - 1.0) < 1e-6);
    return cn_probs;
}

void SimulationConfig::insert_option(string key, string val)
{
    try {
        if (key == "seed") {
            this->seed = stoul(val);
        } else if (key == "num_branches") {
            this->num_branches = stoul(val);
        } else if (key == "max_depth") {
            this->max_depth = stoul(val);
        } else if (key == "min_cell_prev") {
            this->min_cell_prev = stod(val);
        } else if (key == "var_cp_prob") {
            this->var_cp_prob = stod(val);
        } else if (key == "var_allele_copy_prob") {
            this->var_allele_copy_prob = parse_cn_probs(val);
        } else if (key == "ref_allele_copy_prob") {
            this->ref_allele_copy_prob = parse_cn_probs(val);
        } else if (key == "n_sites") {
            this->n_sites = stoul(val);
        } else if (key == "bulk_mean_depth") {
            this->bulk_mean_depth = stoul(val);
        } else if (key == "seq_err") {
            this->seq_err = stod(val);
        } else if (key == "birth_rate") {
            this->birth_rate = stod(val);
        } else if (key == "death_rate") {
            this->death_rate = stod(val);
        } else if (key == "max_cn") {
            this->max_cn = stoul(val);
        } else if (key == "n_cells") {
            this->n_cells = stoul(val);
        } else if (key == "sc_mean_depth") {
            this->sc_mean_depth = stoul(val);
        } else if (key == "dropout_rate") {
            this->dropout_rate = stod(val);
        } else if (key == "randomize_dropout") {
            this->randomize_dropout = stoul(val) == 0? false : true;
        } else if (key == "bursty_prob") {
            this->bursty_prob = stod(val);
        } else if (key == "sc_dropout_alpha0") {
            this->sc_dropout_alpha0 = stod(val);
        } else if (key == "sc_dropout_beta0") {
            this->sc_dropout_beta0 = stod(val);
        } else if (key == "output_path") {
            this->output_path = val;
        } else if (key == "n_reps") {
            this->n_reps = stoul(val);
        } else if (key == "n_sims") {
            this->n_sims = stoul(val);
        } else if (key == "beta_binomial_hp_max") {
            this->beta_binomial_hp_max = stoul(val);
        } else if (key == "test_branching") {
            this->test_branching = stoul(val) == 0? false : true;
        } else {
            throw;
        }
    } catch (const std::exception &e) {
        std::cerr << "Cannot convert key: " << key << " and value: " << val << endl;
        exit(-1);
    }
}
