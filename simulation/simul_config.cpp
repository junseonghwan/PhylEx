//
//  simul_config.cpp
//  lib_simul
//
//  Created by Seong-Hwan Jun on 2019-10-23.
//

#include <cmath>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "simul_config.hpp"

/**
 * Parse the copy number probabilities from a string with
 * comma-separated values. Values must sum to 1.
 *
 * @param val configuration value
 * @return probability vector for var/ref copy number of size `max_cn + 1`
 */
vector<double> parse_cn_probs(string val) {
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

SimulationConfig::SimulationConfig() {
    var_allele_copy_prob.push_back(0);
    var_allele_copy_prob.push_back(1);

    ref_allele_copy_prob.push_back(0);
    ref_allele_copy_prob.push_back(1);
}

void SimulationConfig::insert_option(const string &key, const string &val) {
    try {
        if (key == "seed") {
            this->seed = stoul(val);
        } else if (key == "num_branches") {
            this->num_branches = stoul(val);
        } else if (key == "max_depth") {
            this->max_depth = stoul(val);
        } else if (key == "n_sites") {
            this->n_sites = stoul(val);
        } else if (key == "n_regions") {
            this->n_regions = stoul(val);
        } else if (key == "bulk_mean_depth") {
            this->bulk_mean_depth = stoul(val);
        } else if (key == "seq_err") {
            this->seq_err = stod(val);
        } else if (key == "var_allele_copy_prob") {
            this->var_allele_copy_prob = parse_cn_probs(val);
        } else if (key == "ref_allele_copy_prob") {
            this->ref_allele_copy_prob = parse_cn_probs(val);
        } else if (key == "birth_rate") {
            this->birth_rate = stod(val);
        } else if (key == "death_rate") {
            this->death_rate = stod(val);
        } else if (key == "max_cn") {
            this->max_cn = stoul(val);
        } else if (key == "cn_hmm_smoothing_rate") {
            this->cn_hmm_smoothing_rate = stod(val);
        } else if (key == "var_cp_prob") {
            this->var_cp_prob = stod(val);
        } else if (key == "n_cells") {
            this->n_cells = stoul(val);
        } else if (key == "sc_mean_depth") {
            this->sc_mean_depth = stod(val);
        } else if (key == "dropout_rate") {
            this->dropout_rate = stod(val);
        } else if (key == "randomize_dropout") {
            this->randomize_dropout = stoi(val) != 0;
        } else if (key == "bursty_prob") {
            this->bursty_prob = stod(val);
        } else if (key == "sc_bursty_alpha0") {
            this->sc_bursty_alpha0 = stod(val);
        } else if (key == "sc_bursty_beta0") {
            this->sc_bursty_beta0 = stod(val);
        } else if (key == "sc_error_distn_variance_factor") {
            this->sc_error_distn_variance_factor = stod(val);
        } else if (key == "beta_binomial_hp_max") {
            this->beta_binomial_hp_max = stoul(val);
        } else if (key == "randomize_branching") {
            this->randomize_branching = stoi(val) != 0;
        } else if (key == "randomize_cf") {
            this->randomize_cf = stoi(val) != 0;
        } else if (key == "min_cf") {
            this->min_cf = stod(val);
        } else if (key == "snv_sc_sparsity") {
            this->snv_sc_sparsity = stod(val);
        } else if (key == "n_genes") {
            this->n_genes = stoul(val);
        } else if (key == "zero_inflation_alpha") {
            this->zero_inflation_alpha = stod(val);
        } else if (key == "zero_inflation_beta") {
            this->zero_inflation_beta = stod(val);
        } else if (key == "nb_inv_dispersion_shape") {
            this->nb_inv_dispersion_shape = stod(val);
        } else if (key == "nb_inv_dispersion_scale") {
            this->nb_inv_dispersion_scale = stod(val);
        } else if (key == "size_factor_min") {
            this->size_factor_min = stod(val);
        } else if (key == "gene_copy_expr_prob_alpha") {
            this->gene_copy_expr_prob_alpha = stod(val);
        } else if (key == "gene_copy_expr_prob_beta") {
            this->gene_copy_expr_prob_beta = stod(val);
        } else if (key == "size_factor_max") {
            this->size_factor_max = stod(val);
        } else if (key == "genecode_path") {
            this->genecode_path = val;
        } else if (key == "output_path") {
            this->output_path = val;
        } else {
            std::cerr << "Un-recognized key: " << key << endl;
            exit(-1);
        }
    } catch (const std::exception &e) {
        std::cerr << "Cannot convert key: " << key << " for given value: " << val << endl;
        exit(-1);
    }
}

/**
* Parse configuration file. The format is `key: value` for each line
* and must contain all the necessary parameters.
*
* @param config_file_path where the configuration file is located
*/
SimulationConfig *SimulationConfig::parse_config_file(const string &config_file_path) {
    // check required parameters
    set<string> required_params = {"seed", "num_branches", "max_depth", "n_sites", "n_regions",
                                   "bulk_mean_depth", "seq_err", "birth_rate", "death_rate", "max_cn",
                                   "var_cp_prob", "n_cells", "sc_mean_depth", "dropout_rate",
                                   "randomize_dropout", "bursty_prob", "sc_bursty_alpha0", "sc_bursty_beta0",
                                   "sc_error_distn_variance_factor", "beta_binomial_hp_max",
                                   "randomize_branching", "randomize_cf", "min_cf", "snv_sc_sparsity", "n_genes",
                                   "gene_copy_expr_prob_alpha", "gene_copy_expr_prob_beta",
                                   "size_factor_max", "size_factor_min", "cn_hmm_smoothing_rate", "output_path"};

    // create the config object
    auto config = new SimulationConfig();
    string line;
    ifstream config_file(config_file_path);
    if (!config_file.is_open()) {
        throw invalid_argument("could not open the file `" + config_file_path + "`");
    }

    vector<string> results;

    while (getline(config_file, line)) {
        boost::split(results, line, boost::is_any_of(":"));
        boost::algorithm::trim(results[1]);
        config->insert_option(results[0], results[1]);
        required_params.erase(results[0]);
    }

    if (!required_params.empty()) {
        string missing_param_name = *required_params.begin();
        cerr << "Error: parameter `" << missing_param_name << "` is missing in config file" << endl;
        exit(-1);
    }
    config_file.close();

    // determine the type of the model for gene expression data based on the input parameters
    if (config->nb_inv_dispersion_scale == -1) {
        if (config->zero_inflation_beta == -1) {
            config->exprModel = POISSON;
        } else {
            config->exprModel = ZIP;
        }
    } else {
        if (config->zero_inflation_beta == -1) {
            config->exprModel = NEG_BINOM;
        } else {
            config->exprModel = ZINB;
        }
    }

    return config;
}
