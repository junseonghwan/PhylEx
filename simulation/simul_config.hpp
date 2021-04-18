//
//  simul_config.hpp
//  sc-clustering
//
//  Created by Seong-Hwan Jun on 2019-10-18.
//

#ifndef simul_config_h
#define simul_config_h

#include <string>

#include "single_cell.hpp"

using namespace std;

class SimulationConfig
{
public:
    size_t seed;

    // Tree simulation.
    size_t num_branches;
    size_t max_depth;

    // Clonal copy number probabilities.
    // Default values assigned to these vectors are (0, 1) for both var and ref.
    // This means that there will be one variant and one reference copy.
    vector<double> var_allele_copy_prob;
    vector<double> ref_allele_copy_prob;
    double var_cp_prob;
    
    double birth_rate;
    double death_rate;
    size_t max_cn;

    // Bulk data configuration.
    size_t n_sites;
    size_t n_regions;
    size_t bulk_mean_depth;
    double seq_err;

    // Single cell data configuration.
    size_t n_cells;
    double sc_mean_depth;
    double bursty_prob;
    double dropout_rate;
    double sc_bursty_alpha0;
    double sc_bursty_beta0;
    size_t beta_binomial_hp_max;

    // error rate is parametrized with low variance
    double sc_error_distn_variance_factor = 1.0;
    
    string output_path;
    
    bool randomize_dropout = false;
    bool randomize_branching = false;
    bool randomize_cf = true;
    double min_cf = 0.01;
    
    double snv_sc_sparsity = 1; // ratio of bulk DNA SNVs expressed in scRNA data

    // gene expression configuration
    size_t n_genes;
    double zero_inflation_prob; // rho
    double nb_inverse_dispersion; // r

    void insert_option(const string& key, const string& val);
    SimulationConfig();
};

#endif /* simul_config_h */
