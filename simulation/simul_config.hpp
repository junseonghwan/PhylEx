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
    double min_cell_prev;
    
    // Copy number using genotype.
    vector<double> var_allele_copy_prob;
    vector<double> ref_allele_copy_prob;
    double var_cp_prob;
    
    // Copy numbers using birth and death rates.
    double birth_rate;
    double death_rate;
    size_t max_cn;

    // Bulk data configuration.
    size_t n_sites;
    size_t bulk_mean_depth;
    double seq_err;

    // Single cell data configuration.
    size_t n_cells;
    double sc_mean_depth;
    double bursty_prob;
    double dropout_rate;
    double sc_dropout_alpha0;
    double sc_dropout_beta0;
    size_t beta_binomial_hp_max;
    
    size_t n_reps;
    size_t n_sims;
    string output_path;
    
    bool randomize_dropout = false;
    bool test_branching = false;

    void insert_option(string key, string val);
};

#endif /* simul_config_h */
