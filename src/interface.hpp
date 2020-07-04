//
//  interface.hpp
//  run
//
//  Created by Seong-Hwan Jun on 2020-06-02.
//

#ifndef interface_hpp
#define interface_hpp

#include <string>

#include "bulk_datum.hpp"
#include "clone_node.hpp"
#include "data_util.hpp"
#include "tssb_state.hpp"

const double DEFAULT_ALPHA0_MAX = 10.0;
const double DEFAULT_LAMBDA_MAX = 1;
const double DEFAULT_GAMMA_MAX = 5.0;

using namespace std;

class Config {
    
public:
    size_t seed;
    string bulk_file = "";
    string scRNA_file = "";
    string sc_hyperparam_file = "";
    string cn_prior_path = "";
    string output_path;
    size_t n_mcmc_iter = 2500;
    size_t n_mh_iter = 5000;
    size_t thinning = 10;
    size_t burn_in = 500;
    size_t output_interval = 100; // interval at which to perform file output operation
    double alpha0_max = DEFAULT_ALPHA0_MAX;
    double lambda_max = DEFAULT_LAMBDA_MAX;
    double gamma_max = DEFAULT_GAMMA_MAX;
    double seq_err = 0.01;
    double var_cp_prob_ = 0.25;
    double sc_dropout_alpha0 = 0.01;
    double sc_dropout_beta0 = 0.01;
};

TSSBState *RunSliceSampler(
    const gsl_rng *random,
    ModelParams &params,
    Config &config,
    vector<BulkDatum *> *bulk_data,
    vector<SingleCellData *> *sc_data,
    CopyNumberInputType cn_input_type);
Config parse_config_file(string config_file_path);

void Run(string config_file);

#endif /* interface_hpp */
