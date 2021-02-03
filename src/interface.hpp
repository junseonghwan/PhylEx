//
//  interface.hpp
//  run
//
//  Created by Seong-Hwan Jun on 2020-06-02.
//

#ifndef interface_hpp
#define interface_hpp

#include <string>
#include <unordered_map>

#include "bulk_datum.hpp"
#include "clone_node.hpp"
#include "data_util.hpp"
#include "tssb_state.hpp"

const double DEFAULT_ALPHA0_MAX = 10.0;
const double DEFAULT_LAMBDA_MAX = 1;
const double DEFAULT_GAMMA_MAX = 5.0;

const size_t BULK_WITH_GENOTYPE_COLUMN_COUNT = 5;
const size_t BULK_WITH_TOTAL_CN_COLUMN_COUNT = 4;
// Not supporting CN_PRIOR option for now.
const size_t BULK_WITH_TOTAL_CN_PRIOR_COLUMN_COUNT = -1;

enum CopyNumberInputType {
    TOTAL_CN, GENOTYPE, TOTAL_CN_PROFILE, UNDETERMINED
};

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
    bool use_geometric_mean = false;
};

class Interface {
    Config config_;
    ModelParams model_params_;
    vector<BulkDatum *> bulk_data_;
    vector<SingleCellData *> sc_data_;
    unordered_map<string, size_t> mut_id2bulk_idx_;
    CopyNumberInputType cn_input_type_;

    TSSBState *RunSliceSampler(const gsl_rng *random,
                               ModelParams &params);

    void ReadBulkData();
    void ReadCnPrior();
    void ReadScRnaData();
    void ReadScRnaHyperparams();

public:
    Interface(string config_file);
    void Run();
    void Print();
};


#endif /* interface_hpp */
