//
//  tssb_params.hpp
//  sc-bulk-tumor-tree
//
//  Created by Seong-Hwan Jun on 2018-12-08.
//

#ifndef tssb_params_h
#define tssb_params_h

#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

using namespace std;

double conditional_prob(size_t y, size_t x, double p_var);

class ModelParams
{
    size_t max_depth = 15;
    
    double alpha0_;
    double lambda_;
    double gamma_;
    double seq_error_ = 0.01;
    double var_cp_prob_ = 0.5;

    double dir_conc_multiplicative_factor_ = 10.0;

    double alpha0_max = 50.0;
    double alpha0_min = 1.0;

    double lambda_max = 0.8;
    double lambda_min = 0.05;

    double gamma_max = 10.0;
    double gamma_min = 0.1;

    double sc_dropout_alpha0_ = -1;
    double sc_dropout_beta0_ = -1;

    double sc_bursty_alpha0_ = -1;
    double sc_bursty_beta0_ = -1;

    // sc_mixture_proportions[0]: Dropout for variant.
    // sc_mixture_proportions[1]: Bursty for variant.
    // sc_mixture_proportions[2]: Bi-allelic distribution.
    double sc_mixture_proportions_[3] = {0.05, 0.45, 0.5};

public:
    ModelParams() {}
    ModelParams(double alpha0,
                double gamma,
                double lambda,
                double seq_err);
        
    double GetAlpha0() const { return alpha0_; }
    double GetGamma() const { return gamma_; }
    double GetLambda() const { return lambda_; }
    double GetSequencingError() const { return seq_error_; }
    double GetDirichletConcentrationFactor() const {
        return dir_conc_multiplicative_factor_;
    }
    double GetScBurstyDistributionAlpha0() const {
        return sc_bursty_alpha0_;
    }
    double GetScBurstyDistributionBeta0() const {
        return sc_bursty_beta0_;
    }
    double GetScDropoutDistributionAlpha0() const {
        return sc_dropout_alpha0_;
    }
    double GetScDropoutDistributionBeta0() const {
        return sc_dropout_beta0_;
    }
    double GetScDropoutProportion() const {
        return sc_mixture_proportions_[0];
    }
    double GetScBurstyVariantProportion() const {
        return sc_mixture_proportions_[1];
    }
    double GetScBiallelicProportion() const {
        return sc_mixture_proportions_[2];
    }
    double ComputeAlpha(const vector<size_t> name_vec) const;

    size_t GetMaxDepth() const { return max_depth; }

    double GetAlpha0Bound(bool max);
    double GetGammaBound(bool max);
    double GetLambdaBound(bool max);
    double GetVariantCopyProbability() const { return var_cp_prob_; }

    void SetAlpha0(double val);
    void SetGamma(double val);
    void SetLambda(double val);

    void SetAlpha0Bound(bool max, double val);
    void SetGammaBound(bool max, double val);
    void SetLambdaBound(bool max, double val);
    
    void SetSequencingError(double val);
    
    void SetVariantCopyProbability(double var_cp) { var_cp_prob_ = var_cp; }
    
    void SetSingleCellDropoutAlphaParameter(double sc_dropout_alpha0) {
        sc_dropout_alpha0_ = sc_dropout_alpha0;
    }
    void SetSingleCellDropoutBetaParameter(double sc_dropout_beta0) {
        sc_dropout_beta0_ = sc_dropout_beta0;
    }
    
    void SetSingleCellBurstyAlphaParameter(double sc_bursty_alpha0) {
        sc_bursty_alpha0_ = sc_bursty_alpha0;
    }
    void SetSingleCellBurstyBetaParameter(double sc_bursty_beta0) {
        sc_bursty_beta0_ = sc_bursty_beta0;
    }
    
    void SetDirichletConcentrationFactor(double val);
    
    bool UseSingleCellDropoutDistribution() {
        return (sc_dropout_alpha0_ > 0 && sc_dropout_beta0_ > 0);
    }
    
    void RandomInit(gsl_rng *random);
};

#endif /* tssb_params_h */
