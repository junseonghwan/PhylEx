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
    
    double alpha0;
    double lambda;
    double gamma;
    double seq_error = 0.05;
    double amp_error = 0.001;
    double var_cp_prob_ = 0.5;

    double dir_conc_multiplicative_factor = 10.0;

    double alpha0_max = 50.0;
    double alpha0_min = 1.0;
    double alpha0_sigma = 0.1;

    double lambda_max = 0.8;
    double lambda_min = 0.05;
    double lambda_sigma = 0.1;

    double gamma_max = 10.0;
    double gamma_min = 0.1;
    double gamma_sigma = 0.5;

    double sc_dropout_alpha0_ = 0.01;
    double sc_dropout_beta0_ = 1;

    double sc_bursty_alpha0_ = 1;
    double sc_bursty_beta0_ = 0.01;

    // sc_mixture_proportions[0]: Dropout for variant.
    // sc_mixture_proportions[1]: Bursty for variant.
    // sc_mixture_proportions[2]: Bi-allelic distribution.
    double sc_mixture_proportions[3] = {0.05, 0.4, 0.55};

    double b_rate = 0.1;
    double d_rate = 0.1;
    
    size_t max_cn = 5;
    
    gsl_matrix *P_1 = 0, *Q = 0;
    gsl_matrix *U = 0, *U_inv = 0;
    gsl_vector *d = 0;

public:
    ModelParams(double alpha0,
                double gamma,
                double lambda,
                double seq_err);
        
    bool check_bounds();

    inline double get_alpha0() const { return alpha0; }
    inline double get_gamma() const { return gamma; }
    inline double get_lambda() const { return lambda; }
    double get_alpha0_sigma() const;
    double get_lambda_sigma() const;
    double get_gamma_sigma() const;
    inline double get_seq_error() const { return seq_error; }
    inline double get_amp_error() const { return amp_error; }
    inline double get_dir_conc_mult_factor() const { return dir_conc_multiplicative_factor; }
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
    inline double GetScDropoutProportion() const {
        return sc_mixture_proportions[0];
    }
    inline double GetScBurstyVariantProportion() const {
        return sc_mixture_proportions[1];
    }
    inline double GetScBiallelicProportion() const {
        return sc_mixture_proportions[2];
    }
    inline double get_birth_rate() const { return b_rate; }
    inline double get_death_rate() const { return d_rate; }
    inline size_t get_max_cn() const { return max_cn; }
    double alpha(const vector<size_t> name_vec) const;
    
    inline size_t get_max_depth() const { return max_depth; }

    double get_alpha0_bound(bool max);
    double get_gamma_bound(bool max);
    double get_lambda_bound(bool max);
    inline double get_var_cp_prob() const { return var_cp_prob_; }

    void set_alpha0(double val);
    void set_gamma(double val);
    void set_lambda(double val);
    void set_rho(double val);

    void set_alpha0_bound(bool max, double val);
    void set_gamma_bound(bool max, double val);
    void set_lambda_bound(bool max, double val);
    
    void set_max_cn(size_t val);
    void set_var_cp_prob(double var_cp) { var_cp_prob_ = var_cp; }
    
    void set_sc_dropout_alpha0(double sc_dropout_alpha0) {
        sc_dropout_alpha0_ = sc_dropout_alpha0;
    }
    void set_sc_dropout_beta0(double sc_dropout_beta0) {
        sc_dropout_beta0_ = sc_dropout_beta0;
    }
    
    void set_dir_conc_mult_factor(double val);

    void set_birth_rate(double val);
    void set_death_rate(double val);
    void recompute_transition_matrix();
    void compute_transition_matrix(double b, gsl_matrix *P) const;
    double compute_transition_prob(size_t parent, size_t child) const;
    double compute_initial_prob(const gsl_matrix *P, size_t child_r, size_t child_v, size_t root_cn, double p_var = 0.5) const;
    double compute_initial_prob(double branch_length, size_t child_r, size_t child_v, size_t root_cn, double p_var = 0.5) const;
    inline gsl_matrix *get_Q() { return Q; }
    inline gsl_matrix *get_P_1() { return P_1; }
    
    static ModelParams RandomInit(gsl_rng *random,
                                   double alpha0_max,
                                   double lambda_max,
                                   double gamma_max,
                                   double seq_err,
                                  double alpha0_min = 1.0,
                                  double lambda_min = 0.05,
                                  double gamma_min = 0.1);
};

#endif /* tssb_params_h */
