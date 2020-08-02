//
//  tssb_params.cpp
//  tssb
//
//  Created by Seong-Hwan Jun on 2019-05-10.
//

#include "model_params.hpp"

#include <cmath>
#include <iostream>
#include <vector>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>

#include "numerical_utils.hpp"

ModelParams::ModelParams(double alpha0,
                         double gamma,
                         double lambda,
                         double seq_err) :
alpha0(alpha0), lambda(lambda), gamma(gamma), seq_error(seq_err)
{
//    Q = gsl_matrix_alloc(get_max_cn(), get_max_cn());
//    P_1 = gsl_matrix_alloc(get_max_cn(), get_max_cn());
//    U = gsl_matrix_alloc(get_max_cn(), get_max_cn());
//    U_inv = gsl_matrix_alloc(get_max_cn(), get_max_cn());
//    d = gsl_vector_alloc(get_max_cn());

}

bool ModelParams::check_bounds()
{
    if (alpha0 < alpha0_min || alpha0 > alpha0_max) {
        return false;
    }
    
    if (lambda < lambda_min || lambda > lambda_max) {
        return false;
    }
    
    if (gamma < gamma_min || gamma > gamma_max) {
        return false;
    }
    return true;
}

double ModelParams::alpha(const vector<size_t> name_vec) const
{
    assert(name_vec.size() > 0);
    size_t j = name_vec.size() - 1;
    double ret = pow(this->lambda, j) * this->alpha0;
    return ret;
}

double ModelParams::get_alpha0_bound(bool max)
{
    return (max ? alpha0_max : alpha0_min);
}

double ModelParams::get_gamma_bound(bool max)
{
    return (max ? gamma_max : gamma_min);
}

double ModelParams::get_lambda_bound(bool max)
{
    return (max ? lambda_max : lambda_min);
}

double ModelParams::get_alpha0_sigma() const
{
    return alpha0_sigma;
}

double ModelParams::get_lambda_sigma() const
{
    return lambda_sigma;
}

double ModelParams::get_gamma_sigma() const
{
    return gamma_sigma;
}

void ModelParams::set_alpha0(double val)
{
    alpha0 = val;
}

void ModelParams::set_gamma(double val)
{
    gamma = val;
}

void ModelParams::set_lambda(double val)
{
    lambda = val;
}

void ModelParams::set_alpha0_bound(bool max, double val)
{
    if (max) {
        alpha0_max = val;
    } else {
        alpha0_min = val;
    }
}

void ModelParams::set_gamma_bound(bool max, double val)
{
    if (max) {
        gamma_max = val;
    } else {
        gamma_min = val;
    }
}

void ModelParams::set_lambda_bound(bool max, double val)
{
    if (max) {
        lambda_max = val;
    } else {
        lambda_min = val;
    }
}

void ModelParams::set_birth_rate(double val)
{
    b_rate = val;
}

void ModelParams::set_death_rate(double val)
{
    d_rate = val;
}

void ModelParams::set_max_cn(size_t val)
{
    max_cn = val;
}

void ModelParams::set_dir_conc_mult_factor(double val)
{
    this->dir_conc_multiplicative_factor = val;
}

void ModelParams::recompute_transition_matrix()
{
    construct_rate_matrix(b_rate, d_rate, Q);
    eigen(Q, U, U_inv, d);
    construct_transition_matrix(1, U, U_inv, d, P_1); // Compute P_{ij}(1), where 1 is the branch length
}

void ModelParams::compute_transition_matrix(double b, gsl_matrix *P) const
{
    construct_transition_matrix(b, U, U_inv, d, P); // Compute P_{ij}(1), where 1 is the branch
}

double ModelParams::compute_transition_prob(size_t parent, size_t child) const
{
    // look up transition probability
    double p = gsl_matrix_get(P_1, parent, child);
    return p;
}

double conditional_prob(size_t y, size_t x, double p_var)
{
    if (x == 0) {
        if (y > 0) {
            cerr << "Error: " << y << " > " << x << endl;
            exit(-1);
        }
        return 1;
    }
    
    if (y == 0) {
        return 0;
    }
    double num = gsl_ran_binomial_pdf(y, p_var, x);
    double denom = 1 - gsl_ran_binomial_pdf(0, p_var, x);
    return num/denom;
}

double ModelParams::compute_initial_prob(const gsl_matrix *Pb, size_t child_r, size_t child_v, size_t root_cn, double p_var) const
{
    double sum = 0.0;
    for (size_t x = 0; x < get_max_cn(); x++) {
        for (size_t y = 0; y <= x; y++) {
            double p1 = gsl_matrix_get(Pb, root_cn, x);
            double p2 = gsl_matrix_get(P_1, y, child_v);
            double p3 = gsl_matrix_get(P_1, x-y, child_r);
            sum += p1 * conditional_prob(y, x, p_var) * p2 * p3;
        }
    }
    return sum;
}

double ModelParams::compute_initial_prob(double branch_length, size_t child_r, size_t child_v, size_t root_cn, double p_var) const
{
    gsl_matrix *P = gsl_matrix_alloc(get_max_cn(), get_max_cn());
    construct_transition_matrix(branch_length, U, U_inv, d, P);

    double p = compute_initial_prob(P, child_r, child_v, root_cn);
    gsl_matrix_free(P);
    
    return p;
}

ModelParams ModelParams::RandomInit(gsl_rng *random,
                                     double alpha0_max,
                                     double lambda_max,
                                     double gamma_max,
                                     double seq_err,
                                    double alpha0_min,
                                    double lambda_min,
                                    double gamma_min)
{
    double alpha0 = gsl_ran_flat(random, alpha0_min, alpha0_max);
    double lambda = gsl_ran_flat(random, lambda_min, lambda_max);
    double gamma = gsl_ran_flat(random, gamma_min, gamma_max);
    ModelParams model_params(alpha0, gamma, lambda, seq_err);
    model_params.set_alpha0_bound(true, alpha0_max);
    model_params.set_lambda_bound(true, lambda_max);
    model_params.set_gamma_bound(true, gamma_max);
    return model_params;
}
