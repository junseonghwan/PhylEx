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
alpha0_(alpha0), lambda_(lambda), gamma_(gamma), seq_error_(seq_err)
{
}

double ModelParams::ComputeAlpha(const vector<size_t> name_vec) const
{
    assert(name_vec.size() > 0);
    size_t j = name_vec.size() - 1;
    double ret = pow(this->lambda_, j) * this->alpha0_;
    return ret;
}

double ModelParams::GetAlpha0Bound(bool max)
{
    return (max ? alpha0_max : alpha0_min);
}

double ModelParams::GetGammaBound(bool max)
{
    return (max ? gamma_max : gamma_min);
}

double ModelParams::GetLambdaBound(bool max)
{
    return (max ? lambda_max : lambda_min);
}

void ModelParams::SetAlpha0(double val)
{
    alpha0_ = val;
}

void ModelParams::SetGamma(double val)
{
    gamma_ = val;
}

void ModelParams::SetLambda(double val)
{
    lambda_ = val;
}

void ModelParams::SetAlpha0Bound(bool max, double val)
{
    if (max) {
        alpha0_max = val;
    } else {
        alpha0_min = val;
    }
}

void ModelParams::SetGammaBound(bool max, double val)
{
    if (max) {
        gamma_max = val;
    } else {
        gamma_min = val;
    }
}

void ModelParams::SetLambdaBound(bool max, double val)
{
    if (max) {
        lambda_max = val;
    } else {
        lambda_min = val;
    }
}

void ModelParams::SetSequencingError(double val)
{
    seq_error_ = val;
}

void ModelParams::SetDirichletConcentrationFactor(double val)
{
    this->dir_conc_multiplicative_factor_ = val;
}

void ModelParams::RandomInit(gsl_rng *random)
{
    alpha0_ = gsl_ran_flat(random, alpha0_min, alpha0_max);
    lambda_ = gsl_ran_flat(random, lambda_min, lambda_max);
    gamma_ = gsl_ran_flat(random, gamma_min, gamma_max);
}
