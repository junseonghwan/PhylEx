/*
 * sampling_util.cpp
 *
 *  Created on: 30 Mar 2018
 *      Author: seonghwanjun
 */

#include <algorithm>
#include <limits>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

#include "sampling_utils.hpp"

gsl_rng* generate_random_object(long seed)
{
	// initialize random
    const gsl_rng_type * random_type = gsl_rng_taus;
	gsl_rng *random;

	/* create a generator chosen by the
	     environment variable GSL_RNG_TYPE */
	gsl_rng_env_setup();

	random = gsl_rng_alloc(random_type);
	gsl_rng_set(random, seed);
	return random;
    
}

unsigned int multinomial(const gsl_rng *random, vector<double> &normalized_probs)
{
    double u = gsl_rng_uniform(random);
    double sum = 0.0;
    for (int i = 0; i < normalized_probs.size(); i++)
    {
        sum += normalized_probs[i];
        if (u <= sum) {
            return i;
        }
    }
    
    // probs does not sum to 1: check this and handle the error
    return -1;
}

int sample_birth_death_process(gsl_rng *random, unsigned int num_individuals, double t, double birth_rate, double death_rate)
{
    // first attempt
    unsigned int num_deaths = 0;
    unsigned int num_births = 0;
    for (unsigned int i = 0; i < num_individuals; i++) {
        if (gsl_ran_bernoulli(random, death_rate)) {
            // death event occurred
            num_deaths++;
        } else {
            // sample number of birth events from Poisson distribution
            num_births += gsl_ran_poisson(random, birth_rate);
        }
    }
    return num_births - num_deaths;
}

double weibull(gsl_rng *random, double lambda, double kappa)
{
    return gsl_ran_weibull(random, lambda, kappa);
}

double bounded_beta(const gsl_rng *random, double alpha, double beta)
{
    double x = gsl_ran_beta(random, alpha, beta) - 0.5;
    x = (1 - std::numeric_limits<double>::epsilon()) * x + 0.5;
    return x;
}

double log_beta_binomial_pdf(size_t k, size_t n, double alpha, double beta)
{
    if (k > n)
    {
        return 0;
    }

    double ret = gsl_sf_lnchoose(n, k);
    ret += gsl_sf_lnbeta(k+alpha, n-k+beta);
    ret -= gsl_sf_lnbeta(alpha, beta);
    return ret;
}

// samples an index from uniform distribution over discrete values {0, ..., N-1}
int discrete_uniform(const gsl_rng *random, unsigned long N)
{
    double u = uniform(random);
    double sum = 0.0;
    const double incr = 1./N;
    for (int i = 0; i < N; i++) {
        sum += incr;
        if (u < sum)
            return i;
    }
    return -1; // error
}

void uniform(const gsl_rng *random, unsigned int N, double *ret)
{
    for (int n = 0; n < N; n++)
    {
        ret[n] = gsl_rng_uniform(random);
    }
}

double uniform(const gsl_rng *random)
{
    return uniform(random, 0, 1);
}

double uniform(const gsl_rng *random, double lower, double upper)
{
    double u = gsl_ran_flat(random, lower, upper);
    return u;
}

void multinomial(const gsl_rng *random, unsigned int N, vector<double> &normalized_probs, unsigned int *result)
{
    double probs[normalized_probs.size()];
    std::copy(normalized_probs.begin(), normalized_probs.end(), probs);
    gsl_ran_multinomial(random, N, N, probs, result);
}

int multinomial(const gsl_rng *random, const vector<double> &unnormalized_probs, double norm)
{
    double u = gsl_rng_uniform(random);
    double sum = 0.0;
    for (int i = 0; i < unnormalized_probs.size(); i++)
    {
        sum += (unnormalized_probs[i] / norm);
        if (u <= sum) {
            return i;
        }
    }
    
    // probs does not sum to 1: check this and handle the error
    return -1;
}

void multinomial_sample_indices(const gsl_rng *random, unsigned int N, const vector<double> &normalized_probs, unsigned int *indices)
{
    // draw N uniform values
    double *uvec = new double[N];
    uniform(random, N, uvec); // N
    sort(uvec, uvec + N); // N log N
    double sum = 0.0;
    int idx = 0;
    for (int n = 0; n < N; n++)
    {
        while (uvec[n] > sum + normalized_probs[idx]) {
            sum += normalized_probs[idx];
            idx++;
        }
        if (idx > N) {
            sum = sum + 0;
        }
        indices[n] = idx;
    }
    delete [] uvec;
}

/**
 * Discrete distribution sampler
 *
 * @param random rng object
 * @param probs probability mass function
 * @return outcome of the random event
 */
size_t discrete(const gsl_rng *random, vector<double> probs)
{
    double u = gsl_ran_flat(random, 0, 1);
    double curr = 0.0;
    size_t j = 0;
    for (; j < probs.size(); j++) {
        if (u < curr + probs[j]) {
            break;
        }
        curr += probs[j];
    }
    return j;
}
