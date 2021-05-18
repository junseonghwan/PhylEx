/*
 * sampling_util.hpp
 *
 *  Created on: 30 Mar 2018
 *      Author: seonghwanjun
 */

#ifndef SRC_SAMPLING_UTILS_HPP_
#define SRC_SAMPLING_UTILS_HPP_

#include <vector>

#include <gsl/gsl_rng.h>
#include "numerical_utils.hpp"

using namespace std;

gsl_rng* generate_random_object(long seed);

int sample_birth_death_process(gsl_rng *random, unsigned int num_individuals, double t, double birth_rate, double death_rate);
double weibull(gsl_rng *random, double lambda, double kappa);

double bounded_beta(const gsl_rng *random, double alpha, double beta);
double log_beta_binomial_pdf(size_t k, size_t n, double alpha, double beta);

void uniform(const gsl_rng *random, unsigned int N, double *ret);
double uniform(const gsl_rng *random);
double uniform(const gsl_rng *random, double lower, double upper);

int discrete_uniform(const gsl_rng *random, unsigned long N);
size_t discrete(const gsl_rng *random, vector<double> probs);

// draw one sample and return the index of that sample {1, ..., normalized_probs.size()}
unsigned int multinomial(const gsl_rng *random, vector<double> &normalized_probs);

// ASSERT: the length of result is same as normalized_probs
void multinomial(const gsl_rng *random, unsigned int N, vector<double> &normalized_probs, unsigned int *result);
int multinomial(const gsl_rng *random, const vector<double> &unnormalized_probs, double norm);

// ASSERT: the length of indices is N
// each element of indices taking on values {0, ..., normalized_probs.size()}
void multinomial_sample_indices(const gsl_rng *random, unsigned int N, const vector<double> &normalized_probs, unsigned int *indices);

#endif /* SRC_SAMPLING_UTILS_HPP_ */
