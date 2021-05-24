/*
 * numerical_utils.cpp
 *
 *  Created on: 7 Mar 2018
 *      Author: seonghwanjun
 */

#include "sampling_utils.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits>
#include <algorithm>
#include "numerical_utils.hpp"

#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>

const double DOUBLE_INF = std::numeric_limits<double>::infinity();
const double DOUBLE_NEG_INF = -std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();

double log_dirichlet_pdf(unsigned int K, double *alpha, double *theta)
{
    size_t i;
    double log_p = 0.0;
    double sum_alpha = 0.0;

    for (i = 0; i < K; i++)
    {
        log_p += (alpha[i] - 1.0) * log (theta[i]);
    }

    for (i = 0; i < K; i++)
    {
        sum_alpha += alpha[i];
    }

    log_p += gsl_sf_lngamma (sum_alpha);
    
    for (i = 0; i < K; i++)
    {
        log_p -= gsl_sf_lngamma (alpha[i]);
    }
    return log_p;
}

/**
 * Logarithm of Poisson PMF
 *
 * @param k
 * @param mean
 * @return
 */
double log_poisson_pdf(unsigned int k, double mean) {
    if (k == 0 && mean == 0) {
        return 0;
    }
    return k * log(mean) - mean - gsl_sf_lnfact(k);
}

/**
 * Log of Zero-Inflated Poisson
 *
 *  p = rho * 1(k = 0) + (1 - rho) * poisson(k; mean)
 *
 * @param k
 * @param mean
 * @param rho zero-inflation prob
 * @return
 */
double log_zip_pdf(unsigned int k, double mean, double rho) {
    double log_p = log(1 - rho) + log_poisson_pdf(k, mean);
    if (k == 0) {
        log_p = log_add(log_p, log(rho));
    }
    return log_p;
}

/**
 * Computes the negative binomial pdf of a point in the (mean,r) parametrization
 * @param k number of failures
 * @param mean
 * @param r number of successes or, equivalently, inverse of dispersion
 * @return log probability
 */
double negative_binomial_pdf(size_t k, double mean, double r) {
    double p = r/(mean + r);
    return gsl_ran_negative_binomial_pdf(k, p, r);
}

/**
 * Zero-Inflated Negative Binomial pdf
 *
 * @param k outcome of the stochastic event
 * @param mean mean of the negative binomial distribution
 * @param r inverse of dispersion parameter, or, equivalently, number of failures before k successes
 * @param rho zero-inflation probability
 * @return probability
 */
double zinb_pdf(size_t k, double mean, double r, double rho) {
    double p = (1 - rho) * negative_binomial_pdf(k, mean, r);
    if (k == 0) {
        p += rho;
    }
    return p;
}

double log_zinb_pdf(size_t k, double mean, double r, double rho) {
    double log_p = log(1 - rho) + log_negative_binomial_pdf(k, mean, r);
    if (k == 0) {
        log_p = log_add(log_p, log(rho));
    }
    return log_p;
}

/**
 * Parametrized by mean and inverse dispersion parameter
 *
 * @param k
 * @param mean
 * @param r
 * @return
 */
double log_negative_binomial_pdf(unsigned int k, double mean, double r) {
    double log_p;
    if (mean == 0 && k == 0) {
        log_p = 0;
    } else {
        log_p = gsl_sf_lngamma(k + r) -
          gsl_sf_lngamma(k + 1) -
          gsl_sf_lngamma(r) +
          r * (log(r) - log(r + mean)) +
          k * (log(mean) - log(mean + r));
    }
    return log_p;
}

double log_binomial_pdf(const unsigned int k, const double p, const unsigned int n)
{
    if (k > n)
    {
        return 0;
    }
    else
    {
        double logP;
        
        if (p == 0)
        {
            logP = (k == 0) ? 0 : DOUBLE_NEG_INF;
        }
        else if (p == 1)
        {
            logP = (k == n) ? 0 : DOUBLE_NEG_INF;
        }
        else
        {
            double ln_Cnk = gsl_sf_lnchoose (n, k);
            logP = ln_Cnk + k * log (p) + (n - k) * log1p(-p);
        }

        return logP;
    }
}

double normalize(const vector<double> &log_weights, vector<double> &weights)
{
    // compute the lognorm
    double log_norm = log_add(log_weights);
    normalize(log_weights, weights, log_norm);
    return log_norm;
}

void normalize(const vector<double> &log_weights, vector<double> &weights, double log_norm)
{
    double sum = 0.0;
    for (unsigned int i = 0; i < log_weights.size(); i++)
    {
        weights[i] = exp(log_weights[i] - log_norm);
        sum += weights[i];
    }
    if (abs(sum - 1.0) > 1e-4) {
        cerr << "Error in normalization. Check that log_weights and log_norm are correctly calculated." << endl;
        cerr << ceil(sum*100000)/100000.0 << " != 1.0" << endl;
        cerr << "log_norm: " << log_norm << endl;
        exit(-1);
    }
}

void normalize_destructively(vector<double> &log_weights, double log_norm)
{
    // compute the lognorm
    for (unsigned int i = 0; i < log_weights.size(); i++)
    {
        log_weights[i] = exp(log_weights[i] - log_norm);
    }
}

double normalize_destructively(vector<double> &log_weights)
{
    // compute the lognorm
    double log_norm = log_add(log_weights);
    normalize_destructively(log_weights, log_norm);
    return log_norm;
}

double normalize_destructively(double *log_weights, int size)
{
	// compute the lognorm
	double log_norm = log_add(log_weights, size);
	for (int i = 0; i < size; i++)
	{
		log_weights[i] = exp(log_weights[i] - log_norm);
	}
    return log_norm;
}

double log_add(double x, double y)
{
	// make x the max
	if (y > x) {
		double temp = x;
		x = y;
		y = temp;
	}
	// now x is bigger
	if (x == DOUBLE_NEG_INF) {
		return x;
	}
	double negDiff = y - x;
	if (negDiff < -20) { 
		return x;
	}
	return x + log(1.0 + exp(negDiff));
}

// is this useful? or even, make sense?
double log_subtract(double x, double y)
{
    if (x < y) { // (log(e^x - e^y) = log(neg number) = -Inf
        return DOUBLE_NEG_INF;
    }
    double negDiff = y - x;
    if (negDiff < -20) {
        return x;
    }
    return x + log(1.0 - exp(negDiff));
}


double log_add(vector<double> x)
{
    double max = DOUBLE_NEG_INF;
    double maxIndex = 0;
    for (unsigned int i = 0; i < x.size(); i++)
    {
        if (x[i] > max) {
            max = x[i];
            maxIndex = i;
        }
    }
    if (max == DOUBLE_NEG_INF) return DOUBLE_NEG_INF;
    // compute the negative difference
    double threshold = max - 20;
    double sumNegativeDifferences = 0.0;
    for (unsigned int i = 0; i < x.size(); i++) {
        if (i != maxIndex && x[i] > threshold) {
            sumNegativeDifferences += exp(x[i] - max);
        }
    }
    if (sumNegativeDifferences > 0.0) {
        return max + log(1.0 + sumNegativeDifferences);
    } else {
        return max;
    }
    
}

double log_add(double *x, int size)
{
	double max = DOUBLE_NEG_INF;
	double maxIndex = 0;
	for (int i = 0; i < size; i++)
	{
		if (x[i] > max) {
			max = x[i];
			maxIndex = i;
		}
	}
	if (max == DOUBLE_NEG_INF) return DOUBLE_NEG_INF;
	  // compute the negative difference
	  double threshold = max - 20;
	  double sumNegativeDifferences = 0.0;
	  for (int i = 0; i < size; i++) {
	    if (i != maxIndex && x[i] > threshold) {
	      sumNegativeDifferences += exp(x[i] - max);
	    }
	  }
	  if (sumNegativeDifferences > 0.0) {
	    return max + log(1.0 + sumNegativeDifferences);
	  } else {
	    return max;
	  }
}

double log_prod_beta(vector<double> x, double gamma)
{
    double log_sum = 0.0;
    for (double xx : x) {
        double y = log_beta_pdf(xx, 1.0, gamma);
        log_sum += y;
    }
    return log_sum;
}

double log_beta_pdf(double x, double a, double b)
{
    double ret = (a - 1) * log(x);
    ret += (b - 1) * log(1-x);
    ret += gsl_sf_lngamma(a+b);
    ret -= gsl_sf_lngamma(a);
    ret -= gsl_sf_lngamma(b);
    return ret;
}

void add(double *x, double c, size_t size)
{
    for (size_t s = 0; s < size; s++)
    {
        x[s] += c;
    }
}

void multiply(double *x, double c, double *ret, size_t size)
{
    for (size_t s = 0; s < size; s++)
    {
        ret[s] = x[s] * c;
    }
}

template <typename T>
void print_vector(const vector<T> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " ";
    }
}

void construct_rate_matrix(double birth_rate, double death_rate, gsl_matrix *Q)
{
    size_t N = Q->size1;
    assert(Q->size1 == Q->size2);
    assert(N > 0);

    // fill the entries of Q
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            if (i == 0) {
                gsl_matrix_set(Q, i, j, 0.0);
                continue;
            }

            if (i == j) {
                if (i == N-1) {
                    gsl_matrix_set(Q, i, j, -death_rate);
                } else {
                    gsl_matrix_set(Q, i, j, -(birth_rate + death_rate));
                }
            } else if (i+1 == j) {
                gsl_matrix_set(Q, i, j, birth_rate);
            } else if (i-1 == j) {
                gsl_matrix_set(Q, i, j, death_rate);
            } else {
                gsl_matrix_set(Q, i, j, 0.0);
            }
        }
    }
}

// special rate matrix with q_i0 = 0 to prevent the copy number from reaching 0
void construct_rate_matrix1(double birth_rate, double death_rate, gsl_matrix *Q)
{
    size_t N = Q->size1;
    assert(Q->size1 == Q->size2);
    assert(N > 0);

    // fill the entries of Q
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            if (i == 0) {
                gsl_matrix_set(Q, i, j, 0.0);
                continue;
            } else if (i == 1) {
                gsl_matrix_set(Q, i, 0, 0.0);
                continue;
            }
            
            if (i == j) {
                if (i == N-1) {
                    gsl_matrix_set(Q, i, j, -death_rate);
                } else {
                    gsl_matrix_set(Q, i, j, -(birth_rate + death_rate));
                }
            } else if (i+1 == j) {
                gsl_matrix_set(Q, i, j, birth_rate);
            } else if (i-1 == j) {
                gsl_matrix_set(Q, i, j, death_rate);
            } else {
                gsl_matrix_set(Q, i, j, 0.0);
            }
        }
    }
}

//void svd(const gsl_matrix *Q, gsl_matrix *U, gsl_matrix *V, gsl_matrix *D)
//{
//    size_t N = Q->size1;
//    assert(Q->size1 == Q->size2);
//    assert(U->size1 == U->size2);
//    assert(V->size1 == V->size2);
//    assert(D->size1 == D->size2);
//    assert(U->size1 == N);
//    assert(V->size1 == N);
//    assert(D->size1 == N);
//
//    // perform SVD: Q = UDV^T
//    // if Q is diagonalizable, V^T = U^{-1}
//    // e^{Qt} = U e^D V^{T}
//    gsl_vector *d = gsl_vector_alloc(N);
//    gsl_matrix_memcpy(U, Q);
//    gsl_vector *work = gsl_vector_alloc(N);
//    gsl_linalg_SV_decomp(U, V, d, work); // this function over-writes U
//    
//    // construct diagnoal matrix D
//    for (size_t i = 0; i < N; i++) {
//        for (size_t j = 0; j < N; j++) {
//            // initialize entries to 0
//            gsl_matrix_set(D, i, j, 0);
//        }
//        // set the diagnoal elements of D
//        gsl_matrix_set(D, i, i, gsl_vector_get(d, i));
//    }
//
//    gsl_vector_free(d);
//    gsl_vector_free(work);
//}

void eigen(const gsl_matrix *Q, gsl_matrix *U, gsl_matrix *U_inv, gsl_vector *d)
{
    // perform eigen decomposition of Q = U matrix(d) U_inv
    // store the results to U, U_inv, and d
    size_t N = Q->size1;
    gsl_eigen_nonsymmv_workspace *work_space = gsl_eigen_nonsymmv_alloc (N);
    gsl_matrix_complex *temp_U = gsl_matrix_complex_alloc(N, N);
    gsl_vector_complex *temp_d = gsl_vector_complex_alloc(N);
    
    gsl_matrix *Q_copy = gsl_matrix_alloc(N, N);
    gsl_matrix_memcpy(Q_copy, Q);
    gsl_eigen_nonsymmv(Q_copy, temp_d, temp_U, work_space);

    // order the eigen values by their absolute magnitude
    gsl_eigen_nonsymmv_sort(temp_d, temp_U, GSL_EIGEN_SORT_ABS_DESC);

    // check that temp_d and temp_U are purely real-valued
    gsl_complex val;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            val = gsl_matrix_complex_get(temp_U, i, j);
            assert(val.dat[1] == 0); // imaginary part should equal 0
            gsl_matrix_set(U, i, j, val.dat[0]);
        }
        val = gsl_vector_complex_get(temp_d, i);
        assert(val.dat[1] == 0);
        gsl_vector_set(d, i, val.dat[0]);
    }
    
    // compute the inverse of U
    gsl_permutation *perm = gsl_permutation_alloc(N);
    int s;
    gsl_matrix *LU = gsl_matrix_alloc(N, N);
    gsl_matrix_memcpy(LU, U);
    gsl_linalg_LU_decomp(LU, perm, &s);
    gsl_linalg_LU_invert(LU, perm, U_inv);

    gsl_eigen_nonsymmv_free(work_space);
    gsl_permutation_free(perm);
    gsl_matrix_free(LU);
    gsl_matrix_free(Q_copy);
    gsl_matrix_complex_free(temp_U);
    gsl_vector_complex_free(temp_d);
}

void construct_transition_matrix(double t, const gsl_matrix *Q, gsl_matrix *P)
{
    size_t N = Q->size1;
    gsl_matrix *U = gsl_matrix_alloc(N, N);
    gsl_matrix *U_inv = gsl_matrix_alloc(N, N);
    gsl_vector *d = gsl_vector_alloc(N);
    eigen(Q, U, U_inv, d);
    construct_transition_matrix(t, U, U_inv, d, P);
    gsl_matrix_free(U);
    gsl_matrix_free(U_inv);
    gsl_vector_free(d);
}

void construct_transition_matrix(double t, const gsl_matrix *U, const gsl_matrix *U_inv, const gsl_vector *d, gsl_matrix *P)
{
    size_t N = P->size1;
    assert(P->size1 == P->size2);
    assert(U->size1 == U->size2);
    assert(U_inv->size1 == U_inv->size2);
    assert(U->size1 == N);
    assert(U_inv->size1 == N);
    assert(d->size == N);

    gsl_matrix *exp_D = gsl_matrix_alloc(N, N);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            // initialize entries to 0
            gsl_matrix_set(exp_D, i, j, 0);
        }
        // set the diagnoal elements of D
        gsl_matrix_set(exp_D, i, i, exp(t*gsl_vector_get(d, i)));
    }

    // perform matrix multiplications
    gsl_matrix *temp = gsl_matrix_alloc(N, N);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, exp_D, 0, temp);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp, U_inv, 0, P);
    
    for (size_t i = 0; i < N; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < N; j++) {
            sum += gsl_matrix_get(P, i, j);
        }
        assert(abs(sum - 1) < 1e-3);
    }
    
    gsl_matrix_free(temp);
    gsl_matrix_free(exp_D);
}