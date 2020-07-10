/*
 * numerical_utils.hpp
 *
 *  Created on: 6 Mar 2018
 *      Author: seonghwanjun
 */

#ifndef SRC_NUMERICAL_UTILS_HPP_
#define SRC_NUMERICAL_UTILS_HPP_

#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>
#include "Eigen/Dense"

using namespace std;

using EigenVector = Eigen::VectorXd;
using EigenVectorRef = Eigen::Ref<const EigenVector>;
using EigenMatrix = Eigen::MatrixXd;
using EigenMatrixXcd = Eigen::MatrixXcd;
using EigenMatrixRef = Eigen::Ref<const EigenMatrix>;
using EigenMatrixXcdRef = Eigen::Ref<const EigenMatrixXcd>;

using DoubleMatrix = vector<vector<double> >;

extern const double DOUBLE_INF;
extern const double DOUBLE_NEG_INF;
extern const double NaN;

// returns log sum as by-product
double log_dirichlet_pdf(unsigned int K, double *alpha, double *theta);
double log_binomial_pdf(const unsigned int k, const double p, const unsigned int n);
double normalize(const vector<double> &log_weights, vector<double> &weights);
void normalize(const vector<double> &log_weights, vector<double> &weights, double log_norm);
void normalize_destructively(vector<double> &log_weights, double log_norm);
double normalize_destructively(vector<double> &log_weights);
double normalize_destructively(double *log_weights, int size);
double log_add(double x, double y);
double log_subtract(double x, double y);
double log_add(double *x, int size);
double log_add(vector<double> x);
double log_beta_pdf(double x, double a, double b);
double log_prod_beta(vector<double> x, double gamma);

void add(double *x, double c, size_t size);
void multiply(double *x, double c, double *ret, size_t size);

template <typename T>
void print_vector(const vector<T> &v);

void construct_transition_matrix(double t, const gsl_matrix *Q, gsl_matrix *P);
void construct_rate_matrix(double birth_rate, double death_rate, gsl_matrix *Q);
//void svd(const gsl_matrix *Q, gsl_matrix *U, gsl_matrix *V, gsl_matrix *D);
void eigen(const gsl_matrix *Q, gsl_matrix *U, gsl_matrix *U_inv, gsl_vector *d);
void construct_transition_matrix(double t, const gsl_matrix *U, const gsl_matrix *V, const gsl_vector *D, gsl_matrix *P);

#endif /* SRC_NUMERICAL_UTILS_HPP_ */
