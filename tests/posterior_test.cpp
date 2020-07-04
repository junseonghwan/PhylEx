//
//  posterior_test.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-07-26.
//

#include <gsl/gsl_rng.h>

//double sample_tssb_param(gsl_rng *random, double min_val, double max_val)
//{
//    return gsl_ran_flat(random, min_val, max_val);
//}
//
//// return the mean of the test function and the variance
//pair<double, double> evaluate_test_function1(const vector<double> &dat)
//{
//    double mean = gsl_stats_mean(dat.data(), 1, dat.size());
//    double var = gsl_stats_variance_m(dat.data(), 1, dat.size(), mean);
//    return make_pair(mean, var);
//}
//
//// return the mean of the test function and the variance
//pair<double, double> evaluate_test_function2(const vector<double> &dat)
//{
//    double tss = gsl_stats_tss_m(dat.data(), 1, dat.size(), 0.0);
//    tss /= dat.size();
//    double var = gsl_stats_tss_m(dat.data(), 1, dat.size(), tss);
//    var /= (dat.size() - 1);
//    return make_pair(tss, var);
//}
//
//bool compare(pair<double, double> mc, size_t M1, pair<double, double> sc, size_t M2, double pval)
//{
//    double sd = sqrt(mc.second/M1 + sc.second/M2);
//    double x = (mc.first - sc.first) / sd;
//    double p = gsl_cdf_gaussian_P(x, 1.0);
//    if (p < pval/2 || p > (1-pval/2)) {
//        return false;
//    }
//    return true;
//}
//
//void geweke_posterior_test()
//{
//    // perform posterior test on the TSSB parameters: alpha0, lambda, gamma
//    // use Uniform for all 3 of these parameters with the following bounds:
//    // alpha0: (0, 50)
//    // lambda: (0, 1)
//    // gamma: (0, 10)
//
//    long seed = 1;
//    size_t M1 = 10000;
//    size_t M2 = 100000;
//    size_t n_mcmc_iter = 1;
//    size_t n_mh_iter = 1;
//    size_t n_data = 10;
//    size_t mean_depth = 1000;
//    size_t max_cn = 3;
//    double alpha0_min = 1, alpha0_max = 50;
//    double lambda_min = 0.2, lambda_max = 1.0;
//    double gamma_min = 0.5, gamma_max = 10;
//    double birth_rate = 0; // no copy number variation
//    double death_rate = 0;
//    gsl_rng *random = generate_random_object(seed);
//
//    // marginal-conditional simulator
//    // 1. sample the parameters and the latent variables from the prior
//    // 2. sample the observation given the parameters and the latent variables
//    // 3. repeat steps 1 and 2 for M1 times
//    // 4. evaluate the test function using the sampled values
//    // Note: if the observed data is not involved in the computation of the test function, it need not be generated
//    vector<double> alpha0s;
//    vector<double> lambdas;
//    vector<double> gammas;
//    for (size_t i = 0; i < M1; i++) {
//        alpha0s.push_back(sample_tssb_param(random, alpha0_min, alpha0_max));
//        lambdas.push_back(sample_tssb_param(random, lambda_min, lambda_max));
//        gammas.push_back(sample_tssb_param(random, gamma_min, gamma_max));
//    }
//
//    cout << "Finished marginal conditional" << endl;
//    // evaluate the test functions
//    pair<double, double> g1_alpha0_mc = evaluate_test_function1(alpha0s);
//    pair<double, double> g1_lambda_mc = evaluate_test_function1(lambdas);
//    pair<double, double> g1_gamma_mc = evaluate_test_function1(gammas);
//
//    pair<double, double> g2_alpha0_mc = evaluate_test_function2(alpha0s);
//    pair<double, double> g2_lambda_mc = evaluate_test_function2(lambdas);
//    pair<double, double> g2_gamma_mc = evaluate_test_function2(gammas);
//
//
//    // successive-conditional simulator
//    // 1. sample the parameters and the latent variables from the prior
//    // 2. sample the observation (data) given the param and latent variables
//    // 3. sample the parameters and the latent variables given the data
//    // 4. repeat steps 2 and 3 for M2 times
//    // 5. use the parameters and the data sampled to evaluate test function
//    double alpha0 = sample_tssb_param(random, alpha0_min, alpha0_max);
//    double lambda = sample_tssb_param(random, lambda_min, lambda_max);
//    double gamma = sample_tssb_param(random, gamma_min, gamma_max);
//    ModelParams params(alpha0, gamma, lambda, 0.01, 1);
//    for (size_t i = 0; i < M2; i++) {
//        if (i % 100 == 0) {
//            cout << "Iter: " << i << endl;
//        }
//        // sample the data
//        vector<BulkDatum *> data;
//        sample_tree_from_prior(random, n_data, params, data);
//        // sample the parameters
//        run_slice_sampler(random, n_mcmc_iter, n_mh_iter, params, &data, 0);
//        // store the parameters
//        alpha0s[i] = params.get_alpha0();
//        lambdas[i] = params.get_lambda();
//        gammas[i] = params.get_gamma();
//    }
//
//    // evaluate the test function
//    pair<double, double> g1_alpha0_sc = evaluate_test_function1(alpha0s);
//    pair<double, double> g1_lambda_sc = evaluate_test_function1(lambdas);
//    pair<double, double> g1_gamma_sc = evaluate_test_function1(gammas);
//
//    pair<double, double> g2_alpha0_sc = evaluate_test_function2(alpha0s);
//    pair<double, double> g2_lambda_sc = evaluate_test_function2(lambdas);
//    pair<double, double> g2_gamma_sc = evaluate_test_function2(gammas);
//
//    double pval = 0.05;
//    bool passed = compare(g1_alpha0_mc, M1, g1_alpha0_sc, M2, pval);
//    cout << "g1(alpha0) test passed: " << passed << endl;
//    passed = compare(g2_alpha0_mc, M1, g2_alpha0_sc, M2, pval);
//    cout << "g2(alpha0) test passed: " << passed << endl;
//
//    passed = compare(g1_lambda_mc, M1, g1_lambda_sc, M2, pval);
//    cout << "g1(lambda) test passed: " << passed << endl;
//    passed = compare(g2_lambda_mc, M1, g2_lambda_sc, M2, pval);
//    cout << "g2(lambda) test passed: " << passed << endl;
//
//    passed = compare(g1_gamma_mc, M1, g1_gamma_sc, M2, pval);
//    cout << "g1(gamma) test passed: " << passed << endl;
//    passed = compare(g2_gamma_mc, M1, g2_gamma_sc, M2, pval);
//    cout << "g2(gamma) test passed: " << passed << endl;
//}
//
//void cook_posterior_test()
//{
//    // for i = 1, ..., I
//    // sample theta^0 from prior
//    // generate one dataset
//    // perform posterior inference for l = 1, ..., L
//    // compare theta^0 against posterior samples to get a test statistic
//    // then, aggregate test stats to decide if the code has a bug or not
//    // this approach is not exact way of testing but can yield some insight on obvious bugs
//}

