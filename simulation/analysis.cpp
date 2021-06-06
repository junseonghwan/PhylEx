//
// Created by zemp on 18/05/21.
//

#include "analysis.hpp"

/**
 * Compute likelihood for true copy numbers and random perturbations of one clone
 * Write the results in two files, one containing the copy numbers for all the attempts on each line,
 * another with the log-likelihood values for each copy number signal
 *
 * @param n number of perturbations
 * @param sc_data vector of single cell data
 * @param cell2node cell assignments
 * @param simul_config
 * @param gene_set
 * @param output_path
 * @param cn_node node to which the input copy numbers belong
 */
double geneExprLogLikCopyNumber(const vector<SingleCellData *> &sc_data, const vector<CloneTreeNode *> &cell2node,
                                const vector<size_t> &gene_cn, const SimulationConfig &simulationConfig,
                                const vector<Gene *> &gene_set, const CloneTreeNode &cn_node) {
    double logLik = 0;

    for (int c = 0; c < sc_data.size(); ++c) {
        auto cn = cn_node == *cell2node[c] ? gene_cn : cell2node[c]->getGeneCnProfile();
        auto means = compute_means(sc_data[c]->getSizeFactor(), simulationConfig.depth_sf_ratio,
                                    gene_set, cn, simulationConfig.norm_model);

        for (int g = 0; g < gene_cn.size(); ++g) {
            switch (simulationConfig.exprModel) {
                case POISSON:
                    logLik += log_poisson_pdf(sc_data[c]->getExprReads()[g], means[g]);
                    break;
                case NEG_BINOM:
                    logLik += log_negative_binomial_pdf(sc_data[c]->getExprReads()[g], means[g], gene_set[g]->getNbInvDispersion());
                    break;
                case ZIP:
                    logLik += log_zip_pdf(sc_data[c]->getExprReads()[g], means[g], sc_data[c]->getZeroInflationProbs()[g]);
                    break;
                case ZINB:
                    logLik += log_zinb_pdf(sc_data[c]->getExprReads()[g], means[g], gene_set[g]->getNbInvDispersion(),
                                           sc_data[c]->getZeroInflationProbs()[g]);
                    break;
            }
        }
    }

    return logLik;
}

double* compute_means(double size_factor, double depth_sf_ratio, const vector<Gene *> &gene_set,
                      const vector<size_t> &gene_cn, bool norm_model) {
    auto means = new double[gene_set.size()];
    if (norm_model) {
        double norm_factor = 0;
        auto unnormalized_means = new double[gene_set.size()];
        for (int g = 0; g < gene_set.size(); ++g) {
            unnormalized_means[g] = gene_set[g]->getPerCopyExpr() * gene_cn[g];
            norm_factor += unnormalized_means[g];
        }
        for (int g = 0; g < gene_set.size(); ++g) {
            double depth_size = size_factor * depth_sf_ratio;
            means[g] = depth_size * unnormalized_means[g] / norm_factor;
        }
    } else {
        for (int g = 0; g < gene_set.size(); ++g) {
            means[g] = exp(size_factor * gene_set[g]->getPerCopyExpr() * gene_cn[g]);
        }
    }

    return means;
}

double geneExprLogLikCopyNumber(const vector<SingleCellData *> &sc_data, const vector<CloneTreeNode *> &cell2node,
                                const SimulationConfig &simulationConfig, const vector<Gene *> &gene_set) {
    // call the function using true copy number values
    return geneExprLogLikCopyNumber(sc_data, cell2node, cell2node[0]->getGeneCnProfile(),
                                    simulationConfig, gene_set, *cell2node[0]);
}

void cnSensitivitySimulation(size_t n, const vector<SingleCellData *> &sc_data,
                             const vector<CloneTreeNode *> &cell2node, const vector<CloneTreeNode *> &nodes,
                             const SimulationConfig &simul_config, const vector<Gene *> &gene_set, const string& output_path,
                             gsl_rng *rng) {
    ofstream fTrueCn, fCnVar;
    fTrueCn.open(output_path + "/true_cn.txt", ios::out);
    fCnVar.open(output_path + "/cn_variations.txt", ios::out);

    // header of copy number files
    string header = "clone";
    for (auto g : gene_set) {
        header += "," + g->getEnsemblId();
    }
    fTrueCn << header << endl;
    fCnVar << header << endl;

    // save the true copy numbers for each node
    for (auto node: nodes) {
        fTrueCn << node->GetName();
        for (auto cn: node->getGeneCnProfile()) {
            fTrueCn << "," << cn;
        }
        fTrueCn << endl;
    }
    fTrueCn.close();

    // compute the likelihood for the true copy numbers
    ofstream fTrueLik; // copy numbers file and likelihood file
    fTrueLik.open(output_path + "/true_loglik.txt", ios::out);
    fTrueLik << geneExprLogLikCopyNumber(sc_data, cell2node, simul_config, gene_set) << endl;
    fTrueLik.close();

    // compute the likelihood for the perturbation of some copy numbers
    ofstream fLikVar;
    fLikVar.open(output_path + "/var_loglik.txt", ios::out);

    // change 10% of the gene copy numbers on average
    size_t mean_variation = gene_set.size();
    for (int i = 0; i < n; ++i) {
        // sample the clone in which to change the copy number
        size_t cloneIdx = discrete_uniform(rng, nodes.size());
        // copy the true copy number
        vector<size_t> cn = nodes[cloneIdx]->getGeneCnProfile();
        size_t varying_cn = gsl_ran_poisson(rng, mean_variation * (0.1 + i * 0.005));
        // no more than gene_set.size can change
        varying_cn = varying_cn > gene_set.size() ? gene_set.size() : varying_cn;

        // sample the genes in which to change the copy number
        auto samples = sampleIndices(rng, gene_set.size(), varying_cn);

        for (int j = 0; j < varying_cn; ++j) {
            size_t idx = samples[j];
            cn[idx] = perturbateCn(rng, cn[idx], simul_config.max_cn);
        }

        fCnVar << nodes[cloneIdx]->GetName();
        for (auto c: cn) {
            fCnVar << "," << c;
        }
        fCnVar << endl;

        auto logLik = geneExprLogLikCopyNumber(sc_data, cell2node, cn, simul_config, gene_set, *nodes[cloneIdx]);
        fLikVar << logLik << endl;
    }
    fCnVar.close();
    fLikVar.close();
}

size_t perturbateCn(const gsl_rng *rng, int old_cn, int max_cn) {
    int delta = (int) gsl_ran_binomial(rng, 0.5, 9) - 4;
    // this way, p(delta = -1) = p(delta = 1)
    // make delta != 0
    delta = delta <= 0 ? delta - 1 : delta;
    int new_cn = old_cn + delta;
    // make sure the new copy number stays in the boundaries
    return max(min(new_cn, max_cn), 1);
}

// for each clone and for each bin and for each copy number compute the likelihood
// n_clones * n_bins * 7

double* precomputeLikelihoods(const vector<SingleCellData *> &sc_data,
                              const vector<CloneTreeNode *> &cell2node, const vector<CloneTreeNode *> &nodes,
                              const SimulationConfig &simul_config, const vector<Gene *> &gene_set)
{
    // find the node indices for easy access
    map<CloneTreeNode, int> node2idx;
    size_t V = nodes.size();
    for (int i = 0; i < V; ++i) {
        node2idx[*nodes[i]] = i;
    }
    // define the size of the llArray
    size_t B = nodes[0]->getCnProfile().size();
    size_t C = sc_data.size();
    size_t cnRange = simul_config.max_cn + 1;

    // allocate memory and initialize to 0
    auto llArray = new double[V * B * cnRange]{0};

    // compute norm factors for each node
    auto norm_factors = new double[V];
    for (int v = 0; v < V; ++v) {
        norm_factors[v] = computeNormFactor(nodes[v]);
    }

    // iterate through cells, bins and copy number values
    for (int c = 0; c < C; ++c) {
        int v = node2idx[*cell2node[c]];

        for (int b = 0; b < B; ++b) {
            for (int cn = 0; cn < cnRange; ++cn) {
                size_t idx = v * B * cnRange + b * cnRange + cn; // 3D array indexing
                for (auto g: nodes[v]->getBins()->at(b).getGeneIdxs()) {

                    double mean;
                    if (simul_config.norm_model) {
                        double unnorm_mean = sc_data[c]->getSizeFactor() * simul_config.depth_sf_ratio
                                * gene_set[g]->getPerCopyExpr() * cn;
                        // use approximate norm factor so that it's easier
                        mean = unnorm_mean / norm_factors[v];
                    } else {
                        mean = exp(sc_data[c]->getSizeFactor() * gene_set[g]->getPerCopyExpr() * cn);
                    }

                    // compute the partial likelihood for each gene
                    switch (simul_config.exprModel) {
                        case POISSON:
                            llArray[idx] += log_poisson_pdf(sc_data[c]->getExprReads()[g], mean);
                            break;
                        case NEG_BINOM:
                            llArray[idx] += log_negative_binomial_pdf(sc_data[c]->getExprReads()[g],
                                                                      mean, gene_set[g]->getNbInvDispersion());
                            break;
                        case ZIP:
                            llArray[idx] += log_zip_pdf(sc_data[c]->getExprReads()[g], mean,
                                                        sc_data[c]->getZeroInflationProbs()[g]);
                            break;
                        case ZINB:
                            llArray[idx] += log_zinb_pdf(sc_data[c]->getExprReads()[g], mean,
                                                         gene_set[g]->getNbInvDispersion(),
                                                         sc_data[c]->getZeroInflationProbs()[g]);
                            break;
                    }
                }
            }
        }
    }
    return llArray;

}

void testModel(const vector<SingleCellData *> &sc_data,
                 const vector<CloneTreeNode *> &cell2node, const vector<CloneTreeNode *> &nodes,
                 const SimulationConfig &simul_config, const vector<Gene *> &gene_set)
{
    size_t V = nodes.size();
    size_t G = gene_set.size();
    size_t cnRange = simul_config.max_cn + 1;

    // map nodes to set of cells and compute the norm factors for each node
    vector<vector<size_t>> node2cells(V);
    auto norm_factors = new double[V];
    for (int v = 0; v < V; ++v) {
        vector<size_t> cells;
        for (int c = 0; c < sc_data.size(); ++c) {
            if (cell2node[c] == nodes[v]) {
                cells.push_back(c);
            }
        }
        node2cells[v] = cells;
        norm_factors[v] = computeNormFactor(nodes[v]);
    }

    for (int model = POISSON; model < 4; ++model) {

        double acc = 0;
        double mse = 0;
        double mad = 0;
        double prevChoices = 0;
        double infll_count = 0;
        int count = 0;

        for (int v = 0; v < V; ++v) {
            for (int g = 0; g < G; ++g) {
                int cn_pred = -1;
                double maxll = DOUBLE_NEG_INF;

                for (int cn = 0; cn < cnRange; ++cn) {

                    double ll = 0;
                    for (auto c: node2cells[v]) {

                        double mean;
                        if (simul_config.norm_model) {
                            double unnorm_mean = sc_data[c]->getSizeFactor() * simul_config.depth_sf_ratio
                                                 * gene_set[g]->getPerCopyExpr() * cn;
                            // use approximate norm factor so that it's easier
                            mean = unnorm_mean / norm_factors[v];
                        } else {
                            mean = exp(sc_data[c]->getSizeFactor() * gene_set[g]->getPerCopyExpr() * cn);
                        }

                        // for each cell in node v and for each gene,
                        // compute the likelihood of each copy number and take the argmax
                        // compute the partial likelihood for each gene
                        switch (model) {
                            case POISSON:
                                ll += log_poisson_pdf(sc_data[c]->getExprReads()[g], mean);
                                break;
                            case NEG_BINOM:
                                ll += log_negative_binomial_pdf(sc_data[c]->getExprReads()[g],
                                                                mean, gene_set[g]->getNbInvDispersion());
                                break;
                            case ZIP:
                                ll += log_zip_pdf(sc_data[c]->getExprReads()[g], mean,
                                                  sc_data[c]->getZeroInflationProbs()[g]);
                                break;
                            case ZINB:
                                ll += log_zinb_pdf(sc_data[c]->getExprReads()[g], mean,
                                                   gene_set[g]->getNbInvDispersion(),
                                                   sc_data[c]->getZeroInflationProbs()[g]);
                                break;
                            default:
                                cerr << "Model " << model << " doesn't exist!" << endl;
                        }
                    }
                    infll_count += DOUBLE_NEG_INF == ll;
                    if (ll > maxll) {
                        maxll = ll;
                        cn_pred = cn;
                    } else if (ll == maxll && g > 0 && nodes[v]->getGeneCnProfile()[g - 1] == cn) {
                        cn_pred = cn;
                        prevChoices++;
                    }
                }
                int cn_true = nodes[v]->getGeneCnProfile()[g];
                acc += cn_pred == cn_true;
                mse += pow(cn_pred - cn_true, 2);
                mad += abs(cn_pred - cn_true);
                count++;
            }
        }

        cout << (expr_model) model << " | acc: " << acc/count * 100;
        cout << ", mse: " << mse / count;
        cout << ", mad: " << mad / count << endl;
//        cout << "previous choices (%): " << prevChoices / count * 100 << endl;
//        cout << "negative infinite ll values (%): " << infll_count / (count * cnRange) * 100 << endl;
    }
}