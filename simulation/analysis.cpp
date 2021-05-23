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
        for (int g = 0; g < gene_cn.size(); ++g) {
            double mean = exp(sc_data[c]->getSizeFactor() * gene_set[g]->getPerCopyExpr() * cn[g]);
            switch (simulationConfig.exprModel) {
                case POISSON:
                    logLik += log_poisson_pdf(sc_data[c]->getExprReads()[g], mean);
                    break;
                case NEG_BINOM:
                    logLik += log_negative_binomial_pdf(sc_data[c]->getExprReads()[g], mean, gene_set[g]->getNbInvDispersion());
                    break;
                case ZIP:
                    logLik += log_zip_pdf(sc_data[c]->getExprReads()[g], mean, sc_data[c]->getZeroInflationProbs()[g]);
                    break;
                case ZINB:
                    logLik += log_zinb_pdf(sc_data[c]->getExprReads()[g], mean, gene_set[g]->getNbInvDispersion(),
                                           sc_data[c]->getZeroInflationProbs()[g]);
                    break;
            }
        }
    }

    return logLik;
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
    if (new_cn < 0) {
        return 0;
    } else if (new_cn > max_cn) {
        return max_cn;
    }
    return new_cn;
}