//
// Created by zemp on 18/05/21.
//

#ifndef PHYLEX_ANALYSIS_H
#define PHYLEX_ANALYSIS_H

#include <vector>
#include "single_cell.hpp"
#include "clone_node.hpp"
#include "simul_config.hpp"

double geneExprLogLikCopyNumber(const vector<SingleCellData *> &sc_data, const vector<CloneTreeNode *> &cell2node,
                                const vector<size_t> &gene_cn, const SimulationConfig &simulationConfig,
                                const vector<Gene *> &gene_set, const CloneTreeNode &cn_node);

// true copy number gene expression likelihood
double geneExprLogLikCopyNumber(const vector<SingleCellData *> &sc_data, const vector<CloneTreeNode *> &cell2node,
                                const SimulationConfig &simulationConfig, const vector<Gene *> &gene_set);

void cnSensitivitySimulation(size_t n, const vector<SingleCellData *> &sc_data,
                             const vector<CloneTreeNode *> &cell2node, const vector<CloneTreeNode *> &nodes,
                             const SimulationConfig &simul_config, const vector<Gene *> &gene_set, const string& output_path,
                             gsl_rng *rng);

size_t perturbateCn(const gsl_rng *rng, int old_cn, int max_cn);


#endif //PHYLEX_ANALYSIS_H