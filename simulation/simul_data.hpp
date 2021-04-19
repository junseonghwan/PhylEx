//
//  simul_data.hpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-08-12.
//

#ifndef simul_data_hpp
#define simul_data_hpp

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gene.hpp>

#include "bulk_datum.hpp"
#include "clone_node.hpp"
#include "tssb_state.hpp"
#include "model_params.hpp"
#include "single_cell.hpp"
#include "simul_config.hpp"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

const bool verbose = false;

CloneTreeNode *SampleFromTssbPrior(size_t region_count,
                                   const gsl_rng *random,
                                   size_t n_data,
                                   ModelParams &model_params,
                                   vector<BulkDatum *> &data);
void CreateSNVs(gsl_rng *random,
                const SimulationConfig &simul_config,
                vector<BulkDatum *> &data);
void GenerateLociPairs(gsl_rng *random,
                       const SimulationConfig &simul_config,
                       const unordered_set<Locus> &somatic_loci,
                       vector<LociPair *> &loci_pairs);
void GenerateBulkData(gsl_rng *random,
                      const SimulationConfig &simul_config,
                      vector<BulkDatum *> &data,
                      CloneTreeNode *root_node);
void GenerateBulkDataWithBDProcess(gsl_rng *random,
                                   const SimulationConfig &simul_config,
                                   vector<BulkDatum *> &data,
                                   CloneTreeNode *root_node,
                                   vector<pair<double, double> > &cts_cn);
vector<CloneTreeNode *> GenerateScRnaData(gsl_rng *random, CloneTreeNode *root_node, const vector<BulkDatum *> &data,
                                          const ModelParams &model_params, const SimulationConfig &simul_config,
                                          vector<SingleCellData *> &sc_data,
                                          vector<SingleCellExpression *> &sc_expr_data,
                                          vector<Gene *> &gene_set);
void GenerateScRnaReads(const gsl_rng *random,
                        const SimulationConfig &simul_config,
                        CloneTreeNode *node,
                        const vector<BulkDatum*> &somatic_loci,
                        vector<size_t> &bulk_sc_coverage,
                        SingleCellData &sc);


void GenerateGenes(
        const gsl_rng *rng,
        const SimulationConfig &simul_config,
        vector<Gene *> &gene_set
);

void GenerateScRnaExpression(
        const gsl_rng *rng,
        const SimulationConfig &simul_config,
        CloneTreeNode *node,
        SingleCellExpression &sc,
        vector<Gene *> &gene_set
);

Eigen::MatrixXf EvolveCn(gsl_rng *random,
                         vector<CloneTreeNode *> &nodes,
                         unordered_map<CloneTreeNode *, size_t> &node2idx,
                         CloneTreeNode *assigned_node,
                         const SimulationConfig &config,
                         Eigen::MatrixXf P0,
                         Eigen::MatrixXf P1);

void EvolveCloneSpecificCN(
        const gsl_rng *rng,
        const SimulationConfig &simul_config,
        CloneTreeNode *root,
        Eigen::MatrixXf &P1
);

//void GenerateScRnaReads(const gsl_rng *random,
//                        const SimulationConfig &simul_config,
//                        const ModelParams &model_params,
//                        CloneTreeNode *node,
//                        const unordered_set<Locus> &somatic_loci,
//                        unordered_map<Locus, pair<size_t, size_t> > &beta_binomial_hps,
//                        unordered_map<Locus, LocusDatum *> &single_site_reads);

//void construct_datum2node(vector<CloneTreeNode *> &all_nodes,
//                          unordered_map<BulkDatum *, Node<BulkDatum, CloneTreeNodeParam> *> &datum2node);
//void generate_germline_SNP(gsl_rng *random,
//                           const SimulationConfig &simul_config,
//                           Locus &somatic_locus);
//void generate_scDNA_data(gsl_rng *random,
//                         CloneTreeNode *root_node,
//                         const vector<BulkDatum *> &bulk_data,
//                         const unordered_map<CloneTreeNode *, vector<pair<size_t, size_t> > > &cn_profile,
//                         const ModelParams &model_params,
//                         const SimulationConfig &simul_config,
//                         vector<SingleCellData *> &sc_data);
//void generate_scDNA_reads(const gsl_rng *random,
//                          size_t mean_depth,
//                          CloneTreeNode *node,
//                          const vector<BulkDatum *> &bulk_data,
//                          const ModelParams &model_params,
//                          const unordered_map<CloneTreeNode *, vector<pair<size_t, size_t> > > &cn_profile,
//                          unordered_map<Locus, LocusDatum *> &reads);

//void GenerateLociPairs(gsl_rng *random,
//                       const SimulationConfig &simul_config,
//                       const unordered_set<Locus> &somatic_loci,
//                       vector<LociPair *> &loci_pairs)
//{
//    for (auto it = somatic_loci.begin(); it != somatic_loci.end(); ++it)
//    {
//        const Locus &somatic_locus = *it;
//        size_t germline_pos = somatic_locus.get_pos() + gsl_rng_uniform_int(random, simul_config.max_germline_locus_proximity);
//        Locus germline_locus(somatic_locus.get_mutation_id(), somatic_locus.get_chr(), germline_pos, "");
//        VariantsPhased var_phased = (gsl_ran_binomial(random, 0.5, 1) > 0) ? VariantsPhased::YES : VariantsPhased::NO;
//        LociPair *loci_pair = new LociPair(germline_locus, somatic_locus, var_phased);
//        loci_pairs.push_back(loci_pair);
//    }
//}

#endif /* simul_data_hpp */
