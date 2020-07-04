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

#include "bulk_datum.hpp"
#include "tssb_state.hpp"
#include "model_params.hpp"
#include "single_cell.hpp"
#include "simul_config.hpp"

Node<BulkDatum,CloneTreeNodeParam> *SampleFromTssbPrior(const gsl_rng *random,
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
                      Node<BulkDatum,CloneTreeNodeParam> *root_node);
void GenerateBulkDataWithBDProcess(gsl_rng *random,
                                   const SimulationConfig &simul_config,
                                   vector<BulkDatum *> &data,
                                   Node<BulkDatum,CloneTreeNodeParam> *root_node);
void GenerateScRnaData(gsl_rng *random,
                       Node<BulkDatum,CloneTreeNodeParam> *root_node,
                       const vector<BulkDatum *> &data,
                       const ModelParams &model_params,
                       const SimulationConfig &simul_config,
                       vector<SingleCellData *> &sc_data);
void GenerateScRnaReads2(const gsl_rng *random,
                         const SimulationConfig &simul_config,
                         Node<BulkDatum,CloneTreeNodeParam> *node,
                         const vector<BulkDatum*> &somatic_loci,
                         unordered_map<Locus, LocusDatum *> &single_site_reads);
//void GenerateScRnaReads(const gsl_rng *random,
//                        const SimulationConfig &simul_config,
//                        const ModelParams &model_params,
//                        Node<BulkDatum,CloneTreeNodeParam> *node,
//                        const unordered_set<Locus> &somatic_loci,
//                        unordered_map<Locus, pair<size_t, size_t> > &beta_binomial_hps,
//                        unordered_map<Locus, LocusDatum *> &single_site_reads);

//void construct_datum2node(vector<Node<BulkDatum,CloneTreeNodeParam> *> &all_nodes,
//                          unordered_map<BulkDatum *, Node<BulkDatum, CloneTreeNodeParam> *> &datum2node);
//void generate_germline_SNP(gsl_rng *random,
//                           const SimulationConfig &simul_config,
//                           Locus &somatic_locus);
//void generate_scDNA_data(gsl_rng *random,
//                         Node<BulkDatum,CloneTreeNodeParam> *root_node,
//                         const vector<BulkDatum *> &bulk_data,
//                         const unordered_map<Node<BulkDatum,CloneTreeNodeParam> *, vector<pair<size_t, size_t> > > &cn_profile,
//                         const ModelParams &model_params,
//                         const SimulationConfig &simul_config,
//                         vector<SingleCellData *> &sc_data);
//void generate_scDNA_reads(const gsl_rng *random,
//                          size_t mean_depth,
//                          Node<BulkDatum,CloneTreeNodeParam> *node,
//                          const vector<BulkDatum *> &bulk_data,
//                          const ModelParams &model_params,
//                          const unordered_map<Node<BulkDatum,CloneTreeNodeParam> *, vector<pair<size_t, size_t> > > &cn_profile,
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
