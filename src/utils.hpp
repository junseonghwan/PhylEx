//
//  utils.hpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-07-30.
//

#ifndef utils_hpp
#define utils_hpp

#include <vector>

#include "bulk_datum.hpp"
#include "clone_node.hpp"
#include "single_cell.hpp"
#include "tssb_state.hpp"

double v_measure(double beta,
                 vector<size_t> true_labels,
                 vector<size_t> predicted_labels);

double v_measure(double beta,
                 const vector<BulkDatum *> &data,
                 TSSBState &true_state,
                 TSSBState &sampled_state);

void compute_metrics(const vector<BulkDatum *> &data,
                     TSSBState &true_state,
                     TSSBState &sampled_state);

void call_mutation(const unordered_map<Locus, LocusDatum*> &loci_data, const ModelParams &model_params, unordered_map<Locus, size_t> &mut_map, size_t max_cn = 6);

#endif /* utils_hpp */
