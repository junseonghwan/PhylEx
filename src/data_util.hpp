//
//  data_util.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-02-01.
//

#ifndef data_util_hpp
#define data_util_hpp

#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>

#include "bulk_datum.hpp"
#include "compact_tssb_state.hpp"
#include "single_cell.hpp"
#include "tssb_state.hpp"

using namespace std;

class CompactTSSBState;

class NewickNode
{
public:
  string name;
  double cellular_prevalence;
  vector<NewickNode *> children;
  NewickNode(string str);
};

string convert_chr_to_string(size_t chr);
size_t convert_chr_to_int(string chr);

bool path_exists(string path);
void WriteLogLikToFile(string output_path, double val);
void WriteCtsCopyNumberProfile(string output_path, const vector<BulkDatum *> &bulk,
                               const vector<pair<double, double> > &cts_cn_profile);
void WriteBulkData(string output_path, const vector<BulkDatum *> &bulk, bool output_genotype);
void WriteBetaBinomHp(string output_path,
                      const vector<BulkDatum*> &loci);
void WriteScRnaData(string output_path,
                    const vector<BulkDatum*> &bulk_data,
                    const vector<SingleCellData *> &sc_data);
void WriteCell2NodeAssignment(string output_path,
                              const vector<SingleCellData *> &sc_data,
                              const vector<string> &cell2node);
void write_tree(string output_path,
                const vector<BulkDatum *> &bulk,
                CompactTSSBState &state);
void WriteTreeToFile(string output_path,
                const vector<BulkDatum *> &bulk,
                CloneTreeNode *root_node);

void WriteCopyNumberProfileToFile(string output_path,
                      const vector<BulkDatum *> &bulk_data,
                      CloneTreeNode *root_node,
                      unordered_map<CloneTreeNode *, vector<pair<size_t, size_t> > > &cn_profile);

// output related functions
void write_vector(string path, const vector<unsigned int> &data);
void write_vector(string path, const vector<double> &data);
void write_matrix_as_csv(string path, const gsl_matrix &data);

string write_newick(CloneTreeNode *node);
void fill_node_to_param(CloneTreeNode *node,
                        unordered_map<string, vector<double> > &node2param);
void test_likelihood(string newick_path,
                     string bulk_data_path,
                     string bulk_assignment_path);
void read_data_assignment(string file_path);

#endif /* data_util_hpp */
