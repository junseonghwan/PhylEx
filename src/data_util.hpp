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

enum CopyNumberInputType {
    TOTAL_CN, GENOTYPE, TOTAL_CN_PROFILE, UNDETERMINED
};

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
void WriteBulkData(string output_path, const vector<BulkDatum *> &bulk, bool output_genotype);
void WriteBulkForPhyloWGS(string output_path, const vector<BulkDatum *> &bulk);
void WriteCNVForPhyloWGS(string output_path); // write an empty cnv file for running PhyloWGS
void write_scDNA_data(string output_path,
                      const vector<SingleCellData *> &sc_data,
                      const vector<BulkDatum *> &bulk);
void WriteBetaBinomHp(string output_path,
                      const vector<BulkDatum*> &loci);
void WriteScRnaData(string output_path,
                    const vector<BulkDatum*> &bulk_data,
                    const vector<SingleCellData *> &sc_data);
void write_tree(string output_path,
                const vector<BulkDatum *> &bulk,
                CompactTSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam> &state);
void WriteTreeToFile(string output_path,
                const vector<BulkDatum *> &bulk,
                Node<BulkDatum,CloneTreeNodeParam> *root_node);

void WriteCopyNumberProfileToFile(string output_path,
                      const vector<BulkDatum *> &bulk_data,
                      Node<BulkDatum,CloneTreeNodeParam> *root_node,
                      unordered_map<Node<BulkDatum,CloneTreeNodeParam> *, vector<pair<size_t, size_t> > > &cn_profile);

CopyNumberInputType ReadBulkData(string bulk_data_path,
                                 vector<BulkDatum *> &bulk_data,
                                 unordered_map<string, Locus *> &somatic_loci);
unordered_set<Locus> read_bulk_data_phyloWGS(string bulk_data_path, vector<BulkDatum *> &bulk_data);
//void read_scDNA_data(string sc_data_path, vector<BulkDatum *> &bulk_data, vector<SingleCellData *> &sc_data);
//void read_scRNA_data(string sc_data_path, unordered_set<Locus> &somatic_loci, vector<SingleCellData *> &sc_data);
void ReadCnPrior(string cn_prior_path, vector<BulkDatum *> &bulk_data);
void ReadScRnaData(string scRNA_data_path,
                     unordered_map<string, Locus *> &id2locus,
                     vector<SingleCellData *> &sc_data);
void ReadScRnaHyperparams(string sc_hyperparam_file, unordered_map<string, Locus *> &id2locus);

// output related functions
void write_vector(string path, const vector<unsigned int> &data);
void write_vector(string path, const vector<double> &data);
void write_matrix_as_csv(string path, const gsl_matrix &data);

string write_newick(Node<BulkDatum,CloneTreeNodeParam> *node);
void fill_node_to_param(Node<BulkDatum,CloneTreeNodeParam> *node,
                        unordered_map<string, double> &node2param);
void test_likelihood(string newick_path,
                     string bulk_data_path,
                     string bulk_assignment_path);
void read_data_assignment(string file_path);

#endif /* data_util_hpp */
