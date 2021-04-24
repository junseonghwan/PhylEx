//
//  data_util.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-02-01.
//

#include "data_util.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stack>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <gsl/gsl_randist.h>

#include "tssb_state.hpp"
#include "utils.hpp"

NewickNode::NewickNode(string str) {
  // split the string and set name and cellular prevalence
  vector<string> results;
  boost::split(results, str, boost::is_any_of(":"));
  name = results[0];
  cellular_prevalence = stod(results[1]);
}

bool path_exists(string path)
{
    return boost::filesystem::exists(path);
}


gsl_matrix *read_data(string file_name, string sep, bool header)
{
    unsigned int n_patients = 0;
    unsigned int n_genes = 0;

    string line;
    ifstream dat_file (file_name);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << file_name << endl;
        exit(-1);
    }

    vector<string> results;
    vector<double> dat;

    if (header) {
        // skip the first line
        getline(dat_file, line);
    }
    
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of(sep));
        if (n_genes == 0) {
            n_genes = results.size();
        }
        else if (results.size() != n_genes) {
            cerr << "Error in the input file: " << file_name << endl;
            exit(-1);
        }
        for (unsigned int i = 0; i < n_genes; i++) {
            dat.push_back(stod(results[i]));
        }
        n_patients++;
    }
    dat_file.close();

    gsl_matrix *matrix = gsl_matrix_alloc(n_patients, n_genes);
    unsigned int idx = 0;
    for (unsigned int m = 0; m < n_patients; m++) {
        for (unsigned int n = 0; n < n_genes; n++) {
            gsl_matrix_set(matrix, m, n, dat[idx]);
            idx++;
        }
    }

    return matrix;
}

gsl_matrix *read_csv(string file_name, bool header)
{
    string line;
    ifstream dat_file (file_name);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << file_name << endl;
        exit(-1);
    }
    
    vector<string> results;
    vector<double> dat;
    
    if (header) {
        // skip the first line
        getline(dat_file, line);
    }
    
    unsigned int n_cols = 0;
    unsigned int n_rows = 0;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of(","));
        if (n_cols == 0) {
            n_cols = results.size();
            if (n_cols == 0) {
                cerr << "Error: " << file_name << " contains 0 columns." << endl;
                exit(-1);
            }
        }
        else if (results.size() != n_cols) {
            cerr << "Error: Column numbers do not match in the input file: " << file_name << endl;
            exit(-1);
        }
        for (unsigned int i = 0; i < n_cols; i++) {
            dat.push_back(stod(results[i]));
        }
        n_rows++;
    }
    dat_file.close();
    
    gsl_matrix *matrix = gsl_matrix_alloc(n_rows, n_cols);
    unsigned int idx = 0;
    for (unsigned int m = 0; m < n_rows; m++) {
        for (unsigned int n = 0; n < n_cols; n++) {
            gsl_matrix_set(matrix, m, n, dat[idx]);
            idx++;
        }
    }
    
    return matrix;
}


void write_vector(string path, const vector<unsigned int> &data)
{
    ofstream f;
    f.open(path, ios::out);
    for (unsigned int i = 0; i < data.size(); i++) {
        f << data[i] << endl;
    }
    f.close();
}

void write_vector(string path, const vector<double> &data)
{
    ofstream f;
    f.open(path, ios::out);
    for (unsigned int i = 0; i < data.size(); i++) {
        f << data[i] << endl;
    }
    f.close();
}


void write_matrix_as_csv(string path, const gsl_matrix &data)
{
    ofstream f;
    f.open(path, ios::out);
    if (f.is_open()) {
        for (unsigned int i = 0; i < data.size1; i++) {
            for (unsigned int j = 0; j < data.size2; j++) {
                f << gsl_matrix_get(&data, i, j);
                if (j < data.size2 - 1) {
                    f << ", ";
                }
            }
            f << "\n";
        }
        f.close();
    } else {
        cerr << "Error: cannot open file: " << path << endl;
        exit(-1);
    }
}

template <class T>
string ConvertToCommaSeparatedValues(const vector<T> values) {
    string str = "";
    for (size_t i = 0; i < values.size() - 1; i++) {
        str += to_string(values[i]);
        str += ",";
    }
    str += to_string(values.back());
    return str;
}

//string ConvertToCommaSeparatedValues(EigenVectorRef values) {
//    string str = "";
//    for (size_t i = 0; i < values.size() - 1; i++) {
//        str += to_string(values[i]);
//        str += ",";
//    }
//    str += to_string(values[values.size() - 1]);
//    return str;
//}

void WriteCtsCopyNumberProfile(string output_path, const vector<BulkDatum *> &bulk,
                               const vector<pair<double, double> > &cts_cn_profile)
{
    ofstream f;
    f.open(output_path, ios::out);
    f << "ID\tMajorCN\tMinorCN" << endl;
    for (unsigned int i = 0; i < bulk.size(); i++)
    {
        f << bulk[i]->GetId() << "\t";
        f << cts_cn_profile.at(i).first << "\t";
        f << cts_cn_profile.at(i).second << "\n";
    }
    f.close();
}

void WriteBulkData(string output_path, const vector<BulkDatum *> &bulk, bool output_genotype)
{
    ofstream f;
    f.open(output_path, ios::out);
    if (output_genotype) {
        f << "ID\tCHR\tPOS\tREF\tALT\tb\td\tmajor_cn\tminor_cn" << endl;
    } else {
        f << "ID\tCHR\tPOS\tREF\tALT\tb\td\tcn" << endl;
    }
    for (unsigned int i = 0; i < bulk.size(); i++)
    {
        f << bulk[i]->GetId() << "\t";
        f << bulk[i]->GetLocus().get_chr() << "\t";
        f << bulk[i]->GetLocus().get_pos() << "\t";
        //f << "Gene" << "\t"; // Some generic gene name.
        f << "A" << "\t"; // Write arbitrary nucleotide base.
        f << "C" << "\t"; // Write arbitrary nucleotide base.
        f << ConvertToCommaSeparatedValues(bulk[i]->GetVariantReadCount()) << "\t";
        f << ConvertToCommaSeparatedValues(bulk[i]->GetReadCount()) << "\t";
        if (output_genotype) {
            f << ConvertToCommaSeparatedValues(bulk[i]->GetMajorCN()) << "\t";
            f << ConvertToCommaSeparatedValues(bulk[i]->GetMinorCN()) << "\n";
        } else {
            f << ConvertToCommaSeparatedValues(bulk[i]->GetTotalCN()) << "\n";
        }
    }
    f.close();
}

void WriteBetaBinomHp(string output_path,
                      const vector<BulkDatum*> &bulk_data)
{
    ofstream f;
    f.open(output_path + "/beta_binom_hp.csv", ios::out);
    f << "ID\talpha\tbeta\n";
    for (auto bulk_datum : bulk_data) {
        f << bulk_datum->GetId() << "\t";
        f << bulk_datum->GetLocus().get_alpha() << "\t";
        f << bulk_datum->GetLocus().get_beta() << "\n";
    }
    f.close();
}

void WriteScRnaData(string output_path,
                      const vector<BulkDatum*> &bulk_data,
                      const vector<SingleCellData *> &sc_data)
{
    ofstream f;
    f.open(output_path + "/simul_sc.txt", ios::out);
    f << "ID\tCell\ta\td\tSampleName\n";
    for (unsigned int i = 0; i < sc_data.size(); i++)
    {
        SingleCellData *sc = (SingleCellData *)sc_data[i];
        for (size_t loci_idx = 0; loci_idx < bulk_data.size(); loci_idx++)
        {
            auto bulk_datum = bulk_data[loci_idx];
            size_t b = sc_data[i]->GetVariantReads(loci_idx);
            size_t d = sc_data[i]->GetTotalReads(loci_idx);

            // write to file only if it contains non-zero reads.
            if (d == 0) {
                continue;
            }

            f << bulk_datum->GetLocus().get_mutation_id() << "\t";
            f << sc->GetName() << "\t";
            f << (d - b) << "\t" << d << "\tNA\n";
        }
    }
    f.close();
}

void WriteScRnaExpressionData(const string &output_path, const vector<SingleCellExpression *> &sc_expr_data,
                              const vector<Gene *> &gene_set) {
    ofstream f;
    f.open(output_path + "/simul_fc.txt", ios::out);
    for (int i = 0; i < gene_set.size() + 1; ++i) {
        if (i == 0) {
            // write header
            for (auto c: sc_expr_data) {
                f << c->getCellName() << "\t";
            }
            f << endl;
        } else {
            // write genes expressions
            f << gene_set[i-1]->getEnsemblId() << "\t";
            for (auto c: sc_expr_data) {
                f << c->getExprReads()[i-1] << "\t";
            }
            f << endl;
        }
    }
    f.close();
}

/**
 * Writes copy number for each gene and for each clone in TSV file
 *
 * @param output_path data path
 * @param nodes vector of nodes sorted in the desired way
 */
void WriteClonalCNProfiles(const string &output_path, const vector<CloneTreeNode *> &nodes,
                           const vector<Gene *> &gene_set) {
    ofstream f;
    f.open(output_path + "/simul_clonal_cn.txt", ios::out);
    for (int i = 0; i < gene_set.size() + 1; ++i) {
        if (i == 0) {
            // write header
            for (auto n: nodes) {
                f << n->GetName() << "\t";
            }
            f << endl;
        } else {
            // write copy numbers
            f << gene_set[i-1]->getEnsemblId() << "\t";
            for (auto n: nodes) {
                f << n->getCnProfile()[i-1] << "\t";
            }
            f << endl;
        }
    }
    f.close();
}

void WriteCell2NodeAssignment(const string& output_path,
                              const vector<SingleCellData *> &sc_data,
                              const vector<CloneTreeNode *> &cell2node) {
    ofstream f;
    f.open(output_path + "/cell2node.txt", ios::out);
    for (size_t i = 0; i < cell2node.size(); i++) {
        f << sc_data.at(i)->GetName() << "\t" << cell2node.at(i)->GetName() << "\n";
    }
    f.close();
}

void WriteTreeToFile(string output_path,
                const vector<BulkDatum *> &bulk,
                CloneTreeNode *root_node)
{
    if (!path_exists(output_path)) {
        // create path
        boost::filesystem::create_directories(output_path);
    }

    size_t n_muts = bulk.size();

    // construct datum2node mapping
    vector<CloneTreeNode *> all_nodes;
    CloneTreeNode::BreadthFirstTraversal(root_node, all_nodes, false);
    unordered_map<const BulkDatum *, CloneTreeNode *> datum2node;
    CloneTreeNode::Datum2Node(all_nodes, datum2node);

    // ancestral matrix for mutations for computing accuracy on the ordering
    gsl_matrix *A = CloneTreeNode::GetAncestralMatrix(root_node, bulk, datum2node);
    write_matrix_as_csv(output_path + "/ancestral_matrix.csv", *A);

    // cellular prevalence for each of the mutations: Nx2
    ofstream f;
    f.open(output_path + "/cellular_prev.csv", ios::out);
    CloneTreeNode *node;
    for (size_t i = 0; i < n_muts; i++) {
        node = datum2node[bulk[i]];
        string phi_str = ConvertToCommaSeparatedValues(node->NodeParameter().GetCellularPrevalences());
        f << bulk[i]->GetId() << "\t" << phi_str << endl;
    }
    f.close();

    f.open(output_path + "/clone_freq.csv", ios::out);
    for (size_t i = 0; i < n_muts; i++) {
        node = datum2node[bulk[i]];
        string eta_str = ConvertToCommaSeparatedValues(node->NodeParameter().GetCloneFreqs());
        f << bulk[i]->GetId() << "\t" << eta_str << endl;
    }
    f.close();
    
    // write datum to node string
    f.open(output_path + "/datum2node.tsv", ios::out);
    for (size_t i = 0; i < datum2node.size(); i++) {
        f << bulk[i]->GetId()  << "\t" << datum2node[bulk[i]]->GetName() << "\n";
    }
    f.close();

    // output the tree using Newick format
    string newick = write_newick(root_node);
    newick += ";";
    f.open(output_path + "/tree.newick", ios::out);
    f << newick;
    f.close();
    
    // Cellular prevalence for each of the nodes.
    f.open(output_path + "/node2cellular_prev.csv", ios::out);
    for (auto node : all_nodes) {
        f << node->GetName() << "\t" << ConvertToCommaSeparatedValues(node->NodeParameter().GetCellularPrevalences()) << "\n";
    }
    f.close();
}

double parse_newick(string newick, string data_assigment, vector<BulkDatum *> *data)
{
    stack<string> char_stack;
    vector<NewickNode *> children;
    unordered_map<string, NewickNode*> token2node;
    unordered_map<string, CloneTreeNode*> name2clone_node;

    //auto state = new TSSBState<BulkDatum, SingleCellData, CloneTreeNodeParam>();
    stringstream token;
    bool token_end = false;

    for (size_t i = 0; i < newick.size(); i++) {
      string ch(1, newick.at(i));
      if (ch == "(") {
        // start of token
        char_stack.push(ch);
      } else if (ch == ")" || ch == ",") {
        // end of token
        token_end = true;
      } else {
        token << ch;
      }

      if (token_end) {
        // Push the token to stack
        char_stack.push(token.str());
        
        // Create Newick node
        NewickNode *node = new NewickNode(token.str());
        token2node[token.str()] = node;
        
        // If children is nonempty, add children to this node.
        if (children.size() > 0) {
          for (NewickNode *child : children) {
            node->children.push_back(child);
          }
          // Clear children vector.
          children.clear();
        }

        // If we see ")", we need to fill the children vector.
        if (ch == ")") {
          while (true) {
            string ch = char_stack.top();
            char_stack.pop();
            if (ch == "(")
              break;
            // get Newick node
            NewickNode *node = token2node[ch];
            children.push_back(node);
          }
        }
        
        cout << "Read new token: " << token.str() << endl;
        token.str(string());
        token_end = false;
      }
    }
    cout << "Final token: " << token.str() << endl;
    // Fill the children nodes for root.
    NewickNode *root = new NewickNode(token.str());
    if (children.size() > 0) {
      for (NewickNode *child : children) {
        root->children.push_back(child);
      }
    }

    token2node[token.str()] = root;

    // Initialize TSSBState
    // We need CloneTreeNode *root
    // vector<BulkDatum *> *bulk_data
    CloneTreeNode *clone_root = CloneTreeNode::CreateRootNode(1);
    clone_root->NodeParameter().SetRootParameters();
    name2clone_node[root->name] = clone_root;

    queue<NewickNode *> q1;
    queue<CloneTreeNode *> q2;
    q1.push(root);
    q2.push(clone_root);
    CloneTreeNode *parent;
    vector<string> results;
    while (q1.size() > 0) {
        auto *node2 = q1.front();
        parent = q2.front();
        q1.pop();
        q2.pop();
        for (auto it = node2->children.rbegin(); it != node2->children.rend(); ++it) {
            //      boost::split(results, child->name, boost::is_any_of("_"));
            //      size_t child_idx = stoi(results[results.size()-1]);
            auto *child = *it;
            CloneTreeNode *child_node = (CloneTreeNode*)parent->SpawnChild(0.0);
            child_node->SetCellularPrevalenceAtRegion(0, child->cellular_prevalence);
            q1.push(child);
            q2.push(child_node);
            name2clone_node[child->name] = child_node;
        }
    }
    if (q2.size() > 0) {
        cerr << "Error in constructing the tree structure." << endl;
        exit(-1);
    }
    vector<CloneTreeNode *> ret;
    CloneTreeNode::BreadthFirstTraversal(clone_root, ret);
    // TODO add printing with logger
//    for (size_t i = 0; i < ret.size(); i++) {
//        cout << ret.at(i)->Print() << endl;
//    }
  
    gsl_rng *random = generate_random_object(1);
    ModelParams model_params(1, 0.5, 0.5, 0.001);

    auto *tssb = new TSSBState(random,
                               clone_root,
                               model_params,
                               BulkLogLikWithTotalCopyNumber,
                               ScLikelihood,
                               data,
                               0);

    // convert data to map
    unordered_map<string, BulkDatum *> data_map;
    for (BulkDatum *datum : *data) {
        data_map[datum->GetId()] = datum;
    }

    vector<string> results1, results2;
    boost::split(results1, data_assigment, boost::is_any_of(","));
    for (size_t mut_idx = 0; mut_idx < results1.size(); mut_idx++) {
        boost::split(results2, results1[mut_idx], boost::is_any_of(":"));
        string datum_id = results2[0];
        string node_name = results2[1];
        // find node_name and assign bulk_data[idx] to it.
        if (name2clone_node.count(node_name) > 0) {
            CloneTreeNode *n = name2clone_node.at(node_name);
            //tssb->move_datum(n, data_map.at(datum_id), model_params);
            tssb->move_datum(n, mut_idx, model_params);
            cout << datum_id << " assigned to " << node_name << ": " << n->Print() << endl;
        } else {
            cerr << "Error!!!" << endl;
            exit(-1);
        }
    }
    return tssb->compute_log_likelihood_bulk(model_params);
}

//void test_likelihood(string newick_path,
//                     string bulk_data_path,
//                     string bulk_assignment_path)
//{
//    vector<BulkDatum *> *bulk_data = new vector<BulkDatum *>();
//    vector<TSSBState *> states;
//
//    // Read in bulk data
//    read_bulk_data_phyloWGS(bulk_data_path, *bulk_data);
//
//    // Read in the Newick file.
//    ifstream dat_file (newick_path);
//    if (!dat_file.is_open())
//    {
//        cerr << "Could not open the file: " << newick_path << endl;
//        exit(-1);
//    }
//    ifstream bulk_assignment_file(bulk_assignment_path);
//    if (!bulk_assignment_file.is_open())
//    {
//        cerr << "Could not open the file: " << bulk_assignment_path << endl;
//        exit(-1);
//    }
//
//  string line, line2;
//  while ( getline (dat_file, line) && getline (bulk_assignment_file, line2))
//  {
//      vector<string> results;
//      boost::split(results, line, boost::is_any_of(";"));
//
//      double likelihood1 = parse_newick(results[0], line2, bulk_data);
//      double likelihood2 = stod(results[1]);
//      cout << "Likelihood1: " << likelihood1 << endl;
//      cout << "Likelihood2: " << likelihood2 << endl;
//      assert(fabs(likelihood1 - likelihood2) < 1e-6);
//  }
//
//  dat_file.close();
//  bulk_assignment_file.close();
//}

void read_data_assignment(string file_path)
{
  ifstream dat_file (file_path);
  if (!dat_file.is_open())
  {
      cerr << "Could not open the file: " << file_path << endl;
      exit(-1);
  }

  string line;
  while ( getline (dat_file, line) )
  {
    vector<string> results;
    boost::split(results, line, boost::is_any_of(","));
  }
  dat_file.close();
}

string write_newick(CloneTreeNode *node)
{
    unordered_map<size_t, pair<double, CloneTreeNode *> > &children = node->GetIdx2Child();
    if (children.size() == 0) {
        //return (node->get_name() + ":" + to_string(node->get_node_parameter().get_cellular_prev()));
        return node->GetName() + ":" + to_string(1.0);
    }

    string newick = "(";
    for (size_t i = 0; i < children.size(); i++) {
        CloneTreeNode *child = children[i].second;
        string child_newick = write_newick(child);
        newick += child_newick;
        if (i < children.size() - 1) {
            newick += ",";
        }
    }
    //newick += ")" + (node->get_name() + ":" + to_string(node->get_node_parameter().get_cellular_prev())) + ";";
    newick += ")" + node->GetName() + ":" + to_string(1.0);
    return newick;
}

void fill_node_to_param(CloneTreeNode *node,
                        unordered_map<string, vector<double> > &node2param) {
    vector<CloneTreeNode *> all_nodes;
    CloneTreeNode::BreadthFirstTraversal(node, all_nodes, false);
    for (auto node : all_nodes) {
        node2param[node->GetName()] = node->NodeParameter().GetCellularPrevalences();
    }
}

void WriteLogLikToFile(string output_path, double val)
{
    ofstream f;
    f.open(output_path, ios::out);
    f << val << endl;
    f.close();
}

void write_tree(string output_path,
                const vector<BulkDatum *> &bulk,
                CompactTSSBState &state)
{
    if (!path_exists(output_path)) {
        // create path
        boost::filesystem::create_directories(output_path);
    }

    size_t n_muts = bulk.size();

    // ancestral matrix for mutations for computing accuracy on the ordering
    write_matrix_as_csv(output_path + "/ancestral_matrix.csv", *(state.GetAncestralMatrix()));

    // cellular prevalence for each of the mutations: Nx2
    ofstream f;
    f.open(output_path + "/cellular_prev.tsv", ios::out);
    auto params = state.get_param();
    for (size_t i = 0; i < n_muts; i++) {
        auto p = params[i];
        auto phi = p->GetCellularPrevalences();
        f << bulk[i]->GetId() << "\t" << ConvertToCommaSeparatedValues(phi) << endl;
    }
    f.close();

    // write datum to node string
    f.open(output_path + "/datum2node.tsv", ios::out);
    auto datum2node = state.get_datum2node();
    for (size_t i = 0; i < datum2node.size(); i++) {
        f << bulk[i]->GetId()  << "\t" << datum2node[i] << "\n";
    }
    f.close();
    
    //write_vector(output_path + "/cluster_labels.txt", state.get_cluster_labels());
    f.open(output_path + "/cluster_labels.tsv", ios::out);
    auto cluster_labels = state.get_cluster_labels();
    for (size_t i = 0; i < cluster_labels.size(); i++) {
        f << bulk[i]->GetId() << "\t" << cluster_labels[i] << "\n";
    }
    f.close();

    // output the tree using Newick format
    string newick = state.get_newick();
    newick += ";";
    f.open(output_path + "/tree.newick", ios::out);
    f << newick;
    f.close();

    // Cellular prevalence for each of the nodes.
    f.open(output_path + "/node2cellular_prev.tsv", ios::out);
    auto node2param = state.get_node2param();
    for (auto it = node2param.begin(); it != node2param.end(); ++it) {
        f << it->first << "\t" << ConvertToCommaSeparatedValues(it->second) << "\n";
    }
    f.close();
}

void WriteCopyNumberProfileToFile(string output_path,
                      const vector<BulkDatum *> &bulk_data,
                      CloneTreeNode *root_node,
                      unordered_map<CloneTreeNode *, vector<pair<size_t, size_t> > > &cn_profile)
{
    // cn_profile stores copy number profile for each node as vector, where vector is ordered in the same way as the mutations in the bulk input file
    vector<CloneTreeNode *> nodes;
    CloneTreeNode::BreadthFirstTraversal(root_node, nodes, false);

    size_t n_nodes = nodes.size();
    size_t n_muts = bulk_data.size();
    gsl_matrix *cnv_ref_mat = gsl_matrix_alloc(n_muts, n_nodes);
    gsl_matrix *cnv_var_mat = gsl_matrix_alloc(n_muts, n_nodes);

    for (size_t i = 0; i < nodes.size(); i++)
    {
        vector<pair<size_t, size_t> > &vec = cn_profile[nodes[i]];
        for (size_t j = 0; j < vec.size(); j++)
        {
            gsl_matrix_set(cnv_ref_mat, j, i, vec[j].first);
            gsl_matrix_set(cnv_var_mat, j, i, vec[j].second);
        }
    }

    write_matrix_as_csv(output_path + "/cnv_ref.csv", *cnv_ref_mat);
    write_matrix_as_csv(output_path + "/cnv_var.csv", *cnv_var_mat);
}

string convert_chr_to_string(size_t chr)
{
    if (chr == 22) {
        return "X";
    } else if (chr == 23) {
        return "Y";
    } else {
        return to_string(chr);
    }
}

size_t convert_chr_to_int(string chr)
{
    if (chr == "X") {
        return 22;
    } else if (chr == "Y") {
        return 23;
    } else {
        return stoi(chr);
    }
}

