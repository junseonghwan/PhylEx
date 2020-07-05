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

const size_t BULK_WITH_GENOTYPE_COLUMN_COUNT = 9;
// Not supporting it for now.
const size_t BULK_WITH_TOTAL_CN_COLUMN_COUNT = 8;
const size_t BULK_WITH_TOTAL_CN_PRIOR_COLUMN_COUNT = -1;

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
string ConvertToCommaSeparatedValues(vector<T> values) {
    string str = "";
    for (size_t i = 0; i < values.size() - 1; i++) {
        str += values[i];
        str += ",";
    }
    str += values.back();
    return str;
}

string ConvertToCommaSeparatedValues(EigenVectorRef values) {
    string str = "";
    for (size_t i = 0; i < values.size() - 1; i++) {
        str += values[i];
        str += ",";
    }
    str += values[values.size() - 1];
    return str;
}

void WriteBulkData(string output_path, const vector<BulkDatum *> &bulk, bool output_genotype)
{
    ofstream f;
    f.open(output_path, ios::out);
    if (output_genotype) {
        f << "ID\tCHR\tPOS\tGENE\tREF\tALT\tb\td\tmajor_cn\tminor_cn" << endl;
    } else {
        f << "ID\tCHR\tPOS\tGENE\tREF\tALT\tb\td\tcn" << endl;
    }
    for (unsigned int i = 0; i < bulk.size(); i++)
    {
        f << bulk[i]->GetId() << "\t";
        f << bulk[i]->GetLocus().get_chr() << "\t";
        f << bulk[i]->GetLocus().get_pos() << "\t";
        f << "Gene" << "\t"; // Some generic gene name.
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

void write_scDNA_data(string output_path,
                       const vector<SingleCellData *> &sc_data,
                       const vector<BulkDatum *> &bulk)
{
    ofstream f;
    f.open(output_path, ios::out);
    // header: ID\tChromosome\tPosition\tMutantCount\tReferenceCount\tINFO
    for (unsigned int i = 0; i < sc_data.size(); i++)
    {
        SingleCellData *sc = (SingleCellData *)sc_data[i];
        //auto mut_map = sc->get_mutation_map();
        for (size_t j = 0; j < bulk.size(); j++)
        {
            //size_t mut_status = mut_map[bulk[j]->get_locus()];
            size_t mut_status = sc->get_mutation(bulk[j]->GetLocus());
            f << mut_status;
            if (j < bulk.size() - 1) {
                f << " ";
            }
        }
        f << "\n";
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
        f << bulk_datum->GetLocus().get_mutation_id() << "\t";
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
        for (auto bulk_datum : bulk_data)
        {
            const LocusDatum *locus_datum = sc->get_locus_datum(bulk_datum->GetLocus());

            // write to file only if it contains non-zero reads.
            if (locus_datum == 0 || locus_datum->get_n_total_reads() == 0) {
                continue;
            }
            size_t b = locus_datum->get_n_var_reads();
            size_t d = locus_datum->get_n_total_reads();

            f << bulk_datum->GetLocus().get_mutation_id() << "\t";
            f << sc->get_name() << "\t";
            f << (d - b) << "\t" << d << "\tNA\n";
        }
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
    CloneTreeNode::breadth_first_traversal(root_node, all_nodes, false);
    unordered_map<BulkDatum *, CloneTreeNode *> datum2node;
    CloneTreeNode::construct_datum2node(all_nodes, datum2node);

    // ancestral matrix for mutations for computing accuracy on the ordering
    gsl_matrix *A = CloneTreeNode::GetAncestralMatrix(root_node, bulk, datum2node);
    write_matrix_as_csv(output_path + "/ancestral_matrix.csv", *A);

    // cellular prevalence for each of the mutations: Nx2
    ofstream f;
    f.open(output_path + "/cellular_prev.csv", ios::out);
    CloneTreeNode *node;
    for (size_t i = 0; i < n_muts; i++) {
        node = datum2node[bulk[i]];
        string phi_str = ConvertToCommaSeparatedValues(node->get_node_parameter().get_cellular_prevs());
        f << bulk[i]->GetId() << ", " << phi_str << endl;
    }
    f.close();

    f.open(output_path + "/clone_freq.csv", ios::out);
    for (size_t i = 0; i < n_muts; i++) {
        node = datum2node[bulk[i]];
        string eta_str = ConvertToCommaSeparatedValues(node->get_node_parameter().get_clone_freqs());
        f << bulk[i]->GetId() << ", " << eta_str << endl;
    }
    f.close();
    
    // output the tree using Newick format
//    string newick = write_newick(root_node);
//    newick += ";";
//    f.open(output_path + "/tree.newick", ios::out);
//    f << newick;
//    f.close();
    
    // Cellular prevalence for each of the nodes.
    f.open(output_path + "/node2cellular_prev.csv", ios::out);
    for (auto node : all_nodes) {
        f << node->get_name() << ", " << ConvertToCommaSeparatedValues(node->get_node_parameter().get_cellular_prevs()) << "\n";
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
    CloneTreeNode *clone_root = CloneTreeNode::create_root_node(1);
    clone_root->get_node_parameter().SetRootParameters();
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
            CloneTreeNode *child_node = (CloneTreeNode*)parent->spawn_child(0.0);
            child_node->set_cellular_prev(0, child->cellular_prevalence);
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
    CloneTreeNode::breadth_first_traversal(clone_root, ret);
    for (size_t i = 0; i < ret.size(); i++) {
        cout << ret.at(i)->print() << endl;
    }
  
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
            cout << datum_id << " assigned to " << node_name << ": " << n->print() << endl;
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
    unordered_map<size_t, pair<double, CloneTreeNode *> > &children = node->get_idx2child();
    if (children.size() == 0) {
        //return (node->get_name() + ":" + to_string(node->get_node_parameter().get_cellular_prev()));
        return node->get_name() + ":" + to_string(1.0);
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
    newick += ")" + node->get_name() + ":" + to_string(1.0);
    return newick;
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
    //write_matrix_as_csv(output_path + "/ancestral_matrix.csv", *state.get_ancestral_matrix());

    // cellular prevalence for each of the mutations: Nx2
    ofstream f;
    f.open(output_path + "/cellular_prev.csv", ios::out);
    auto params = state.get_param();
    for (size_t i = 0; i < n_muts; i++) {
        auto p = params[i];
        auto phi = p->get_cellular_prevs();
        f << bulk[i]->GetId() << ", " << ConvertToCommaSeparatedValues(phi) << endl;
    }
    f.close();

    // write datum to node string
    f.open(output_path + "/datum2node.txt", ios::out);
    auto datum2node = state.get_datum2node();
    for (size_t i = 0; i < datum2node.size(); i++) {
        f << bulk[i]->GetId()  << "," << datum2node[i] << "\n";
    }
    f.close();
    
    //write_vector(output_path + "/cluster_labels.txt", state.get_cluster_labels());
    f.open(output_path + "/cluster_labels.txt", ios::out);
    auto cluster_labels = state.get_cluster_labels();
    for (size_t i = 0; i < cluster_labels.size(); i++) {
        f << bulk[i]->GetId() << "," << cluster_labels[i] << "\n";
    }
    f.close();

    // output the tree using Newick format
    string newick = state.get_newick();
    newick += ";";
    f.open(output_path + "/tree.newick", ios::out);
    f << newick;
    f.close();
    
    // Cellular prevalence for each of the nodes.
    f.open(output_path + "/node2cellular_prev.csv", ios::out);
    auto node2param = state.get_node2param();
    for (auto it = node2param.begin(); it != node2param.end(); ++it) {
        f << it->first << ", " << it->second << "\n";
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
    CloneTreeNode::breadth_first_traversal(root_node, nodes, false);

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

vector<size_t> ParseRegionalData(string line) {
    vector<string> result;
    boost::split(result, line, boost::is_any_of(","));
    vector<size_t> data(result.size());
    for (size_t i = 0; i < result.size(); i++) {
        data[i] = stoul(result[i]);
    }
    return data;
}

vector<double> ParseRegionalCNData(string line) {
    vector<string> result;
    vector<double> data;
    boost::split(result, line, boost::is_any_of(","));
    for (size_t i = 0; i < result.size(); i++) {
        data[i] = stod(result[i]);
    }
    return data;
}

void ProcessBulkWithTotalCopyNumberProfile(ifstream &dat_file,
                                           vector<BulkDatum *> &bulk_data,
                                           unordered_map<string, Locus *> &somatic_loci) {
    vector<string> results;
    string line;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of("\t"));
        string mut_id = results[0];
        string chr = results[1];
        size_t pos = stol(results[2]);
        auto n_vars = ParseRegionalData(results[6]);
        auto n_reads = ParseRegionalData(results[7]);
        
        if (somatic_loci.count(mut_id) > 0) {
            cerr << "Error: " << mut_id << " already exists!" << endl;
            exit(-1);
        }
        
        Locus *locus = new Locus(mut_id, chr, pos, "");
        BulkDatum *datum = new BulkDatum(mut_id, *locus, n_vars, n_reads);
        bulk_data.push_back(datum);
        somatic_loci[mut_id] = locus;
    }
}

void ProcessBulkWithTotalCopyNumber(ifstream &dat_file,
                                    vector<BulkDatum *> &bulk_data,
                                    unordered_map<string, Locus *> &somatic_loci) {
    vector<string> results;
    string line;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of("\t"));
        string mut_id = results[0];
        string chr = results[1];
        size_t pos = stol(results[2]);
        auto n_vars = ParseRegionalData(results[5]);
        auto n_reads = ParseRegionalData(results[6]);
        auto total_cn = ParseRegionalData(results[7]);

        if (somatic_loci.count(mut_id) > 0) {
            cerr << "Error: " << mut_id << " already exists!" << endl;
            exit(-1);
        }
        
        Locus *locus = new Locus(mut_id, chr, pos, "");
        BulkDatum *datum = new BulkDatum(mut_id, *locus, n_vars, n_reads, total_cn);
        bulk_data.push_back(datum);
        somatic_loci[mut_id] = locus;
    }
}

void ProcessBulkWithGenotype(ifstream &dat_file,
                             vector<BulkDatum *> &bulk_data,
                             unordered_map<string, Locus *> &somatic_loci)
{
    vector<string> results;
    string line;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of("\t"));
        string mut_id = results[0];
        string chr = results[1];
        size_t pos = stol(results[2]);

        vector<size_t> var_reads = ParseRegionalData(results[5]);
        vector<size_t> total_reads = ParseRegionalData(results[6]);
        
        vector<size_t> major_cns = ParseRegionalData(results[7]);
        vector<size_t> minor_cns = ParseRegionalData(results[8]);
        
        if (somatic_loci.count(mut_id) > 0) {
            cerr << "Error: " << mut_id << " already exists!" << endl;
            exit(-1);
        }

        Locus *locus = new Locus(mut_id, chr, pos, "");
        BulkDatum *datum = new BulkDatum(mut_id, *locus,
                                         var_reads, total_reads,
                                         major_cns, minor_cns);
        bulk_data.push_back(datum);
        somatic_loci[mut_id] = locus;
    }
}

CopyNumberInputType ReadBulkData(string bulk_data_path,
                                 vector<BulkDatum *> &bulk_data,
                                 unordered_map<string, Locus *> &somatic_loci)
{
    string line;
    ifstream dat_file (bulk_data_path);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << bulk_data_path << endl;
        exit(-1);
    }

    CopyNumberInputType cn_input_type;
    vector<string> results;

    // Retrieve the first line to determine the input format.
    getline(dat_file, line);
    boost::split(results, line, boost::is_any_of("\t"));

    if (results.size() == BULK_WITH_GENOTYPE_COLUMN_COUNT) {
        ProcessBulkWithGenotype(dat_file, bulk_data, somatic_loci);
        cn_input_type = CopyNumberInputType::GENOTYPE;
    } else if (results.size() == BULK_WITH_TOTAL_CN_COLUMN_COUNT) {
        ProcessBulkWithTotalCopyNumber(dat_file, bulk_data, somatic_loci);
        cn_input_type = CopyNumberInputType::TOTAL_CN;
    } else if (results.size() == BULK_WITH_TOTAL_CN_PRIOR_COLUMN_COUNT) {
        ProcessBulkWithTotalCopyNumberProfile(dat_file, bulk_data, somatic_loci);
        cn_input_type = CopyNumberInputType::UNDETERMINED;
    } else {
        cerr << "Error: invalid bulk input format.\n";
        exit(-1);
    }
    
    dat_file.close();

    return cn_input_type;
}

// cn_prior_path points to a file with the following format:
// First column is the mutation ID used to look up BulkDatum from bulk_data.
// Each of the following column is comma separated, one for each region.
// The second column denotes the probability of copy number being 0.
// The third column denotes the probability of copy number being 1 and so on.
void ReadCnPrior(string cn_prior_path, vector<BulkDatum *> &bulk_data)
{
    ifstream dat_file (cn_prior_path);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << cn_prior_path << endl;
        exit(-1);
    }

    unordered_map<string, BulkDatum *> id2bulk;
    for (auto bulk_datum : bulk_data) {
        id2bulk[bulk_datum->GetId()] = bulk_datum;
    }

    vector<string> results;
    string line;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of(" "));
        if (!id2bulk.count(results[0])) {
            continue;
        }
        auto bulk_datum = id2bulk[results[0]];

        vector<vector<double> > cn_profile;
        for (size_t i = 1; i < results.size(); i++) {
            auto cn_probs = ParseRegionalCNData(results[i]);
            cn_profile.push_back(cn_probs);
        }
        bulk_datum->SetCopyNumberProbs(cn_profile);
    }

    dat_file.close();
}

//void read_scDNA_data(string sc_data_path, vector<BulkDatum *> &bulk_data, vector<SingleCellData *> &sc_data)
//{
//    // matrix of {0, 1, 3}
//    // first column is header
//    string line;
//    ifstream dat_file (sc_data_path);
//    if (!dat_file.is_open())
//    {
//        cerr << "Could not open the file: " << sc_data_path << endl;
//        exit(-1);
//    }
//
//    vector<string> results;
//    vector<double> dat;
//
//    size_t line_idx = 0;
//    while ( getline (dat_file, line) )
//    {
//        boost::split(results, line, boost::is_any_of(" "));
//        if (results.size() != bulk_data.size()) {
//            cerr << "Num mutations for single cell must match the number of mutations in bulk." << endl;
//            exit(-1);
//        }
//        unordered_map<Locus, size_t> mutation_map;
//        for (size_t i = 0; i < results.size(); i++) {
//            size_t mut = stoi(results[i]);
//            const Locus &locus = bulk_data[i]->get_locus();
//            mutation_map[locus] = mut;
//        }
//        SingleCellData *sc = new SingleCellData("c" + to_string(line_idx), mutation_map);
//        line_idx++;
//        sc_data.push_back(sc);
//    }
//    dat_file.close();
//}

//void read_scRNA_data(string scRNA_data_path,
//                     unordered_set<Locus> &somatic_loci,
//                     vector<SingleCellData *> &sc_data)
//{
//    // first column is header
//    // file contains 12 columns:
//    // cell name, somatic chr, somatic pos, germline chr, germline pos,
//    // n_var, n_reads, var_var, var_ref, ref_var, ref_ref, expression level
//    string line;
//    ifstream dat_file (scRNA_data_path);
//    if (!dat_file.is_open())
//    {
//        cerr << "Could not open the file: " << scRNA_data_path << endl;
//        exit(-1);
//    }
//
//    unordered_map<string, SingleCellData *> dat;
//    vector<string> cell_name_orders;
//
//    vector<string> results;
//    size_t line_idx = 0;
//    SingleCellData *sc = 0;
//    while ( getline (dat_file, line) )
//    {
//        if (line_idx == 0) {
//            line_idx++;
//            continue;
//        }
//        boost::split(results, line, boost::is_any_of(","));
//        for (size_t i = 0; i < results.size(); i++) {
//            boost::algorithm::trim(results[i]);
//        }
//        if (dat.count(results[0]) > 0) {
//            sc = dat[results[0]];
//        } else {
//            sc = new SingleCellData(results[0]);
//            sc->set_data_type(SC_DATA_TYPE::RNA);
//            dat[results[0]] = sc;
//            cell_name_orders.push_back(results[0]);
//        }
//
//        Locus somatic_locus(results[1], stoul(results[2]));
//        // check that locus is also found in the bulk data
//        if (somatic_loci.count(somatic_locus) == 0) {
//            cerr << "Single cell RNA data contains loci that is not found in the bulk." << endl;
//            exit(-1);
//        }
//
//        Locus germline_locus(results[3], stoul(results[4]));
//        LociPair *loci_pair = new LociPair(germline_locus, somatic_locus, VariantsPhased::UNKNOWN);
//        size_t n_var_reads = stoul(results[5]);
//        size_t n_reads = stoul(results[6]);
//        LocusDatum *locus_datum = new LocusDatum(n_reads, n_var_reads);
//        size_t var_var = stoul(results[7]);
//        size_t var_ref = stoul(results[8]);
//        size_t ref_var = stoul(results[9]);
//        size_t ref_ref = stoul(results[10]);
//        vector<size_t> paired_reads;
//        paired_reads.push_back(var_var);
//        paired_reads.push_back(var_ref);
//        paired_reads.push_back(ref_var);
//        paired_reads.push_back(ref_ref);
//        LociPairDatum *loci_pair_datum = new LociPairDatum(loci_pair, paired_reads);
//
//        double expr_level = stod(results[11]);
//        sc->insert_read(somatic_locus, locus_datum);
//        sc->set_paired_read(somatic_locus, loci_pair_datum);
//        sc->set_expr_level(somatic_locus, expr_level);
//
//        line_idx++;
//    }
//    dat_file.close();
//
//    for (string cell_name : cell_name_orders) {
//        sc_data.push_back(dat[cell_name]);
//    }
//}

void ReadScRnaData(string scRNA_data_path,
                     unordered_map<string, Locus *> &id2locus,
                     vector<SingleCellData *> &sc_data)
{
    // first column is header
    // file contains 4 columns:
    // mutation id, cell name, a, d
    string line;
    ifstream dat_file (scRNA_data_path);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << scRNA_data_path << endl;
        exit(-1);
    }
    
    vector<string> cell_name_order;
    unordered_map<string, SingleCellData *> dat;
    
    vector<string> results;
    size_t line_idx = 0;
    SingleCellData *sc = 0;
    string mutation_id, cell_name;
    size_t ref_reads, total_reads;
    while ( getline (dat_file, line) )
    {
        if (line_idx == 0) {
            line_idx++;
            continue;
        }
        boost::split(results, line, boost::is_any_of("\t"));
        for (size_t i = 0; i < results.size(); i++) {
            boost::algorithm::trim(results[i]);
        }
        mutation_id = results[0];
        cell_name = results[1];
        ref_reads = stol(results[2]);
        total_reads = stol(results[3]);

        if (dat.count(cell_name) > 0) {
            sc = dat[cell_name];
        } else {
            sc = new SingleCellData(cell_name);
            sc->set_data_type(SC_DATA_TYPE::RNA);
            dat[cell_name] = sc;
            cell_name_order.push_back(cell_name);
        }

        // Find locus
        if (id2locus.count(mutation_id) == 0) {
            cerr << "Single cell RNA data contains loci " << mutation_id << " not found in the bulk." << endl;
            exit(-1);
        }
        const Locus &locus = *id2locus[mutation_id];
        LocusDatum *locus_datum = new LocusDatum(total_reads, total_reads - ref_reads);
        sc->insert_read(locus, locus_datum);
        
        line_idx++;
    }
    dat_file.close();

    for (string cell_name : cell_name_order) {
        sc_data.push_back(dat[cell_name]);
    }
}

void ReadScRnaHyperparams(string sc_hyperparam_file, unordered_map<string, Locus *> &id2locus)
{
    // first column is header
    // file contains 3 columns:
    // mutation id, alpha, beta
    string line;
    ifstream dat_file (sc_hyperparam_file);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << sc_hyperparam_file << endl;
        exit(-1);
    }
    
    vector<string> results;
    size_t line_idx = 0;
    string mutation_id;
    double alpha, beta, delta0;
    while ( getline (dat_file, line) )
    {
        if (line_idx == 0) {
            line_idx++;
            continue;
        }
        boost::split(results, line, boost::is_any_of("\t"));
        for (size_t i = 0; i < results.size(); i++) {
            boost::algorithm::trim(results[i]);
        }
        mutation_id = results[0];
        alpha = stod(results[1]);
        beta = stod(results[2]);
        delta0 = stod(results[3]);
        if (fabs(delta0 - 1.0) < 1e-9) {
            delta0 = 0.999;
        } else if (fabs(delta0) < 1e-9) {
            delta0 = 0.001;
        }
        // Find locus
        if (id2locus.count(mutation_id) == 0) {
            cerr << "Hyper parameter contains loci that is not found in the bulk." << endl;
            exit(-1);
        }
        Locus &locus = *id2locus[mutation_id];
        locus.set_alpha(alpha);
        locus.set_beta(beta);
        locus.set_dropout_prob(delta0);
        
        line_idx++;
    }
    dat_file.close();

}
