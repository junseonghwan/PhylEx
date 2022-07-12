
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>

#include <gsl/gsl_rng.h>

#include "clone_node.hpp"
#include "node.hpp"

using namespace std;

class Config {
public:
    size_t seed;
    string tree_file = "";
    string cell_assignment_file = "";
    string gexpr_file = "";
    string gene_info_file = "";
    string output_path;
    size_t max_iter = 100;
};

void ParseConfigFile(string config_file_path, Config &config)
{
    string line;
    ifstream config_file(config_file_path);
    if (!config_file.is_open())
    {
        cerr << "Could not open the file: " << config_file_path << endl;
        exit(-1);
    }
    
    vector<string> results;
    vector<double> dat;
    
    while ( getline (config_file, line) )
    {
        boost::split(results, line, boost::is_any_of(":"));
        boost::algorithm::trim(results[1]);
        if (results[0] == "seed") {
            config.seed = stoul(results[1]);
        } else if (results[0] == "tree_file") {
            config.tree_file = results[1];
        } else if (results[0] == "cell_assignment_file") {
            config.cell_assignment_file = results[1];
        } else if (results[0] == "gexpr_file") {
            config.gexpr_file = results[1];
        } else if (results[0] == "gene_info_file") {
            config.gene_info_file = results[1];
        } else if (results[0] == "output_path") {
            config.output_path = results[1];
        } else if (results[0] == "max_iter") {
            config.max_iter = stoul(results[1]);
        } else {
            cerr << "Unknown option: " << results[0] << endl;
        }
    }
    config_file.close();
}

void ReadTree(const string &tree_file_path)
{
    // Each line contains a node name.
    // We will read in thee nodes, and sort it, then build a tree.
    // Return the root of the tree, which represents the healthy population.
    string line;
    ifstream tree_file(tree_file_path);
    if (!tree_file.is_open())
    {
        cerr << "Could not open the file: " << tree_file_path << endl;
        exit(-1);
    }
    
    vector<string> nodes;
    while ( getline (tree_file, line) )
    {
        boost::algorithm::trim(line);
        nodes.push_back(line);
    }
    // Sort the nodes by length and then lexicographically.
    sort(nodes.begin(), nodes.end(), [](const string &node1, const string &node2) -> bool {
        if (node1.size() < node2.size()) {
            return true;
        } else if (node1.size() > node2.size()) {
            return false;
        } else {
            vector<string> ret1, ret2;
            boost::algorithm::split(ret1, node1, boost::is_any_of("_"));
            boost::algorithm::split(ret2, node2, boost::is_any_of("_"));
            for (size_t i = 0; i < ret1.size(); i++) {
                return (ret1.at(i) < ret2.at(i));
            }
            // if same, return true.
            return true;
        }
    });
    
    Node *healthy_clone = new Node(0, "");
    Node *root = new Node(healthy_clone, nodes.at(0));
    unordered_map<string, Node *> node_str2node_map;
    node_str2node_map[nodes.at(0)] = root;
    for (size_t i = 1; i < nodes.size(); i++) {
        auto &node_str = nodes.at(i);
        auto parent_node_str = CloneTreeNode::RetrieveParentString(node_str);
        Node *node = new Node(node_str2node_map.at(parent_node_str), node_str);
        node_str2node_map[node_str] = node;
    }
    return healthy_clone;
}

void ReadCellAssignment(const string &cell_assignment_file_path)
{
    string line;
    ifstream cell_assignment_file(cell_assignment_file_path);
    if (!cell_assignment_file.is_open())
    {
        cerr << "Could not open the file: " << cell_assignment_file_path << endl;
        exit(-1);
    }

    // Two columns separated by tab: Cell name and Node.
    vector<string> results;
    while ( getline (cell_assignment_file, line) )
    {
        boost::algorithm::split(results, line, "\\s+");
        string cell_name = results[0];
        string node = results[1];
    }
}

void ReadGexpr(const string &gexpr_file)
{
    // N x C matrix with column names given by Cell name that matches what is given in the cell assignment file.
// We assume that the rows are not named but the row index matches the gene info file.
}

void ReadGeneInfo(const string &gene_info_file)
{
    // Four columns: index, chr, start, end.
}

int main(int argc, char *argv[])
{
    string config_file;
    
    namespace po = boost::program_options;
    po::options_description desc("Program options");
    desc.add_options()
    ("help", "Put a help message here.")
    ("config_file,c", po::value<string>(&config_file)->required(), "path to configuration file.")
    ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    
    cout << "Starting program..." << endl;
    Config config;
    ParseConfigFile(config_file, config);
    ReadTree(config.tree_file);
    ReadCellAssignment(config.cell_assignment_file);
    ReadGexpr(config.gexpr_file);
    ReadGeneInfo(config.gene_info_file);
    
    //RunInference();

    return 0;
}
