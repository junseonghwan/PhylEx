
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <fstream>
#include <iostream>
#include <string>

#include "bulk_datum.hpp"
#include "clone_node.hpp"
#include "data_util.hpp"
#include "model_params.hpp"
#include "single_cell.hpp"
#include "simul_config.hpp"
#include "tssb_state.hpp"
#include "utils.hpp"

#include "simul_data.hpp"

using namespace std;

void parse_config_file(string config_file_path, SimulationConfig &config)
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
        config.insert_option(results[0], results[1]);
    }
    config_file.close();
}

void CreateLinearTree(gsl_rng *random, CloneTreeNode *root, size_t max_depth) {
    vector<Node<BulkDatum,CloneTreeNodeParam> *> nodes;
    nodes.push_back(root);
    auto node = root;
    for (size_t i = 0; i < max_depth; i++) {
        node = (CloneTreeNode*)node->spawn_child(1);
        nodes.push_back(node);
    }
    // Set clone frequency.
    size_t n_nodes = max_depth + 1;
    double *cell_prev = new double[n_nodes];
    double *clone_freq = new double[n_nodes];
    cell_prev[0] = 1.0;
    cell_prev[1] = gsl_ran_flat(random, 0.5, 1.0);
    for (size_t i = 2; i < n_nodes; i++) {
        cell_prev[i] = cell_prev[i-1]/2;
    }
    for (size_t i = 0; i < n_nodes; i++) {
        clone_freq[i] = cell_prev[i] - cell_prev[i+1];
        ((CloneTreeNode *)nodes[i])->set_cellular_prev(cell_prev[i]);
        ((CloneTreeNode *)nodes[i])->set_clone_freq(clone_freq[i]);
    }

    for (size_t i = 0; i < n_nodes; i++) {
        cout << nodes[i]->print() << endl;
    }
}

void CreateNaryTree(gsl_rng *random,
                    CloneTreeNode *root,
                    size_t max_depth,
                    size_t num_branches,
                    double min_cell_prev) {
    deque<CloneTreeNode *> queue;
    queue.push_back(root);
    root->set_clone_freq(1.0);
    root->set_cellular_prev(1.0);
    while (!queue.empty()) {
        auto node = queue.front();
        queue.pop_front();
        if (node->get_node_parameter().get_cellular_prev() < min_cell_prev ||
            node->get_depth() >= max_depth) {
            continue;
        }

        for (size_t i = 0; i < num_branches; i++) {
            auto child_node = (CloneTreeNode*)node->spawn_child(0.5);
            // Set cellular prevalence.
            double parent_clone_freq = node->get_node_parameter().get_clone_freq();
            double u = 0.5;
            if (random != 0) {
                u = gsl_ran_flat(random, 0, 1);
            }
            double cell_prev = u * parent_clone_freq;
            child_node->set_clone_freq(cell_prev);
            child_node->set_cellular_prev(cell_prev);
            node->set_clone_freq(parent_clone_freq - cell_prev);
            queue.push_back(child_node);
        }
    }
    vector<Node<BulkDatum, CloneTreeNodeParam> *> nodes;
    Node<BulkDatum, CloneTreeNodeParam>::breadth_first_traversal(root,
                                                                 nodes);
    for (size_t i = 0; i < nodes.size(); i++) {
        cout << nodes[i]->print() << endl;
    }
}

int main(int argc, char *argv[])
{
    string config_file_path;
    
    namespace po = boost::program_options;
    po::options_description desc("Program options");
    desc.add_options()
    ("help", "Put a help message here.")
    ("config,c", po::value<string>(&config_file_path)->required(), "path to config file.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    // Parse the config file.
    SimulationConfig simul_config;
    parse_config_file(config_file_path, simul_config);

    gsl_rng *rand = generate_random_object(simul_config.seed);
    auto model_params = ModelParams::RandomInit(rand, 10, 1, 10, simul_config.seq_err);
    model_params.set_var_cp_prob(simul_config.var_cp_prob);
    
    bool bd_process = false;
    if (simul_config.birth_rate > 0 && simul_config.max_cn > 2) {
        bd_process = true;
        model_params.set_birth_rate(simul_config.birth_rate);
        model_params.set_death_rate(simul_config.death_rate);
    } else if (simul_config.ref_allele_copy_prob.size() >= 1 &&
               simul_config.var_allele_copy_prob.size() >= 2) {
        bd_process = false;
    } else {
        cerr << "Error in copy number simulation configuration." << endl;
        cerr << "For BD simulation, specify birth, death rates, and max_cn." << endl;
        cerr << "For clonal copy number, specify prior over copy number with at least 3 entries." << endl;
    }
    
    for (size_t n = 0; n < simul_config.n_sims; n++) {
        string sim_path = simul_config.output_path + "/sim" + to_string(n);
        gsl_rng *random = generate_random_object(gsl_rng_get(rand));

        auto root_node = CloneTreeNode::create_root_node();
        if (simul_config.num_branches == 1) {
            CreateLinearTree(rand, root_node, simul_config.max_depth);
        } else {
            if (simul_config.test_branching) {
                CreateNaryTree(0,
                               root_node,
                               simul_config.max_depth,
                               simul_config.num_branches,
                               simul_config.min_cell_prev);
            } else {
                CreateNaryTree(rand,
                               root_node,
                               simul_config.max_depth,
                               simul_config.num_branches,
                               simul_config.min_cell_prev);
            }
        }
        
        vector<BulkDatum *> data;
        
        // Create somatic SNVs: chr and pos.
        CreateSNVs(random, simul_config, data);
        // Incorporate copy number variation using BD (non-clonal) or genotype (clonal).
        if (bd_process) {
            GenerateBulkDataWithBDProcess(random,
                                          simul_config,
                                          data,
                                          root_node);
        } else {
            GenerateBulkData(random,
                             simul_config,
                             data,
                             root_node);
        }

        // Generate single cells data.
        // Output all simulation data.
        for (size_t n = 0; n < simul_config.n_reps; n++) {

            string output_path = sim_path + "/rep" + to_string(n) + "/";

            // Create directories for output.
            boost::filesystem::path outpath(output_path);
            boost::filesystem::create_directories(outpath);

            WriteBulkData(output_path + "/simul_ssm.txt", data, false);
            WriteBulkData(output_path + "/genotype_ssm.txt", data, true);

            // Generate single cell reads.
            vector<SingleCellData *> sc_data;
            GenerateScRnaData(random,
                              root_node,
                              data,
                              model_params,
                              simul_config,
                              sc_data);
            WriteScRnaData(output_path, data, sc_data);
            WriteBetaBinomHp(output_path, data);

            // Output information needed for evaluation.
            vector<Node<BulkDatum,CloneTreeNodeParam> *> all_nodes;
            Node<BulkDatum,CloneTreeNodeParam>::breadth_first_traversal(root_node,
                                                                        all_nodes,
                                                                        false);
            unordered_map<BulkDatum *, Node<BulkDatum, CloneTreeNodeParam> *> datum2node;
            Node<BulkDatum,CloneTreeNodeParam>::construct_datum2node(all_nodes, datum2node);
            double bulk_log_lik = 0.0;
            for (size_t i = 0; i < simul_config.n_sites; i++) {
                bulk_log_lik += BulkLogLikWithGenotype(datum2node[data[i]], data[i], model_params);
            }
            double sc_log_lik =
                TSSBState<BulkDatum,SingleCellData,CloneTreeNodeParam>::compute_log_likelihood_sc(
                          root_node, data, sc_data, ScLikelihood, model_params);
            WriteTreeToFile(output_path, data, root_node);
            WriteLogLikToFile(output_path + "/log_lik_bulk.txt", bulk_log_lik);
            WriteLogLikToFile(output_path + "/log_lik_sc.txt", sc_log_lik);

            // output cluster labels
            vector<unsigned int> cluster_labels;
            Node<BulkDatum,CloneTreeNodeParam>::get_cluster_labels(root_node,
                                                                   data,
                                                                   cluster_labels);
            ofstream f;
            f.open(output_path + "/cluster_labels.txt", ios::out);
            for (size_t i = 0; i < cluster_labels.size(); i++) {
                f << data[i]->GetId() << "," << cluster_labels[i] << "\n";
            }
            f.close();
        }

        gsl_rng_free(random);

        // copy the simul.config to simul.output_path
        boost::filesystem::copy_file(config_file_path, simul_config.output_path + "/simul.config", boost::filesystem::copy_option::overwrite_if_exists);
    }

    cout << "Done!" << endl;

    return 0;
}
