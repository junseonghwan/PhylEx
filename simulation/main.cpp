
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <string>

#include "bulk_datum.hpp"
#include "clone_node.hpp"
#include "data_util.hpp"
#include "model_params.hpp"
#include "single_cell.hpp"
#include "simul_config.hpp"
#include "tssb_state.hpp"

#include "simul_data.hpp"

using namespace std;

/**
 * Parse configuration file. The format is `key: value` for each line
 * and must contain all the necessary parameters.
 *
 * @param[in] config_file_path where the configuration file is located
 * @param[out] config object saving the parameters
 */
void parse_config_file(const string &config_file_path, SimulationConfig &config) {
    string line;
    ifstream config_file(config_file_path);
    if (!config_file.is_open()) {
        cerr << "Could not open the file: " << config_file_path << endl;
        exit(-1);
    }

    vector<string> results;

    while (getline(config_file, line)) {
        boost::split(results, line, boost::is_any_of(":"));
        boost::algorithm::trim(results[1]);
        config.insert_option(results[0], results[1]);
    }
    config_file.close();
}

void CreateLinearTree(size_t region_count,
                      gsl_rng *random,
                      CloneTreeNode *root,
                      size_t max_depth) {
    vector<CloneTreeNode *> nodes;
    nodes.push_back(root);
    auto node = root;
    for (size_t i = 0; i < max_depth; i++) {
        node = (CloneTreeNode *) node->SpawnChild(1);
        nodes.push_back(node);
    }
    size_t n_nodes = max_depth + 1;
    for (size_t region = 0; region < region_count; region++) {
        // Sample cellular prevalences.
        double *cell_prev = new double[n_nodes];
        double *clone_freq = new double[n_nodes];
        cell_prev[0] = 1.0;
        cell_prev[1] = gsl_ran_flat(random, 0.5, 1.0);
        for (size_t i = 2; i < n_nodes; i++) {
            cell_prev[i] = cell_prev[i - 1] / 2;
        }
        for (size_t i = 0; i < n_nodes; i++) {
            clone_freq[i] = cell_prev[i] - cell_prev[i + 1];
            nodes[i]->SetCellularPrevalenceAtRegion(region, cell_prev[i]);
            nodes[i]->SetCloneFrequencyAtRegion(region, clone_freq[i]);
        }
    }

    for (size_t i = 0; i < n_nodes; i++) {
        cout << nodes[i]->Print() << endl;
    }
}

bool BelowMinimumCellFraction(CloneTreeNode *node, double min) {
    auto cell_fractions = node->NodeParameter().GetCloneFreqs();
    for (auto cf : cell_fractions) {
        if (cf <= min) {
            return true;
        }
    }
    return false;
}

void CreateNaryTree(size_t region_count,
                    gsl_rng *random,
                    CloneTreeNode *root,
                    size_t max_depth,
                    size_t num_branches,
                    bool randomize_branching,
                    bool randomize_cf,
                    double min_cell_fraction) {
    deque<CloneTreeNode *> queue;
    queue.push_back(root);
    root->NodeParameter().SetRootParameters();

    while (!queue.empty()) {
        auto node = queue.front();
        queue.pop_front();
        if (node->GetDepth() >= max_depth ||
            BelowMinimumCellFraction(node, min_cell_fraction)) {
            // develop other branches in the queue
            continue;
        }

        // define number of children
        size_t branch_count;
        if (node->IsRoot()) {
            // root always has exactly one child
            branch_count = 1;
        } else {
            if (randomize_branching) {
                branch_count = gsl_rng_uniform_int(random, num_branches - 1) + 1;
            } else {
                branch_count = num_branches;
            }
        }

        for (size_t i = 0; i < branch_count; i++) {
            auto child_node = (CloneTreeNode *) node->SpawnChild(0.5);
            // Set cellular prevalence.
            for (size_t region = 0; region < region_count; region++) {
                double parent_clone_freq = node->NodeParameter().GetCloneFreqs()[region];
                double u = 0.5;
                if (randomize_cf) {
                    u = gsl_ran_flat(random, 0, 1);
                }
                double cell_prev = u * parent_clone_freq;
                // clone frequency is equal to cell prevalence for leaf nodes
                child_node->SetCloneFrequencyAtRegion(region, cell_prev);

                child_node->SetCellularPrevalenceAtRegion(region, cell_prev);
                node->SetCloneFrequencyAtRegion(region, parent_clone_freq - cell_prev);
            }
            queue.push_back(child_node);
        }
    }
    vector<CloneTreeNode *> nodes;
    CloneTreeNode::BreadthFirstTraversal(root, nodes);
    for (auto &node : nodes) {
        cout << node->Print() << endl;
    }
}

int main(int argc, char *argv[]) {
    string config_file_path;

    namespace po = boost::program_options;
    po::options_description desc("Program options");
    desc.add_options()
            ("help", "Put a help message here.")
            ("config,c", po::value<string>(&config_file_path)->required(), "path to config file.");

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

    ModelParams model_params;
    model_params.SetAlpha0Bound(true, 10.0);
    model_params.SetLambdaBound(false, 1.0);
    model_params.SetGammaBound(true, 10.0);
    model_params.SetSequencingError(simul_config.seq_err);
    model_params.RandomInit(rand);
    model_params.SetVariantCopyProbability(simul_config.var_cp_prob);
    model_params.SetSingleCellBurstyAlphaParameter(simul_config.sc_bursty_alpha0);
    model_params.SetSingleCellBurstyBetaParameter(simul_config.sc_bursty_beta0);

    bool bd_process = false;
    if (simul_config.birth_rate > 0 && simul_config.max_cn > 2) {
        bd_process = true;
    } else if (simul_config.ref_allele_copy_prob.size() >= 1 &&
               simul_config.var_allele_copy_prob.size() >= 2) {
        bd_process = false;
    } else {
        cerr << "Error in copy number simulation configuration." << endl;
        cerr << "For BD simulation, specify birth, death rates, and max_cn." << endl;
        cerr << "For clonal copy number, specify prior over copy number with at least 3 non-zero entries." << endl;
    }

    string output_path = simul_config.output_path;

    // Create directories for output.
    boost::filesystem::path outpath(output_path);
    boost::filesystem::create_directories(outpath);

    auto root_node = CloneTreeNode::CreateRootNode(simul_config.n_regions);
    root_node->SampleNodeParameters(rand, model_params);
    if (simul_config.num_branches == 1) {
        CreateLinearTree(simul_config.n_regions, rand, root_node, simul_config.max_depth);
    } else {
        CreateNaryTree(simul_config.n_regions,
                       rand,
                       root_node,
                       simul_config.max_depth,
                       simul_config.num_branches,
                       simul_config.randomize_branching,
                       simul_config.randomize_cf,
                       simul_config.min_cf);
    }

    vector<BulkDatum *> data;

    // Create somatic SNVs: chr and pos.
    // Each BulkDatum instance has a reference to the Locus instance.
    // We need to keep these alive in the memory during the simulation.
    // TODO: Consider a different approach? Shared pointer?
    // Note that the inference program depends on Locus being a reference
    // in the BulkDatum to update the hyper parameters for each locus from file.
    // So changing it to a copy of Locus for each BulkDatum is not an option.
    CreateSNVs(rand, simul_config, data);
    vector<pair<double, double> > cts_cn_profile;
    if (bd_process) {
        GenerateBulkDataWithBDProcess(rand, simul_config, data, root_node, cts_cn_profile);
    } else {
        GenerateBulkData(rand, simul_config, data, root_node);
    }

    WriteBulkData(output_path + "/simul_ssm.txt", data, false);
    WriteBulkData(output_path + "/genotype_ssm.txt", data, true);
    if (cts_cn_profile.size() > 0) {
        WriteCtsCopyNumberProfile(output_path + "/cts_cn_profile.txt",
                                  data,
                                  cts_cn_profile);
    }
    // generate or sample a gene set for expression data
    // TODO sample `gene_set` from real genes
    vector<Gene *> gene_set;
    GenerateGenes(rand, simul_config, gene_set);

    // Generate single cell reads.
    vector<SingleCellData *> sc_data;
    vector<SingleCellExpression *> sc_expr_data;
    auto cell2node = GenerateScRnaData(rand, root_node, data, model_params, simul_config,
                                       sc_data,
                                       sc_expr_data,
                                       gene_set);
    WriteScRnaData(output_path, data, sc_data);
    WriteCell2NodeAssignment(output_path, sc_data, cell2node);
    WriteBetaBinomHp(output_path, data);
    WriteScRnaExpressionData(output_path, sc_expr_data, gene_set);

    // Output information needed for evaluation.
    vector<CloneTreeNode *> sorted_nodes;
    CloneTreeNode::BreadthFirstTraversal(root_node,
                                         sorted_nodes,
                                         false);

    WriteClonalCNProfiles(output_path, sorted_nodes, gene_set);

    unordered_map<const BulkDatum *, CloneTreeNode *> datum2node;
    CloneTreeNode::Datum2Node(sorted_nodes, datum2node);
    double bulk_log_lik = 0.0;
    for (size_t region = 0; region < simul_config.n_regions; region++) {
        for (size_t i = 0; i < simul_config.n_sites; i++) {
            bulk_log_lik += BulkLogLikWithGenotype(region,
                                                   datum2node[data[i]],
                                                   data[i],
                                                   model_params);
        }
    }
    double sc_log_lik = ComputeSingleCellLikelihood(root_node,
                                                    data,
                                                    sc_data,
                                                    model_params);
    WriteTreeToFile(output_path, data, root_node);
    WriteLogLikToFile(output_path + "/log_lik_bulk.txt", bulk_log_lik);
    WriteLogLikToFile(output_path + "/log_lik_sc.txt", sc_log_lik);

    // output cluster labels
    vector<unsigned int> cluster_labels;
    CloneTreeNode::GetClusterLabels(root_node, data, cluster_labels);
    ofstream f;
    f.open(output_path + "/cluster_labels.txt", ios::out);
    for (size_t i = 0; i < cluster_labels.size(); i++) {
        f << data[i]->GetId() << "," << cluster_labels[i] << "\n";
    }
    f.close();

    gsl_rng_free(rand);
    cout << "Done!" << endl;

    return 0;
}
