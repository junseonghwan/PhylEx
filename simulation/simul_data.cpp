//
//  simul_data.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-08-12.
//

#include "simul_data.hpp"

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

#include "data_util.hpp"
#include "utils.hpp"
#include "tssb_state.hpp"

size_t SampleCnProfile(const gsl_rng *random,
                       size_t curr,
                       Eigen::Ref<Eigen::MatrixXf> P)
{
    if (curr == 0) {
        return curr;
    }
    
    // sample a new state
    vector<double> probs;
    for (size_t j = 0; j < P.cols(); j++) {
        probs.push_back(P(curr, j));
    }
    size_t new_state = multinomial(random, probs);
    return new_state;
}

// Returns matrix Q of dimension NxN where N = max_cn - min_cn.
// If min_cn = 0, then since 0 is an absorbing state, Q(0,j) = 0.
Eigen::MatrixXf GetCnRateMatrix(double birth_rate,
                                double death_rate,
                                size_t min_cn,
                                size_t max_cn)
{
    size_t n_states = max_cn - min_cn + 1;
    Eigen::MatrixXf Q = Eigen::MatrixXf::Zero(n_states, n_states);
    for (size_t i = 0; i < n_states; i++) {
        if (i >= 1) {
            Q(i,i-1) = death_rate;
        }
        if (i < n_states - 1) {
            Q(i,i+1) = birth_rate;
        }
    }
    Q.diagonal() = -Q.rowwise().sum();
    
    // We check if cn state 0 is in the rate matrix.
    // If so, set it as the absorbing state.
    if (min_cn == 0) {
        Q.row(0).setZero();
    }

    assert(abs(Q.sum()) < 1e-6);
    return Q;
}

Eigen::MatrixXf ExponentiateMatrix(Eigen::Ref<Eigen::MatrixXf> M)
{
    std::cout << M << std::endl;
    Eigen::EigenSolver<Eigen::MatrixXf> diagnoalization;
    diagnoalization.compute(M);
    auto evalues = diagnoalization.eigenvalues();
    Eigen::MatrixXcf D = Eigen::MatrixXcf::Zero(evalues.rows(), evalues.rows());
    for (size_t i = 0; i < evalues.rows(); i++) {
        D(i,i) = exp(evalues(i).real());
    }
    Eigen::MatrixXcf V = diagnoalization.eigenvectors();
    Eigen::MatrixXcf ret = V * D * V.inverse();
    Eigen::MatrixXf expM = Eigen::MatrixXf::Zero(M.rows(), M.cols());
    for (size_t i = 0; i < expM.rows(); i++) {
        for (size_t j = 0; j < expM.cols(); j++) {
            if (ret.coeffRef(i, j).imag() > 0) {
                cerr << "Imaginary number found in the transition matrix." << endl;
                cerr << "(" << ret.coeffRef(i, j).imag() << ", " << ret.coeffRef(i, j).real() << ")" << std::endl;
                //exit(-1);
            }
            expM(i,j) = ret(i, j).real();
        }
    }
    auto rowsums = expM.rowwise().sum();
    for (size_t i = 0; i < rowsums.size(); i++) {
        assert(abs(1 - rowsums(i)) < 1e-3);
    }

    return expM;
}

vector<Locus> CreateSNVs(gsl_rng *random,
                         const SimulationConfig &simul_config,
                         vector<BulkDatum *> &data)
{
    string chr;
    size_t pos;
    vector<Locus> loci;
    for (size_t i = 0; i < simul_config.n_sites; i++)
    {
        chr = convert_chr_to_string(gsl_rng_uniform_int(random, 23) + 1);
        pos = gsl_rng_uniform_int(random, 10000000);
        size_t alpha = gsl_rng_uniform_int(random, simul_config.beta_binomial_hp_max - 1) + 1;
        size_t beta = gsl_rng_uniform_int(random, simul_config.beta_binomial_hp_max - 1) + 1;
        BulkDatum *datum = new BulkDatum("s" + to_string(i), chr, pos);
        datum->SetLocuHyperParameters(alpha, beta, 0.5);
        data.push_back(datum);
    }
    return loci;
}

CloneTreeNode *SampleFromTssbPrior(size_t region_count,
                                   const gsl_rng *random,
                                   size_t n_data,
                                   ModelParams &model_params,
                                   vector<BulkDatum *> &data)
{
    CloneTreeNode *root = CloneTreeNode::CreateRootNode(region_count);
    root->SetNuStick(0.0);
    root->SampleNodeParameters(random, model_params, 0);

    // sample tree and latent assignment of datum to node
    for (size_t i = 0; i < n_data; i++) {
        double u = gsl_ran_flat(random, 0, 1);
        auto node = CloneTreeNode::FindNode(random, u, root, model_params);
        node->AddDatum(data[i]);
    }

    cull(root);

    return root;
}

void GenerateBulkData(gsl_rng *random,
                      const SimulationConfig &simul_config,
                      vector<BulkDatum *> &data,
                      CloneTreeNode *root_node)
{
    // Sample a node for each SNV.
    // When assigning SNV to a node, sample clonal copy number profile.
    size_t n_data = data.size();
    vector<CloneTreeNode *> all_nodes;
    CloneTreeNode::BreadthFirstTraversal(root_node, all_nodes, false);
    double seq_err = simul_config.seq_err;
    size_t b_alleles, depth;
    BulkDatum *datum;
    for (size_t i = 0; i < n_data; i++) {
        // Sample a node -- +1 so that we don't sample the root node
        size_t node_id = discrete_uniform(random, all_nodes.size()-1) + 1;
        auto node = all_nodes[node_id];
        datum = data[i];
        node->AddDatum(datum);
        
        // Generate bulk data for each region.
        for (size_t region = 0; region < simul_config.n_regions; region++) {
            double xi = 0.0;
            double total_cn = 2;
            size_t total_allele_count, var_allele_count, ref_allele_count, var_cn;

            if (node->GetParentNode() == 0) {
                // This is an error b/c,
                // We ensured when sampling node_id to not choose the root.
                cerr << "Error: we sampled the root node." << endl;
                exit(-1);
            } else {
                // Sample clonal copy number.
                var_allele_count = discrete(random, simul_config.var_allele_copy_prob);
                ref_allele_count = discrete(random, simul_config.ref_allele_copy_prob);
                total_allele_count = var_allele_count + ref_allele_count;
                double phi = node->NodeParameter().GetCellularPrevalences()[region];
                total_cn = phi * total_allele_count + (1 - phi) * 2;
                depth = gsl_ran_poisson(random, simul_config.bulk_mean_depth * total_cn/2);
                // Ensure var_cn >= 1.
                var_cn = gsl_ran_binomial(random, simul_config.var_cp_prob, var_allele_count - 1) + 1;
                if (var_cn == total_allele_count) {
                    xi = (1 - phi) * seq_err + phi * (1 - seq_err);
                } else {
                    xi = (1 - phi) * seq_err + phi * (double)var_cn/total_allele_count;
                }
            }
            b_alleles = gsl_ran_binomial(random, xi, depth);
            if (ref_allele_count > var_allele_count) {
                datum->AddRegionData(b_alleles, depth, ref_allele_count, var_allele_count);
            } else {
                datum->AddRegionData(b_alleles, depth, var_allele_count, ref_allele_count);
            }
            cout << "======" << endl;
            cout << "Datum " << i << endl;
            cout << datum->GetVariantReadCount()[region] << "/" << datum->GetReadCount()[region] << endl;
            cout << "xi: " << xi << endl;
            cout << "======" << endl;
        }
    }
}

// Assumes nodes are in breadth first traversal order.
Eigen::MatrixXf EvolveCn(gsl_rng *random,
                     vector<CloneTreeNode *> &nodes,
                     unordered_map<CloneTreeNode *, size_t> &node2idx,
                     CloneTreeNode *assigned_node,
                     const SimulationConfig &config,
                     Eigen::MatrixXf P0,
                     Eigen::MatrixXf P1)
{
    size_t n_nodes = nodes.size();

    // Indexing:
    // Rows -> nodes.
    // First column -> ref_cn.
    // Second column -> var_cn.
    Eigen::MatrixXf cn_profile(n_nodes, 2);

    // Get ancestor nodes.
    unordered_set<CloneTreeNode *> ancestors;
    auto node = assigned_node;
    while (true) {
        if (node->IsRoot()) {
            break;
        }
        ancestors.insert(node);
        node = node->GetParentNode();
    }

    size_t ref_cn = 0, var_cn = 0;
    for (size_t i = 0; i < nodes.size(); i++) {
        auto node = nodes[i];
        auto parent_node = node->GetParentNode();
        if (node->IsRoot()) {
            cn_profile(i, 0) = 2;
            cn_profile(i, 1) = 0;
            continue;
        }

        // Retrieve the copy number profile of the parent
        auto parent_idx = node2idx[parent_node];
        ref_cn = cn_profile(parent_idx, 0);
        var_cn = cn_profile(parent_idx, 1);
        if (node == assigned_node) {
            assert(var_cn == 0 && ref_cn >= 1);
            // Evolve ref_cn using P1.
            SampleCnProfile(random, ref_cn, P1);
            // Sample var_cn, ensure that it is at least 1.
            var_cn = gsl_ran_binomial(random, config.var_cp_prob, ref_cn - 1) + 1;
            assert(var_cn >= 1);
            ref_cn -= var_cn;
        } else if (ancestors.count(node)) {
            ref_cn = SampleCnProfile(random, ref_cn, P1);
            assert(var_cn == 0 && ref_cn >= 1);
        } else {
            ref_cn = SampleCnProfile(random, ref_cn, P0);
            var_cn = SampleCnProfile(random, var_cn, P0);
        }
        assert(ref_cn <= config.max_cn && var_cn <= config.max_cn);

        cn_profile(i, 0) = ref_cn;
        cn_profile(i, 1) = var_cn;
    }

    return cn_profile;
}

void GenerateBulkDataWithBDProcess(gsl_rng *random,
                                   const SimulationConfig &simul_config,
                                   vector<BulkDatum *> &data,
                                   CloneTreeNode *root_node)
{
    unordered_map<CloneTreeNode*, vector<pair<size_t, size_t> > > cn_profile;

    size_t n_data = data.size();
    double birth_rate = simul_config.birth_rate;
    double death_rate = simul_config.death_rate;

    // Construct rate matrix for copy number profile and perform uniformization.
    Eigen::MatrixXf Q1 = GetCnRateMatrix(birth_rate, death_rate, 1, simul_config.max_cn);
    Eigen::MatrixXf P1 = ExponentiateMatrix(Q1);

    Eigen::MatrixXf Q0 = GetCnRateMatrix(birth_rate, death_rate, 0, simul_config.max_cn);
    Eigen::MatrixXf P0 = ExponentiateMatrix(Q0);

    vector<CloneTreeNode *> nodes;
    CloneTreeNode::BreadthFirstTraversal(root_node, nodes, false);
    unordered_map<CloneTreeNode *, size_t> node2idx;
    for (size_t i = 0; i < nodes.size(); i++) {
        auto node = nodes[i];
        node2idx[node] = i;
    }

    // Assign data to nodes.
    size_t b_alleles, depth;
    BulkDatum *datum;
    for (size_t i = 0; i < n_data; i++) {
        // Sample a node -- +1 so that we don't sample the root node
        size_t node_id = discrete_uniform(random, nodes.size()-1) + 1;
        auto assigned_node = nodes[node_id];
        datum = data[i];
        assigned_node->AddDatum(datum);

        // Evolve copy number.
        auto cn_profile = EvolveCn(random,
                                   nodes,
                                   node2idx,
                                   assigned_node,
                                   simul_config,
                                   P0,
                                   P1);

        auto node_cns = cn_profile.rowwise().sum();
        double total_cn = 0.0;
        double variant_cn = 0.0;
        double reference_cn = 0.0;
        double xi = 0.0;
        double sum = 0.0;
        // TODO: Incorporate sequencing error.
        for (size_t j = 0; j < nodes.size(); j++) {
            double eta = nodes[j]->GetCloneFreqs(0);
            total_cn += eta * node_cns(j);
            reference_cn += eta * cn_profile(j,0);
            variant_cn += eta * cn_profile(j,1);
            if (node_cns(j) > 0) {
                xi += eta * (double)cn_profile(j,1)/node_cns(j);
            }
            sum += eta;
        }
        depth = gsl_ran_poisson(random, simul_config.bulk_mean_depth * total_cn/2);
        cout << xi << ", " << depth << endl;
        b_alleles = gsl_ran_binomial(random, xi, depth);
        size_t int_var_cn = (size_t)round(variant_cn);
        size_t int_ref_cn = (size_t)round(reference_cn);
        if (int_var_cn > int_ref_cn) {
            datum->AddRegionData(b_alleles, depth, int_var_cn, int_ref_cn);
        } else {
            datum->AddRegionData(b_alleles, depth, int_ref_cn, int_var_cn);
        }
        cout << "======" << endl;
        cout << "Datum " << i << endl;
        cout << datum->GetVariantReadCount(0) << "/" << datum->GetReadCount(0) << endl;
        cout << "xi: " << xi << endl;
        cout << "======" << endl;
    }
}

vector<CloneTreeNode *> GenerateScRnaData(gsl_rng *random,
                                          CloneTreeNode *root_node,
                                          const vector<BulkDatum *> &data,
                                          const ModelParams &model_params,
                                          const SimulationConfig &simul_config,
                                          vector<SingleCellData *> &sc_data)
{
    vector<size_t> bulk_sc_coverage(data.size(), 0);
    // We subselect bulk to have single cell coverage.
    for (size_t i = 0; i < data.size(); i++) {
        if (gsl_ran_flat(random, 0, 1) < simul_config.snv_sc_sparsity) {
            bulk_sc_coverage[i] = 1;
        }
    }

    // Retrieve all nodes/clones with at least one SNV assigned.
    vector<CloneTreeNode *> non_empty_nodes;
    CloneTreeNode::BreadthFirstTraversal(root_node, non_empty_nodes, true);

    vector<CloneTreeNode *> cell2node;
    CloneTreeNode *node;
    for (size_t c = 0; c < simul_config.n_cells; c++) {
        size_t idx = gsl_rng_uniform_int(random, non_empty_nodes.size());
        node = non_empty_nodes[idx];
        cell2node.push_back(node);

        cout << "Cell " << c << " assigned to " << node->GetName() << endl;
        SingleCellData *sc = new SingleCellData("c" + to_string(c), data.size());
        GenerateScRnaReads(random,
                            simul_config,
                            node,
                            data,
                            bulk_sc_coverage,
                            *sc);
        
        sc_data.push_back(sc);
        cout << "=====" << endl;
    }
    
    return cell2node;
}

void GenerateScRnaReads(const gsl_rng *random,
                        const SimulationConfig &simul_config,
                        CloneTreeNode *node,
                        const vector<BulkDatum*> &somatic_loci,
                        vector<size_t> &bulk_sc_coverage,
                        SingleCellData &sc)
{
    size_t n_expr_loci = 0;
    if (simul_config.randomize_dropout) {
        // In this simulation, we randomly select number of sites to be
        // expressed.
        double realized_dropout_rate = gsl_ran_flat(random, simul_config.dropout_rate, 1);
        n_expr_loci = (size_t) round(somatic_loci.size() * (1 - realized_dropout_rate));
    } else {
        // In this simulation, we fix the number of sites to be
        // expressed.
        n_expr_loci = (size_t)round(somatic_loci.size() * (1 - simul_config.dropout_rate));
    }

    // Randomly select n_expr_loci to be expressed.
    unordered_set<size_t> sampled_loci_idx;
    size_t locus_idx = 0;
    for (size_t i = 0; i < n_expr_loci; i++) {
        do {
            locus_idx = discrete_uniform(random, somatic_loci.size());
        } while(sampled_loci_idx.count(locus_idx) != 0);
        if (bulk_sc_coverage[locus_idx] == 1) {
            sampled_loci_idx.insert(locus_idx);
        }
    }
    
    cout << "Number of sites expressed: " << sampled_loci_idx.size() << "\n";

    unordered_set<Locus> snvs;
    CloneTreeNode::RetrieveLoci(node, snvs);
    size_t var_loci_expr_count = 0;
    size_t var_read_observed_count = 0;
    for (auto locus_idx : sampled_loci_idx) {
        size_t total_read_count = gsl_ran_poisson(random, simul_config.sc_mean_depth);
        size_t var_read_count = 0;
        auto locus = somatic_loci[locus_idx]->GetLocus();
        bool carries_snv = (snvs.count(locus) > 0);
        double var_expr_prob;
        if (carries_snv) {
            var_loci_expr_count++;
            // Sample the number of variants by first sampling either bursty or non-bursty distn.
            bool bursty = gsl_ran_bernoulli(random, simul_config.bursty_prob);
            if (bursty) {
                //cout << "Bursty." << endl;
                var_expr_prob = gsl_ran_beta(random,
                                             simul_config.sc_bursty_alpha0,
                                             simul_config.sc_bursty_beta0);
            } else {
                //cout << "Not Bursty." << endl;
                double alpha = locus.get_alpha();
                double beta = locus.get_beta();
                var_expr_prob = gsl_ran_beta(random, alpha, beta);
            }
        } else {
            // We expect to see variant read only in error.
            //cout << "SNV absent." << endl;
            var_expr_prob = gsl_ran_beta(random, simul_config.seq_err, 1 - simul_config.seq_err);
        }
        var_read_count = gsl_ran_binomial(random, var_expr_prob, total_read_count);
        var_read_observed_count += var_read_count > 0 ? 1 : 0;
        sc.InsertDatum(locus_idx, var_read_count, total_read_count);
    }
    cout << var_loci_expr_count << "/" << snvs.size() << " variants expressed.\n";
    cout << var_read_observed_count << "/" << var_loci_expr_count << " sites with variant reads observed.\n";
}

