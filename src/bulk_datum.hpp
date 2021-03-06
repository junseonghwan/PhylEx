//
//  bulk_datum.hpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-06-27.
//

#ifndef bulk_datum_hpp
#define bulk_datum_hpp

#include <string>
#include <unordered_map>
#include <vector>

#include "loci.hpp"

using namespace std;

class BulkDatum
{
    string  id_;

    Locus locus_;
    
    vector<size_t> variant_reads_, total_reads_;
    vector<size_t> major_cns_, minor_cns_;
    vector<size_t> total_cns_;
    // Each entry stores prior belief on the total copy number.
    // e.g., cn_prob_matrix[0] is the prior belief on total copy number for the
    // first region. cn_prob_matrix[0][0] is the probability that total copy
    // number is 0, cn_prob_matrix[0][n], where n is the total possible copy number.
    vector<vector<double> > cn_prob_matrix_;
public:
    BulkDatum(string name, string chr, size_t pos);
    BulkDatum(string name, string chr, size_t pos,
              vector<size_t> &n_variants, vector<size_t> &n_reads);
    BulkDatum(string name, string chr, size_t pos,
              vector<size_t> &n_variants, vector<size_t> &n_reads,
              vector<size_t> &total_cns);
    BulkDatum(string name, string chr, size_t pos,
              vector<size_t> &n_variants, vector<size_t> &n_reads,
              vector<size_t> &major_cn, vector<size_t> &minor_cn);

    string GetId() const { return id_; }
    inline const vector<size_t> &GetVariantReadCount() const { return variant_reads_; }
    inline const vector<size_t> &GetReadCount() const { return total_reads_; }
    inline const vector<size_t> &GetMajorCN() const { return major_cns_; }
    inline const vector<size_t> &GetMinorCN() const { return minor_cns_; }
    inline const vector<size_t> &GetTotalCN() const { return total_cns_; }
    
    inline size_t GetVariantReadCount(size_t region) const { return variant_reads_[region]; }
    inline size_t GetReadCount(size_t region) const { return total_reads_[region]; }
    inline const size_t GetMajorCN(size_t region) const { return major_cns_[region]; }
    inline const size_t GetMinorCN(size_t region) const { return minor_cns_[region]; }
    inline const size_t GetTotalCN(size_t region) const { return total_cns_[region]; }

    const vector<double> &GetCNProbs(size_t idx) const { return cn_prob_matrix_[idx]; }
    const Locus &GetLocus() const { return locus_; }

    void AddRegionData(size_t var_reads, size_t total_reads,
                       size_t major_cn, size_t minor_cn);
    void SetCopyNumberProbs(vector<vector<double> > cn_prob_matrix) {
        cn_prob_matrix_ = cn_prob_matrix;
    }
    void SetLocuHyperParameters(double alpha, double beta, double delta0);
    size_t GetRegionCount() const {
        return variant_reads_.size();
    }
};

#endif /* bulk_datum_hpp */
