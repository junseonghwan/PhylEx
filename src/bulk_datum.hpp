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

    Locus &locus_;

    size_t  read_count_ = 0;
    size_t  variant_read_count_ = 0;
    
    // Three different inputs for copy number.
    // 1. Total copy number.
    // 2. Major and minor copy numbers.
    // 3. Prior on total copy number.
    size_t  total_cn_ = 0;
    pair<size_t, size_t> genotype_;
    vector<double> cn_profile_;
public:
    BulkDatum(string name, Locus &locus);
    BulkDatum(string name, Locus &locus, size_t n_variants, size_t n_reads);
    BulkDatum(string name, Locus &locus, size_t n_variants, size_t n_reads,
              size_t total_cn);
    BulkDatum(string name, Locus &locus, size_t n_variants, size_t n_reads,
              size_t major_cn, size_t minor_cn);
    
    string GetId() const { return id_; }
    size_t GetVariantReadCount() const { return variant_read_count_; }
    size_t GetReadCount() const { return read_count_; }
    size_t GetTotalCN() const { return total_cn_; }
    size_t GetMajorCN() const { return genotype_.first; }
    size_t GetMinorCN() const { return genotype_.second; }
    const Locus &GetLocus() const { return locus_; }
    const vector<double> &GetCopyNumberProfile() const {
        return cn_profile_;
    }

    void SetVariantReadCount(size_t val);
    void SetReadCount(size_t n_reads);
    void SetTotalCopyNumber(size_t val);
    void SetGenotype(pair<size_t, size_t> val);
    void SetCopyNumberPrior(vector<double> &cn_profile);
};

#endif /* bulk_datum_hpp */
