//
//  single_cell_dna.hpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-08-01.
//

#ifndef single_cell_dna_hpp
#define single_cell_dna_hpp

#include <string>
#include <vector>
#include <unordered_map>

#include "loci.hpp"
#include "gene.hpp"

using namespace std;

// Single cell data has reads at each site of interest, this is common to both DNA and RNA data
// if DNA uses mutation matrix (0/1/3) then, reads will be set to 0 (null pointer)
// for RNA, in addition to these reads, it may have reads that cover germline and somatic loci
class SingleCellData
{
    string cell_name;

    // allelic imbalance data
    vector<size_t> var_reads_;
    vector<size_t> total_reads_;
    vector<size_t> loci_idxs_;

public:
    SingleCellData(string cell_name, size_t loci_count);
    SingleCellData(string cell_name,
                   vector<size_t> &var_reads,
                   vector<size_t> &total_reads);
    ~SingleCellData();

    string GetName() const;

    //const LocusDatum *get_locus_datum(const Locus &locus) const;
    //const bool has_locus_datum(const Locus &locus) const;

    inline size_t GetVariantReads(size_t loci_idx) const {
        return var_reads_[loci_idx];
    }
    inline size_t GetTotalReads(size_t loci_idx) const {
        return total_reads_[loci_idx];
    }
    inline const vector<size_t> &GetLoci() const {
        return loci_idxs_;
    }
    void InsertDatum(size_t loci_idx, size_t var_reads, size_t total_reads);
    string Print() const;
};

class SingleCellExpression {
    string cell_name;

    // gene expression data
    vector<size_t> gene_idxs; // maps expression read to specific gene
    vector<size_t> expr_reads;

    // parameters
    double depth_size; // s_c, or library size
    vector<double> zero_inflation_probs; // rho_gc

public:

    explicit SingleCellExpression(string cellName);

    SingleCellExpression(string cellName, double depth_size);

    const string &getCellName() const;

    void setCellName(const string &cellName);

    const vector<size_t> &getGeneIdxs() const;

    void setGeneIdxs(const vector<size_t> &geneIdxs);

    const vector<size_t> &getExprReads() const;

    void setExprReads(const vector<size_t> &exprReads);

    const vector<double> &getZeroInflationProbs() const;

    void setZeroInflationProbs(const vector<double> &zeroInflationProbs);

    double getDepthSize() const;

    void setDepthSize(double depthSize);

    void print();
};

#endif /* single_cell_dna_hpp */
