//
//  single_cell_dna.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-08-01.
//

#include "single_cell.hpp"

#include <iostream>

// SingleCellData class definitions
SingleCellData::SingleCellData(string cell_name, size_t loci_count) :
        cell_name(cell_name), var_reads_(loci_count), total_reads_(loci_count) {
}

SingleCellData::SingleCellData(string cell_name,
                               vector<size_t> &var_reads,
                               vector<size_t> &total_reads) :
        cell_name(cell_name), var_reads_(var_reads), total_reads_(total_reads) {
    size_t loci_count = total_reads.size();
    for (size_t i = 0; i < loci_count; i++) {
        if (total_reads[i] > 0) {
            loci_idxs_.push_back(i);
        }
    }
}

SingleCellData::~SingleCellData() {
}

string SingleCellData::GetName() const {
    return cell_name;
}

string SingleCellData::Print() const {
    string ret = cell_name + "\n";
    return ret;
}

void SingleCellData::InsertDatum(size_t loci_idx,
                                 size_t var_reads,
                                 size_t total_reads) {
    loci_idxs_.push_back(loci_idx);
    var_reads_[loci_idx] = var_reads;
    total_reads_[loci_idx] = total_reads;
}

const vector<size_t> &SingleCellData::getExprReads() const {
    return expr_reads;
}

void SingleCellData::setExprReads(const vector<size_t> &exprReads) {
    expr_reads = exprReads;
}

const vector<double> &SingleCellData::getZeroInflationProbs() const {
    return zero_inflation_probs;
}

void SingleCellData::setZeroInflationProbs(const vector<double> &zeroInflationProbs) {
    zero_inflation_probs = zeroInflationProbs;
}

double SingleCellData::getDepthSize() const {
    return depth_size;
}

void SingleCellData::setDepthSize(double depthSize) {
    depth_size = depthSize;
}

void SingleCellData::print() {
    // TODO extend with single cell allelic imbalance data
    cout << "[sc] cell " << cell_name << ", sn " << depth_size << endl;
    cout << "- expr: [";
    for (auto it: expr_reads) {
        cout << it << ", ";
    }
    cout << "]" << endl << "- ZI probs: [";
    for (auto it: zero_inflation_probs) {
        cout << it << ", ";
    }
    cout << "]" << endl << "- loci: [";
    for (auto it: loci_idxs_) {
        cout << it << ", ";
    }
    cout << "]" << endl << "- var|ref: [";
    for (int i = 0; i < var_reads_.size(); ++i) {
        cout << var_reads_[i] << "|" << total_reads_[i] << ", ";
    }
    cout << "]" << endl;
}

