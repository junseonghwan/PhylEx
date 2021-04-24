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

const string &SingleCellExpression::getCellName() const {
    return cell_name;
}

void SingleCellExpression::setCellName(const string &cellName) {
    cell_name = cellName;
}

const vector<size_t> &SingleCellExpression::getGeneIdxs() const {
    return gene_idxs;
}

void SingleCellExpression::setGeneIdxs(const vector<size_t> &geneIdxs) {
    gene_idxs = geneIdxs;
}

const vector<size_t> &SingleCellExpression::getExprReads() const {
    return expr_reads;
}

void SingleCellExpression::setExprReads(const vector<size_t> &exprReads) {
    expr_reads = exprReads;
}

SingleCellExpression::SingleCellExpression(string cellName, double depth_size) : cell_name(std::move(cellName)),
                                                                                 depth_size(depth_size) {}

const vector<double> &SingleCellExpression::getZeroInflationProbs() const {
    return zero_inflation_probs;
}

void SingleCellExpression::setZeroInflationProbs(const vector<double> &zeroInflationProbs) {
    zero_inflation_probs = zeroInflationProbs;
}

SingleCellExpression::SingleCellExpression(string cellName) : cell_name(std::move(cellName)) {}

double SingleCellExpression::getDepthSize() const {
    return depth_size;
}

void SingleCellExpression::setDepthSize(double depthSize) {
    depth_size = depthSize;
}

void SingleCellExpression::print() {
    cout << "[sce] cell " << cell_name << ", sn " << depth_size << endl;
    cout << "- expr: [";
    for (auto it: expr_reads) {
        cout << it << ", ";
    }
    cout << "]" << endl << "- ZI probs: [";
    for (auto it: zero_inflation_probs) {
        cout << it << ", ";
    }
    cout << "]" << endl;
}

