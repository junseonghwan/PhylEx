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
cell_name(cell_name), var_reads_(loci_count), total_reads_(loci_count)
{
}

SingleCellData::SingleCellData(string cell_name,
                               vector<size_t> &var_reads,
                               vector<size_t> &total_reads) :
cell_name(cell_name), var_reads_(var_reads), total_reads_(total_reads)
{
    size_t loci_count = total_reads.size();
    for (size_t i = 0; i < loci_count; i++) {
        if (total_reads[i] > 0) {
            loci_idxs_.push_back(i);
        }
    }
}

SingleCellData::~SingleCellData()
{
}

string SingleCellData::GetName() const
{
    return cell_name;
}

string SingleCellData::Print() const
{
    string ret = cell_name + "\n";
    return ret;
}

void SingleCellData::InsertDatum(size_t loci_idx,
                                 size_t var_reads,
                                 size_t total_reads)
{
    loci_idxs_.push_back(loci_idx);
    var_reads_[loci_idx] = var_reads;
    total_reads_[loci_idx] = total_reads;
}
