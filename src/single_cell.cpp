//
//  single_cell_dna.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-08-01.
//

#include "single_cell.hpp"

#include <iostream>

// SingleCellData class definitions
SingleCellData::SingleCellData(string cell_name) :
cell_name(cell_name)
{
}

SingleCellData::SingleCellData(string cell_name,
                               unordered_map<Locus, LocusDatum*> &reads) :
cell_name(cell_name), reads(reads)
{
}

SingleCellData::~SingleCellData()
{
}

string SingleCellData::get_name() const
{
    return cell_name;
}

void SingleCellData::insert_read(const Locus &locus, LocusDatum *locus_datum)
{
    if (reads.count(locus) == 0)
        reads[locus] = locus_datum;
}

string SingleCellData::print() const
{
    string ret = cell_name + "\n";
    return ret;
}

void SingleCellData::remove_locus(const Locus &locus)
{
    reads.erase(locus);
    mutation_map.erase(locus);
    //paired_reads.erase(locus);
    expression_levels.erase(locus);
}

const LocusDatum *SingleCellData::get_locus_datum(const Locus &locus) const
{
    if (reads.count(locus) == 0) {
        return 0;
    }
    return reads.at(locus);

}

const bool SingleCellData::has_locus_datum(const Locus &locus) const
{
    return reads.count(locus) > 0;
}
