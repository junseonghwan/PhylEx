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
                               unordered_map<Locus, size_t> &mutation_map) :
cell_name(cell_name), mutation_map(mutation_map)
{
    sc_data_type = SC_DATA_TYPE::DNA;
}

SingleCellData::SingleCellData(string cell_name,
                               unordered_map<Locus, size_t> &mutation_map,
                               unordered_map<Locus, LocusDatum*> &reads) :
cell_name(cell_name), reads(reads), mutation_map(mutation_map)
{
    sc_data_type = SC_DATA_TYPE::DNA;
}

SingleCellData::SingleCellData(string cell_name,
                               unordered_map<Locus, LocusDatum*> &reads) :
cell_name(cell_name), reads(reads)
{
    sc_data_type = SC_DATA_TYPE::RNA;
}

SingleCellData::~SingleCellData()
{
}

string SingleCellData::get_name() const
{
    return cell_name;
}

//const unordered_map<Locus, LocusDatum*> &SingleCellData::get_reads() const
//{
//    return reads;
//}

void SingleCellData::insert_read(const Locus &locus, LocusDatum *locus_datum)
{
    if (reads.count(locus) == 0)
        reads[locus] = locus_datum;
}

size_t SingleCellData::get_mutation(const Locus &locus) const
{
    if (mutation_map.count(locus) == 0) {
        cerr << "Error: locus does not exist." << endl;
        exit(-1);
    }
    return mutation_map.at(locus);
}

string SingleCellData::print() const
{
    string ret = cell_name + "\n";
    return ret;
}

//const unordered_map<Locus, LociPairDatum*> &SingleCellData::get_paired_reads() const
//{
//    return paired_reads;
//}


//const unordered_map<Locus, double> &SingleCellData::get_expression_levels() const
//{
//    return expression_levels;
//}
//
//const unordered_map<Locus, size_t> &SingleCellData::get_mutation_map() const
//{
//    return mutation_map;
//}

//void SingleCellData::set_paired_read(Locus &locus, LociPairDatum *loci_pair_datum)
//{
//    paired_reads[locus] = loci_pair_datum;
//}

void SingleCellData::set_expr_level(Locus &locus, double expr_level)
{
    expression_levels[locus] = expr_level;
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

//const LociPairDatum *SingleCellData::get_loci_pair_datum(const Locus &locus) const
//{
//    if (paired_reads.count(locus) == 0) {
//        return 0;
//    } return paired_reads.at(locus);
//}

double SingleCellData::get_expression_level(const Locus &locus) const
{
    if (expression_levels.count(locus) == 0) {
        return -1;
    }
    return expression_levels.at(locus);
}
