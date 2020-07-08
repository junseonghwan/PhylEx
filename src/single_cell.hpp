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

using namespace std;

// Single cell data has reads at each site of interest, this is common to both DNA and RNA data
// if DNA uses mutation matrix (0/1/3) then, reads will be set to 0 (null pointer)
// for RNA, in addition to these reads, it may have reads that cover germline and somatic loci
class SingleCellData
{
    string cell_name;
    unordered_map<Locus, LocusDatum*> reads;
    unordered_map<Locus, size_t> mutation_map;
    unordered_map<Locus, double> expression_levels;

public:
    SingleCellData(string cell_name);
    SingleCellData(string cell_name,
                   unordered_map<Locus, LocusDatum*> &reads);

    string get_name() const;

    const LocusDatum *get_locus_datum(const Locus &locus) const;
    const bool has_locus_datum(const Locus &locus) const;

    void insert_read(const Locus &, LocusDatum *);
    void remove_locus(const Locus &locus);
    
    string print() const;

    virtual ~SingleCellData();
};

#endif /* single_cell_dna_hpp */
