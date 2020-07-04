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

enum SC_DATA_TYPE {
    DNA = 0, RNA = 1
};

// Single cell data has reads at each site of interest, this is common to both DNA and RNA data
// if DNA uses mutation matrix (0/1/3) then, reads will be set to 0 (null pointer)
// for RNA, in addition to these reads, it may have reads that cover germline and somatic loci
class SingleCellData
{
    string cell_name;
    unordered_map<Locus, LocusDatum*> reads;
    unordered_map<Locus, size_t> mutation_map;
    unordered_map<Locus, double> expression_levels;
    SC_DATA_TYPE sc_data_type;

public:
    SingleCellData(string cell_name);
    SingleCellData(string cell_name, unordered_map<Locus, size_t> &mutation_map);
    SingleCellData(string cell_name, unordered_map<Locus, size_t> &mutation_map, unordered_map<Locus, LocusDatum*> &reads);
    SingleCellData(string cell_name,
                  unordered_map<Locus, LocusDatum*> &reads);

    string get_name() const;
    inline SC_DATA_TYPE get_data_type() const { return sc_data_type; }
    inline void set_data_type(SC_DATA_TYPE val) { sc_data_type = val; }

    const LocusDatum *get_locus_datum(const Locus &locus) const;
    const bool has_locus_datum(const Locus &locus) const;
    const LociPairDatum *get_loci_pair_datum(const Locus &locus) const;
    double get_expression_level(const Locus &locus) const;

    size_t get_mutation(const Locus &locus) const;

    void set_expr_level(Locus&, double);
    void insert_read(const Locus &, LocusDatum *);
    void remove_locus(const Locus &locus);
    
    string print() const;

    virtual ~SingleCellData();
};

#endif /* single_cell_dna_hpp */
