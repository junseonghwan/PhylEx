//
//  loci.hpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-08-03.
//

#ifndef loci_hpp
#define loci_hpp

#include <string>
#include <vector>

using namespace std;

enum PairedReadType
{
    VAR_VAR = 0,
    VAR_REF = 1,
    REF_VAR = 2,
    REF_REF = 3,
    size = 4
};

enum VariantsPhased
{
    NO = 0,
    YES = 1,
    UNKNOWN = 2
};

class Locus
{
    string chr = "";
    size_t pos = 0;
    string gene_name = "";
    string mutation_id;
    string unique_str;
    size_t hash_code;
    // Hyper parameter for Beta prior on probability of for bi-allelic locus.
    double alpha_ = 1.0;
    double beta_ = 1.0;
    double dropout_prob_ = 0.5;
public:
    Locus(string mutation_id, string chr, size_t pos, string gene_name);
    Locus(const Locus &other);
    inline string get_chr() const { return chr; }
    inline size_t get_pos() const { return pos; }
    inline size_t get_hash_code() const { return hash_code; }
    inline string get_mutation_id() const { return mutation_id; }
    inline string get_gene_name() const { return gene_name; }
    bool operator==(const Locus &other) const;
    string to_str() const;
    double get_alpha() const;
    double get_beta() const;
    double get_dropout_prob() const { return dropout_prob_; }
    void set_alpha(double alpha);
    void set_beta(double beta);
    void set_dropout_prob(double dropout_prob) {
        dropout_prob_ = dropout_prob;
    }
};

class LocusDatum
{
    size_t n_var_reads;
    size_t n_total_reads;
public:
    LocusDatum(size_t n_total_reads, size_t n_var_reads);
    inline size_t get_n_total_reads() const { return n_total_reads; }
    inline size_t get_n_var_reads() const { return n_var_reads; }
};

class LociPair
{
    Locus germline_locus;
    Locus somatic_locus;
    VariantsPhased var_phased;

public:
    LociPair(const Locus &germline_locus, const Locus &somatic_locus, VariantsPhased var_phased);
    inline const Locus &get_germline_locus() const { return germline_locus; }
    inline const Locus &get_somatic_locus() const { return somatic_locus; }
    inline const VariantsPhased is_var_phased() const { return var_phased; }
    size_t distance() const;
    bool operator==(const LociPair &other) const;
};

// stores reads that cover both loci
class LociPairDatum
{
    LociPair *loci_pair;
    vector<size_t> reads;
public:
    LociPairDatum(LociPair *loci_pair, vector<size_t> reads);
    const size_t get_n_reads(PairedReadType read_type) const;
    const size_t get_total_reads() const;
    const LociPair *get_loci_pair() const;
};

namespace std
{
    template <>
    struct hash<Locus> {
    public:
        size_t operator()(const Locus& obj) const {
            return obj.get_hash_code();
        }
    };
}

#endif /* loci_hpp */
