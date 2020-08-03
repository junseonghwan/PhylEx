//
//  loci.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-08-03.
//

#include "loci.hpp"

#include <iostream>
#include <cstdlib>
#include <limits.h>

LocusDatum::LocusDatum(size_t n_total_reads, size_t n_var_reads) :
n_var_reads(n_var_reads), n_total_reads(n_total_reads)
{
    if (n_var_reads > n_total_reads) {
        cerr << "Error: n_var_reads > n_total_reads.\n";
        cerr << "n_var_reads: " << n_var_reads << ".\n";
        cerr << "n_total_reads: " << n_total_reads << ".\n";
        exit(-1);
    }
}

Locus::Locus(string mutation_id, string chr, size_t pos) :
chr(chr), pos(pos), mutation_id(mutation_id)
{
    unique_str = chr + ":" + to_string(pos);
    hash_code = hash<string>()(mutation_id);
}

Locus::Locus(const Locus &other) :
chr(other.get_chr()),
pos(other.get_pos()),
mutation_id(other.mutation_id),
alpha_(other.alpha_),
beta_(other.beta_),
bursty_prob_(other.bursty_prob_)
{
    unique_str = chr + ":" + to_string(pos);
    //hash_code = hash<string>()(unique_str);
//    cout << other.get_chr() << ":" << other.get_pos() << endl;
    hash_code = hash<string>()(mutation_id);
}

string Locus::to_str() const
{
    return unique_str;
}

bool Locus::operator==(const Locus &other) const
{
    //return (this->chr == other.get_chr() && this->pos == other.get_pos());
    return (this->get_mutation_id() == other.get_mutation_id());
}

double Locus::get_alpha() const
{
    return alpha_;
}

double Locus::get_beta() const
{
    return beta_;
}

void Locus::set_alpha(double alpha)
{
    this->alpha_ = alpha;
}

void Locus::set_beta(double beta)
{
    this->beta_ = beta;
}

LociPair::LociPair(const Locus &germline_locus,
                   const Locus &somatic_locus,
                   VariantsPhased var_phased) :
germline_locus(germline_locus), somatic_locus(somatic_locus), var_phased(var_phased)
{
}

bool LociPair::operator==(const LociPair &other) const
{
    return (this->germline_locus == other.germline_locus && this->somatic_locus == other.somatic_locus);
}

size_t LociPair::distance() const
{
    if (this->germline_locus.get_chr() != this->somatic_locus.get_chr()) {
        return INT_MAX;
    }
    long dist = this->germline_locus.get_pos() - this->somatic_locus.get_pos();
    return abs(dist);
}

LociPairDatum::LociPairDatum(LociPair *loci_pair, vector<size_t> reads) :
loci_pair(loci_pair), reads(reads)
{
}

const size_t LociPairDatum::get_n_reads(PairedReadType read_type) const
{
    if (read_type >= PairedReadType::size) {
        //error
        cerr << "PairedReadType must be a value between 0-3." << endl;
        exit(-1);
    }
    return reads[read_type];
}

const size_t LociPairDatum::get_total_reads() const
{
    size_t n = reads[PairedReadType::VAR_VAR];
    n += reads[PairedReadType::VAR_REF];
    n += reads[PairedReadType::REF_VAR];
    n += reads[PairedReadType::REF_REF];
    return n;
}

const LociPair *LociPairDatum::get_loci_pair() const
{
    return loci_pair;
}
