//
//  bulk_datum.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-06-27.
//

#include "bulk_datum.hpp"

#include <iostream>

BulkDatum::BulkDatum(string name, Locus &locus) :
id_(name), locus_(locus)
{
}

BulkDatum::BulkDatum(string name,
                     Locus &locus,
                     size_t n_variants,
                     size_t n_reads) :
id_(name),
locus_(locus),
read_count_(n_reads),
variant_read_count_(n_variants)
{
}
BulkDatum::BulkDatum(string name,
                     Locus &locus,
                     size_t n_variants,
                     size_t n_reads,
                     size_t total_cn) :
id_(name),
locus_(locus),
read_count_(n_reads),
variant_read_count_(n_variants),
total_cn_(total_cn)
{
}

BulkDatum::BulkDatum(string name, Locus &locus, size_t n_variants, size_t n_reads,
                     size_t major_cn, size_t minor_cn) :
id_(name),
locus_(locus),
read_count_(n_reads),
variant_read_count_(n_variants)
{
    SetGenotype(make_pair(major_cn, minor_cn));
}

void BulkDatum::SetTotalCopyNumber(size_t val)
{
    this->total_cn_ = val;
}
void BulkDatum::SetVariantReadCount(size_t val)
{
    this->variant_read_count_ = val;
}
void BulkDatum::SetReadCount(size_t val)
{
    this->read_count_ = val;
}
void BulkDatum::SetGenotype(pair<size_t, size_t> genotype)
{
    if (genotype.second > genotype.first) {
        std::cerr << "Error: invalid genotype passed in.\n";
        std::cerr << "MajorCN: " << genotype.first << ".\n";
        std::cerr << "MinorCN: " << genotype.second << ".\n";
        exit(-1);
    }
    genotype_.first = genotype.first;
    genotype_.second = genotype.second;
}
void BulkDatum::SetCopyNumberPrior(vector<double> &cn_profile)
{
    cn_profile_ = cn_profile;
}
