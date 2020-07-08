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

BulkDatum::BulkDatum(string name, Locus &locus,
                     vector<size_t> &n_variants, vector<size_t> &n_reads) :
id_(name), locus_(locus),
variant_reads_(n_variants), total_reads_(n_reads)
{
    
}


BulkDatum::BulkDatum(string name, Locus &locus,
                     vector<size_t> &n_variants, vector<size_t> &n_reads,
                     vector<size_t> &total_cns) :
id_(name), locus_(locus),
variant_reads_(n_variants), total_reads_(n_reads),
total_cns_(total_cns)
{
    
}

BulkDatum::BulkDatum(string name, Locus &locus,
                     vector<size_t> &variant_reads, vector<size_t> &total_reads,
                     vector<size_t> &major_cns, vector<size_t> &minor_cns) :
id_(name),
locus_(locus),
variant_reads_(variant_reads),
total_reads_(total_reads),
major_cns_(major_cns),
minor_cns_(minor_cns)
{
}


void BulkDatum::AddRegionData(size_t var_reads, size_t total_reads,
                          size_t major_cn, size_t minor_cn)
{
    variant_reads_.push_back(var_reads);
    total_reads_.push_back(total_reads);
    major_cns_.push_back(major_cn);
    minor_cns_.push_back(minor_cn);
    total_cns_.push_back(major_cn + minor_cn);
}
