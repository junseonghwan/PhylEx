//
// Created by zemp on 17/04/21.
//

#ifndef PHYLEX_GENE_HPP
#define PHYLEX_GENE_HPP

#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>

#include "bulk_datum.hpp"

using namespace std;

/**
 * Nucleic Acid Sequence: class to identify both genes and bins
 */
class NASequence {
protected:
    NASequence(const string &chr, size_t start_pos, size_t end_pos);
    NASequence(size_t numChr, size_t start_pos, size_t end_pos);

    string chr;
    size_t numChr;
    size_t start_pos;
    size_t end_pos;
public:
    bool operator<(const NASequence &rhs) const;

    bool contains(const NASequence &rhs);

    bool containsSNV(BulkDatum *snv);

    const string &getChr() const;

    void setChr(const string &chr);

    size_t getStartPos() const;

    void setStartPos(size_t startPos);

    size_t getEndPos() const;

    void setEndPos(size_t endPos);

    size_t getNumChr() const;

    void setNumChr(size_t numChr);

    string toString() {
        return chr + "[" + to_string(start_pos) + "-" + to_string(end_pos) + "]";
    }
};


class Gene : public NASequence {

    string ensembl_id;
    string name;

    double per_copy_expr; // mu_g
    double nb_inv_dispersion; // r_g
    double gene_copy_prob; // delta_g

public:

    Gene(string ensemblId, string chr, size_t startPos, size_t endPos, string name);

    Gene(string ensembl_id, const string& chr, size_t start_pos, size_t end_pos, double per_copy_expr,
         double nb_inv_dispersion, string name);

    const string &getEnsemblId() const;

    void setEnsemblId(const string &ensemblId);

    double getPerCopyExpr() const;

    void setPerCopyExpr(double perCopyExpr);

    const string &getName() const;

    void setName(const string &name);

    double getNbInvDispersion() const;

    void setNbInvDispersion(double nbInvDispersion);

    double getGeneCopyProb() const;

    void setGeneCopyProb(double geneCopyProb);

    static vector<Gene *> readGeneCodeFromFile(const string &path);
};

class Bin : public NASequence {
    vector<Gene *> genes;
    vector<size_t> geneIdxs;

public:
    Bin(const string& chr, size_t start_pos, size_t end_pos);
    Bin(size_t numChr, size_t start_pos, size_t end_pos);

    static vector<Bin> generateBinsFromGenes(const vector<Gene *> &gene_set, size_t bin_size);
    void insertGene(Gene *gene, size_t geneIdx);

    const vector<Gene *> &getGenes() const;

    void setGenes(const vector<Gene *> &genes);

    const vector<size_t> &getGeneIdxs() const;
};

bool comparePtrToNASeq(NASequence *lhs, NASequence *rhs);

#endif //PHYLEX_GENE_HPP
