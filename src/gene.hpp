//
// Created by zemp on 17/04/21.
//

#ifndef PHYLEX_GENE_HPP
#define PHYLEX_GENE_HPP

#include <string>

using namespace std;

class Gene {
    string chr;
    size_t start_pos;
    size_t end_pos;
    string ensembl_id;

    double per_copy_expr{}; // mu

public:

    Gene(string ensembl_id, string chr, size_t start_pos, size_t end_pos);
    Gene(string ensembl_id, string chr, size_t start_pos, size_t end_pos, double per_copy_expr);

    const string &getChr() const;

    void setChr(const string &chr);

    size_t getStartPos() const;

    void setStartPos(size_t startPos);

    size_t getEndPos() const;

    void setEndPos(size_t endPos);

    const string &getEnsemblId() const;

    void setEnsemblId(const string &ensemblId);

    double getPerCopyExpr() const;

    void setPerCopyExpr(double perCopyExpr);
};


#endif //PHYLEX_GENE_HPP
