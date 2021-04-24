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

using namespace std;

class Gene {
    string chr;
    size_t start_pos;
    size_t end_pos;
    string ensembl_id;
    string name;

    double per_copy_expr; // mu_g
    double nb_inv_dispersion; // r_g

public:

    Gene(string ensemblId, string chr, size_t startPos, size_t endPos, string name);

    Gene(string ensembl_id, string chr, size_t start_pos, size_t end_pos, double per_copy_expr,
         double nb_inv_dispersion, string name);

    const string &getChr() const;

    void setChr(const string &chromosome);

    size_t getStartPos() const;

    void setStartPos(size_t startPos);

    size_t getEndPos() const;

    void setEndPos(size_t endPos);

    const string &getEnsemblId() const;

    void setEnsemblId(const string &ensemblId);

    double getPerCopyExpr() const;

    void setPerCopyExpr(double perCopyExpr);

    const string &getName() const;

    void setName(const string &name);

    double getNbInvDispersion() const;

    void setNbInvDispersion(double nbInvDispersion);

    static vector<Gene *> readGeneCodeFromFile(const string &path);
};


#endif //PHYLEX_GENE_HPP
