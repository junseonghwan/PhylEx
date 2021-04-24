//
// Created by zemp on 17/04/21.
//

#include "gene.hpp"
#include <boost/algorithm/string.hpp>

Gene::Gene(string ensembl_id, string chr, size_t start_pos, size_t end_pos, double per_copy_expr = -1.,
           double nb_inv_dispersion = -1., string name = "") : chr(std::move(chr)),
                                                               start_pos(start_pos),
                                                               end_pos(end_pos),
                                                               ensembl_id(std::move(ensembl_id)),
                                                               per_copy_expr(per_copy_expr),
                                                               nb_inv_dispersion(nb_inv_dispersion),
                                                               name(std::move(name)) {}

const string &Gene::getChr() const {
    return chr;
}

void Gene::setChr(const string &chromosome) {
    chr = chromosome;
}

size_t Gene::getStartPos() const {
    return start_pos;
}

void Gene::setStartPos(size_t startPos) {
    start_pos = startPos;
}

size_t Gene::getEndPos() const {
    return end_pos;
}

void Gene::setEndPos(size_t endPos) {
    end_pos = endPos;
}

const string &Gene::getEnsemblId() const {
    return ensembl_id;
}

void Gene::setEnsemblId(const string &ensemblId) {
    ensembl_id = ensemblId;
}

double Gene::getPerCopyExpr() const {
    return per_copy_expr;
}

void Gene::setPerCopyExpr(double perCopyExpr) {
    per_copy_expr = perCopyExpr;
}

const string &Gene::getName() const {
    return name;
}

void Gene::setName(const string &name) {
    Gene::name = name;
}

double Gene::getNbInvDispersion() const {
    return nb_inv_dispersion;
}

void Gene::setNbInvDispersion(double nbInvDispersion) {
    nb_inv_dispersion = nbInvDispersion;
}

vector<Gene *> Gene::readGeneCodeFromFile(const string &path) {
    ifstream file(path);

    if (!file.is_open()) {
        cerr << "Could not open the file: " << path << endl;
        exit(-1);
    }

    string line;
    vector<Gene *> genecode;
    while(getline(file, line)) {
        vector<string> ids;
        vector<string> res;

        boost::split(res, line, boost::is_any_of("\t"));
        boost::split(ids, res[0], boost::is_any_of("|"));

        genecode.push_back(new Gene(ids[1], res[1], stoul(res[2]), stoul(res[3]),
                                    ids[0]));
    }

    file.close();

    return genecode;
}

Gene::Gene(string ensemblId, string chr, size_t startPos, size_t endPos, string name) : chr(std::move(chr)),
                                                                                        start_pos(
                                                                                                startPos),
                                                                                        end_pos(endPos),
                                                                                        ensembl_id(std::move(
                                                                                                ensemblId)),
                                                                                        name(std::move(name)) {}

