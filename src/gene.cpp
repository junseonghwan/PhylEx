//
// Created by zemp on 17/04/21.
//

#include "gene.hpp"

#include <utility>

Gene::Gene(string ensembl_id, string chr, size_t start_pos, size_t end_pos) : chr(std::move(chr)), start_pos(start_pos),
                                                                              end_pos(end_pos),
                                                                              ensembl_id(std::move(ensembl_id)) {}

Gene::Gene(string ensembl_id, string chr, size_t start_pos, size_t end_pos, double per_copy_expr) : chr(std::move(chr)),
                                                                                                    start_pos(
                                                                                                            start_pos),
                                                                                                    end_pos(end_pos),
                                                                                                    ensembl_id(
                                                                                                            std::move(
                                                                                                                    ensembl_id)),
                                                                                                    per_copy_expr(
                                                                                                            per_copy_expr) {}

const string &Gene::getChr() const {
    return chr;
}

void Gene::setChr(const string &chr) {
    Gene::chr = chr;
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

