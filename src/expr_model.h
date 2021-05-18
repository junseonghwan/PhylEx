//
// Created by zemp on 17/05/21.
//

#ifndef PHYLEX_EXPR_MODEL_H
#define PHYLEX_EXPR_MODEL_H

#include <iostream>

enum expr_model {
    POISSON, NEG_BINOM, ZIP, ZINB
};

inline std::ostream& operator<<(std::ostream &os, expr_model model) {
    switch (model) {
        case POISSON:
            os << "POISSON";
            break;
        case ZIP:
            os << "ZIP";
            break;
        case ZINB:
            os << "ZINB";
            break;
        case NEG_BINOM:
            os << "NEG_BINOM";
            break;
    }
    return os;
}

#include <string>
#include <set>
#include <fstream>
#include <boost/program_options.hpp>
#include "single_cell.hpp"

#endif //PHYLEX_EXPR_MODEL_H
