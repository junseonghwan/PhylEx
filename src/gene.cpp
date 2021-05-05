//
// Created by zemp on 17/04/21.
//

#include "gene.hpp"
#include <boost/algorithm/string.hpp>
#include <utility>
#include <queue>

Gene::Gene(string ensemblId, string chr, size_t startPos, size_t endPos, string name) :
        NASequence(chr, startPos, endPos),
        ensembl_id(std::move(ensemblId)),
        name(std::move(name)) {}

Gene::Gene(string ensembl_id, const string &chr, size_t start_pos, size_t end_pos, double per_copy_expr = -1.,
           double nb_inv_dispersion = -1., string name = "") : NASequence(chr, start_pos, end_pos),
                                                               ensembl_id(std::move(ensembl_id)),
                                                               per_copy_expr(per_copy_expr),
                                                               nb_inv_dispersion(nb_inv_dispersion),
                                                               name(std::move(name)) {}


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
    while (getline(file, line)) {
        vector<string> ids;
        vector<string> res;
        vector<string> ens_id; // remove the dot at the end

        boost::split(res, line, boost::is_any_of("\t"));
        boost::split(ids, res[0], boost::is_any_of("|"));
        boost::split(ens_id, ids[1], boost::is_any_of("."));
        res[1].erase(res[1].begin(), res[1].begin() + 3); // drop "chr" from chromosome id
        genecode.push_back(new Gene(ens_id[0], res[1], stoul(res[2]), stoul(res[3]),
                                    ids[0]));
    }

    file.close();

    return genecode;
}

double Gene::getGeneCopyProb() const {
    return gene_copy_prob;
}

void Gene::setGeneCopyProb(double geneCopyProb) {
    gene_copy_prob = geneCopyProb;
}

// TODO implement other smart binning methods
/**
 * Creates bins with a fixed size given a set of genes
 *
 * @param sortedGeneSet ordered gene set
 * @param bin_size size of each bin
 * @return vector of Bin
 */
vector<Bin> Bin::generateBinsFromGenes(const vector<Gene *> &sortedGeneSet, size_t bin_size) {

    // find max position for each chromosome
    vector<size_t> chrMaxPos(25, 0);
    for (auto gene: sortedGeneSet) {
        if (chrMaxPos[gene->getNumChr() - 1] < gene->getEndPos()) {
            chrMaxPos[gene->getNumChr() - 1] = gene->getEndPos();
        }
    }

    // create bins
    vector<Bin> bin_set;
    queue<Gene *> sortedGeneQueue(deque(sortedGeneSet.begin(), sortedGeneSet.end()));
    size_t geneIdx = 0;
    // chrMaxPos is ordered by chromosomes
    for (size_t chr = 0; chr < chrMaxPos.size(); ++chr) {
        int n_bins = chrMaxPos[chr] / bin_size + 1;
        for (int i = 0; i < n_bins; ++i) {
            size_t start_pos = bin_size * i + 1;
            Bin bin = Bin(chr + 1, start_pos, start_pos + bin_size);
            // only save bins in which there is at least one gene
            bool containsGene = false;

            while (!sortedGeneQueue.empty() && bin.contains(*sortedGeneQueue.front())) {
                bin.insertGene(sortedGeneQueue.front(), geneIdx++);
                sortedGeneQueue.pop();
                containsGene = true;
            }
            if (containsGene) {
                bin_set.push_back(bin);
            }
        }
    }

    assert(geneIdx == sortedGeneSet.size());
    return bin_set;
}

Bin::Bin(const string &chr, size_t start_pos, size_t end_pos) :
        NASequence(chr, start_pos, end_pos) {}

void Bin::insertGene(Gene *gene, size_t geneIdx) {
    genes.push_back(gene);
    geneIdxs.push_back(geneIdx);
}

const vector<Gene *> &Bin::getGenes() const {
    return genes;
}

void Bin::setGenes(const vector<Gene *> &genes) {
    Bin::genes = genes;
}

Bin::Bin(size_t numChr, size_t start_pos, size_t end_pos) : NASequence(numChr, start_pos, end_pos) {}

const vector<size_t> &Bin::getGeneIdxs() const {
    return geneIdxs;
}

NASequence::NASequence(const string &chr, size_t start_pos, size_t end_pos) : chr(chr),
                                                                              start_pos(start_pos),
                                                                              end_pos(end_pos) {
    // convert chr to number
    if (chr == "X") {
        numChr = 23;
    } else if (chr == "Y") {
        numChr = 24;
    } else if (chr == "M") {
        numChr = 25;
    } else {
        numChr = stoul(chr);
    }
}

bool NASequence::operator<(const NASequence &rhs) const {
    // FIXME genes can overlap, this might be an issue
    return numChr < rhs.numChr || (numChr == rhs.numChr && start_pos < rhs.start_pos);
}

const string &NASequence::getChr() const {
    return chr;
}

void NASequence::setChr(const string &chr) {
    NASequence::chr = chr;
}

size_t NASequence::getStartPos() const {
    return start_pos;
}

void NASequence::setStartPos(size_t startPos) {
    start_pos = startPos;
}

size_t NASequence::getEndPos() const {
    return end_pos;
}

void NASequence::setEndPos(size_t endPos) {
    end_pos = endPos;
}

bool NASequence::contains(const NASequence &seq) {
    return numChr == seq.numChr && start_pos <= seq.start_pos && end_pos > seq.start_pos;
}

size_t NASequence::getNumChr() const {
    return numChr;
}

void NASequence::setNumChr(size_t numChr) {
    NASequence::numChr = numChr;
}

NASequence::NASequence(size_t numChr, size_t start_pos, size_t end_pos) : numChr(numChr),
                                                                          start_pos(start_pos),
                                                                          end_pos(end_pos) {
    assert(numChr > 0 && numChr < 26);
    // convert numChr to string
    switch (numChr) {
        case 23:
            chr = "X";
            break;
        case 24:
            chr = "Y";
            break;
        case 25:
            chr = "M";
            break;
        default:
            chr = to_string(numChr);
    }
    // FIXME the number for X, Y and M. X should be 22 and 22 should not be a chromosome number
}

bool NASequence::containsSNV(BulkDatum *snv) {
    return chr == snv->GetLocus().get_chr()
        && start_pos <= snv->GetLocus().get_pos() && end_pos >= snv->GetLocus().get_pos();
}

bool comparePtrToNASeq(NASequence *lhs, NASequence *rhs) { return *lhs < *rhs; }

