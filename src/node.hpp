//
//  Node.hpp
//  Abstract class representing a node.
//  It must be overriden and provide a custom implementation of the likelihood
//
//  Created by Seong-Hwan Jun on 2018-05-28.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef Node_hpp
#define Node_hpp

#include <iostream>
#include <limits>
#include <map>
#include <memory.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <gsl/gsl_randist.h>

#include "model_params.hpp"
#include "sampling_utils.hpp"
#include "tssb_state.hpp"

using namespace std;

template <class S, class C, class P> class TSSBState;

template <class S, class P>
class Node
{
    // change the last part of the node's name to j -- used extensively by reorder_sticks
    void edit_name(size_t j);

protected:
    P param;
    vector<size_t> name;
    double nu = 0.0;
    Node<S,P> *parent_node = 0;
    unordered_set<S *> data;

    // map: child branch idx -> pair<psi_stick, Node *>
    unordered_map<size_t, pair<double, Node<S,P> *> > idx2child;
    bool operator==(const Node &other) const
    {
        if (this->parent_node != other.parent)
            return false;
        if (name.size() != other.name.size())
            return false;
        size_t n = name.size();
        return (name[n-1] == other.name[n-1]);
    }

public:
    Node(size_t child_idx, Node<S,P> *parent);
    virtual ~Node();

    // getters
    string get_name() const;
    inline const vector<size_t> &get_name_vec() const { return name; };
    inline size_t get_depth() const { return name.size() - 1; }
    size_t get_num_data() const;
    double get_nu_stick() const;
    Node<S,P> *get_parent_node() const;
    const unordered_set<S*> &get_data() const;
    const pair<double, Node<S,P> *> &get_child(size_t child_idx) const;
    size_t get_num_children() const;
    bool is_root() const;
    bool is_leaf() const;
    unordered_map<size_t, pair<double, Node<S,P> *> > &get_idx2child();

    // setters
    void set_nu_stick(double nu);
    void set_psi_stick(size_t child_idx, double psi);

    // add/remove datum
    void add_datum(S *datum);
    void remove_datum(S *datum);
    bool contains_datum(S *datum) const;

    // operations on sticks
    void cull(unordered_set<size_t> &cull_list);
    void reorder_sticks(const gsl_rng *random, const ModelParams &params);
    void reset_children_names();

    // identify the branch that contains u, return the corresponding child Node
    Node<S,P> *locate_child(const gsl_rng *random, double &u, const ModelParams &hyper_params);
    void InitializeChild(const gsl_rng *random,
                         const ModelParams &params);
    
    // Virtual functions to be overridden
    // get parameters
    virtual const P &get_node_parameter() const = 0;
    // sample node parameters
    virtual void sample_node_parameters(const gsl_rng *random, const ModelParams &params, Node<S,P> *parent) = 0;
    virtual Node<S,P> *spawn_child(double psi) = 0;
    // print the contents of this Node
    virtual string print() = 0;

    static string form_node_string(string curr_node_str, size_t branch);
    static string get_parent_string(string curr_node_str);

    static void breadth_first_traversal(Node<S,P> *root, vector<Node<S,P> *> &ret, bool non_empty = false);
    static Node<S,P> *find_node(const gsl_rng *random, double u, Node<S,P> *root, const ModelParams &params);
    // return all data from ancestors including node itself
    static void get_dataset(Node<S,P> *node, unordered_set<const S*> &dataset);
    static void get_cluster_labels(Node<S,P> *root,
                                   const vector<S *> &data,
                                   vector<unsigned int> &cluster_labels);
    static void construct_datum2node(vector<Node<S,P> *> &all_nodes,
                                                unordered_map<S*, Node<S,P> *> &datum2node);
    static gsl_matrix *GetAncestralMatrix(Node<S,P> *root,
                                            const vector<S*> &bulk_data,
                                            const unordered_map<S*, Node<S,P> *> &datum2node);

    static bool less(Node<S,P> *l, Node<S,P> *r);
    friend bool operator<(const Node<S,P>& lhs, const Node<S,P>& rhs) {
        vector<string> lhs_arr, rhs_arr;
        boost::split(lhs_arr, lhs.name, boost::is_any_of("_"));
        boost::split(rhs_arr, rhs.name, boost::is_any_of("_"));
        if (lhs_arr.size() < rhs_arr.size())
            return true;
        else if (lhs_arr.size() > rhs_arr.size())
            return false;
        else {
            // same length, compare one by one
            for (size_t i = 0; i < lhs_arr.size(); i++) {
                if (stoi(lhs_arr[i]) < stoi(rhs_arr[i])) {
                    return true;
                } else if (stoi(lhs_arr[i]) > stoi(rhs_arr[i])) {
                    return false;
                }
            }
            return false;
        }
        //return tie(lhs.name) < tie(rhs.name);
    }
    friend bool operator> (const Node<S,P>& lhs, const Node<S,P>& rhs){ return rhs < lhs; }
    friend bool operator<=(const Node<S,P>& lhs, const Node<S,P>& rhs){ return !(lhs > rhs); }
    friend bool operator>=(const Node<S,P>& lhs, const Node<S,P>& rhs){ return !(lhs < rhs); }
};

template <class S, class P>
string Node<S,P>::get_parent_string(string curr_node_str)
{
    // split the string using "_"
    vector<string> fields;
    boost::split(fields, curr_node_str, boost::is_any_of("_"));
    string ret = "";
    size_t num_iter = fields.size() - 1;
    for (size_t i = 0; i < num_iter; i++)
    {
        ret += fields[i];
        if (i < num_iter - 1) {
            ret += "_";
        }
    }
    return ret;
}

template <class S, class P>
string Node<S,P>::form_node_string(string curr_node_str, size_t branch)
{
    if (curr_node_str == "")
        return to_string(branch);
    
    return curr_node_str + "_" + to_string(branch);
}

template <class S, class P>
Node<S,P>::Node(size_t child_idx, Node<S,P> *parent) :
parent_node(parent)
{
    if (parent != 0) {
        this->name = parent->name; // copy the name vector
    }
    this->name.push_back(child_idx);
}

template <class S, class P>
Node<S,P>::~Node()
{
}

template <class S, class P>
Node<S,P> *Node<S,P>::find_node(const gsl_rng *random, double u, Node<S,P> *root, const ModelParams &params)
{
    Node<S,P> *node = root;
    double nu = 0.0;
    size_t depth = 0;
    while (true) {
        depth = node->get_name_vec().size();
        if (depth >= params.get_max_depth()) {
            break;
        }
        nu = node->get_nu_stick();

        if (u < nu) {
            // found the node!!
            break;
        }
        
        // shrink u relative to the remaining stick (1 - nu)
        u = (u - nu) / (1 - nu);
        
        // find sub branch by enumerating over the branching sticks
        node = node->locate_child(random, u, params);
    }
    return node; // TODO: consider returning the path for the purposes of better estimation of the parameters
}

template <class S, class P>
string Node<S,P>::get_name() const
{
    // concatenate vector into string
    string str = "";
    for (size_t i = 0; i < name.size(); i++) {
        str = this->form_node_string(str, name[i]);
    }
    return str;
}

template <class S, class P>
void Node<S,P>::edit_name(size_t j)
{
    size_t n = this->name.size();
    for (size_t i = 0; i < n - 1; i++) {
        name[i] = parent_node->name[i];
    }
    name[n-1] = j;
}

// getters
template <class S, class P>
Node<S,P> *Node<S,P>::get_parent_node() const
{
    return parent_node;
}

template <class S, class P>
const unordered_set<S*> &Node<S,P>::get_data() const
{
    return data;
}

template <class S, class P>
void Node<S,P>::get_dataset(Node<S,P> *node,
                            unordered_set<const S*> &dataset)
{
    // trace up to the root node to get all SNVs
    while (node != 0) {
        dataset.insert(node->get_data().begin(), node->get_data().end());
        node = node->get_parent_node();
    }
}

template <class S, class P>
size_t Node<S,P>::get_num_data() const
{
    return data.size();
}

template <class S, class P>
unordered_map<size_t, pair<double, Node<S,P> *> > &Node<S,P>::get_idx2child()
{
    return idx2child;
}

template <class S, class P>
double Node<S,P>::get_nu_stick() const
{
    return this->nu;
}

template <class S, class P>
size_t Node<S,P>::get_num_children() const
{
    return idx2child.size();
}

template <class S, class P>
bool Node<S,P>::is_root() const
{
    return parent_node == 0;
}

template <class S, class P>
bool Node<S,P>::is_leaf() const
{
    return (get_num_children() == 0);
}

template <class S, class P>
const pair<double, Node<S,P> *> &Node<S,P>::get_child(size_t child_idx) const
{
    if (child_idx < idx2child.size()) {
        return idx2child.at(child_idx);
    }
    cerr << "Error: child index out of bounds." << endl;
    exit(-1);
}

template <class S, class P>
void Node<S,P>::set_nu_stick(double nu)
{
    this->nu = nu;
}

template <class S, class P>
void Node<S,P>::set_psi_stick(size_t child_idx, double psi)
{
    pair<double, Node*> &child = idx2child.at(child_idx);
    child.first = psi;
}

template <class S, class P>
void Node<S,P>::add_datum(S *datum)
{
    data.insert(datum);
}

template <class S, class P>
void Node<S,P>::remove_datum(S *datum)
{
    data.erase(datum);
}

template <class S, class P>
bool Node<S,P>::contains_datum(S *datum) const
{
    return (data.count(datum) > 0);
}

template <class S, class P>
void Node<S,P>::cull(unordered_set<size_t> &cull_list)
{
    // need to update the cullular prevalence as the nodes are culled (unfortunately, this needs to happen outside by the one that calls cull operation)
    unordered_map<size_t, pair<double, Node<S,P> *> > new_idx2child;
    for (size_t j = 0; j < idx2child.size(); j++)
    {
        Node<S,P> *node = idx2child[j].second;
        if (cull_list.count(j) == 0) {
            new_idx2child[j] = idx2child[j];
        } else {
            //cout << get_name() << " cull child " << j << endl;
            idx2child[j].second = 0;
            delete node;
        }
    }
    idx2child = new_idx2child;
}

template <class S, class P>
void Node<S,P>::reset_children_names()
{
    // clear the children info
    unordered_map<size_t, pair<double, Node<S,P> *> > new_map;

    Node<S,P> *child = 0;
    size_t n_children = idx2child.size();
    size_t idx = 0;
    for (size_t i = 0; i < n_children; i++) {
        child = idx2child[i].second;
        if (child == 0) {
            //idx2child.erase(i); // this is not necessarily since we will overwrite idx2child
        } else {
            child->edit_name(idx);
            new_map[idx] = idx2child[i];
            idx++;
        }
    }
    idx2child = new_map;
}


template <class S, class P>
void Node<S,P>::reorder_sticks(const gsl_rng *random, const ModelParams &params)
{
    if (idx2child.size() == 0)
        return;
    
    //vector<double> unnorm_w(idx2child.size());
    vector<double> weights(idx2child.size());
    vector<double> intervals(idx2child.size());
    unordered_set<int> represented;
    double cum_prod = 1.0;
    for (size_t i = 0; i < idx2child.size(); i++)
    {
        weights[i] = idx2child[i].first * cum_prod;
        cum_prod *= (1 - idx2child[i].first);
        represented.insert(i);
        if (i == 0)
            intervals[i] = weights[i];
        else
            intervals[i] = intervals[i-1] + weights[i];
    }

    unordered_map<size_t, pair<double, Node<S,P> *> > new_map;

    // throw uniform darts until all represented sticks are sampled
    double represented_stick_length = intervals[idx2child.size()-1];
    size_t n_children = intervals.size();
    while (represented.size() > 0) {
        double u = gsl_ran_flat(random, 0, 1);
        while (u > represented_stick_length) {
            // need to represent new children: i.e., draw psi sticks
            double psi_j = bounded_beta(random, 1, params.get_gamma());
            Node<S,P> *child = this->spawn_child(psi_j);
            double nu_stick = bounded_beta(random, 1.0, params.alpha(child->get_name_vec()));
            child->set_nu_stick(nu_stick);
            child->sample_node_parameters(random, params, this);
            
            double ww = psi_j * cum_prod;
            weights.push_back(ww);
            intervals.push_back(intervals[n_children - 1] + ww);
            represented.insert(n_children);
            represented_stick_length += ww;
            n_children++;
        }

        for (size_t i = 0; i < n_children; i++) {
            if (u < intervals[i]) {
                // found the child
                pair<double, Node<S,P> *> &ret = idx2child[i];
                ret.second->edit_name(i);
                new_map[i] = ret;
                represented.erase(i);
                break;
            }
        }
        
        // renormalize?
    }
    
    idx2child = new_map;
}

//template <class S, class P>
//void Node<S,P>::reorder_sticks(const gsl_rng *random)
//{
//    if (idx2child.size() <= 1)
//        return;
//
//    vector<double> unnorm_w(idx2child.size());
//    double norm = 0.0;
//    double cum_prod = 1.0;
//    for (size_t i = 0; i < idx2child.size(); i++)
//    {
//        double ww = idx2child[i].first;
//        cum_prod *= (1 - ww);
//        unnorm_w[i] = (1 - cum_prod) - norm;
//        if (unnorm_w[i] > 0) {
//            unnorm_w[i] = std::numeric_limits< double >::min();
//        }
//        assert(unnorm_w[i] > 0);
//        norm += unnorm_w[i];
//    }
//    size_t N = unnorm_w.size();
//    //this->children_node_str.clear(); // remove all names
//    unordered_map<size_t, pair<double, Node<S,P> *> > new_map;
//    int idx;
//    for (size_t i = 0; i < N; i++)
//    {
//        if (norm == 0.0) {
//            // the remaining sticks have 0 measure -- just keep the current ordering
//            idx = i;
//        } else {
//            idx = multinomial(random, unnorm_w, norm);
//        }
//        pair<double, Node<S,P> *> &ret = idx2child[idx];
//        ret.second->edit_name(i);
//        //this->children_node_str.insert(ret.second->get_name_vec());
//        new_map[i] = ret;
//        norm -= unnorm_w[idx];
//        unnorm_w[idx] = 0.0;
//    }
//    idx2child = new_map;
//}

template <class S, class P>
void Node<S,P>::InitializeChild(const gsl_rng *random,
                                const ModelParams &params)
{
    double psi_j = bounded_beta(random, 1, params.get_gamma());
    Node<S,P> *child = this->spawn_child(psi_j);
    double nu_stick = bounded_beta(random, 1.0, params.alpha(child->get_name_vec()));
    child->set_nu_stick(nu_stick);
    child->sample_node_parameters(random, params, this);
}

template <class S, class P>
Node<S,P> *Node<S,P>::locate_child(const gsl_rng *random, double &u, const ModelParams &params)
{
    unsigned int j = 0;
    double cum_prod = 1.0;
    double begin = 0;
    while (true) {
        if (get_num_children() <= j) { // draw new psi-stick and create a new child
            InitializeChild(random, params);
            /*
            double psi_j = bounded_beta(random, 1, params.get_gamma());
            Node<S,P> *child = this->spawn_child(psi_j);
            double nu_stick = bounded_beta(random, 1.0, params.alpha(child->get_name_vec()));
            child->set_nu_stick(nu_stick);
            child->sample_node_parameters(random, params, this);
             */
        }
        pair<double, Node*> pair = idx2child[j];
        double interval_length = cum_prod * pair.first;
        if (u >= begin && u < (begin + interval_length)) {
            // found the corresponding interval
            u = (u - begin) / interval_length;
            break;
        }
        begin += interval_length;
        cum_prod *= (1 - pair.first);
        j++;
    }
    return idx2child[j].second;
}

template <class S, class P>
bool Node<S,P>::less(Node<S,P> *lhs, Node<S,P> *rhs)
{
    size_t lhs_size = lhs->name.size();
    size_t rhs_size = rhs->name.size();
    if (rhs_size == 0 && lhs_size == 0) {
        return false;
    } else if (lhs_size == 0) {
        return true;
    } else if (rhs_size == 0) {
        return false;
    }
    
    size_t len = min(lhs_size, rhs_size);
    for (size_t i = 0; i < len; i++) {
        if (lhs->name[i] < rhs->name[i]) {
            return true;
        } else if (lhs->name[i] > rhs->name[i]) {
            return false;
        }
    }
    
    // if we are here, the two strings are same up to len
    // then, we determine the order one is shorter than the other
    if (lhs_size < rhs_size) {
        return true;
    } else {
        return false;
    }
}

template <class S, class P>
void Node<S,P>::breadth_first_traversal(Node<S,P> *root_node, vector<Node<S,P> *> &ret, bool non_empty)
{
    // get nodes in the subtree rooted at root
    if (ret.size() > 0) {
        cout << "Warning: get_all_nodes() is going to clear vector ret." << endl;
        ret.clear();
    }
    if (root_node == 0) {
        cout << "Warning: root node is null." << endl;
        return;
    }

    set<Node<S,P> *> nodes;
    Node<S,P> *node = 0;

    queue<Node<S,P> *> q;
    q.push(root_node);
    while (!q.empty())
    {
        node = q.front();
        q.pop();
        if (!nodes.count(node)) {
            if (!non_empty) {
                ret.push_back(node);
                nodes.insert(node);
            } else {
                if (node->get_num_data() > 0) {
                    ret.push_back(node);
                    nodes.insert(node);
                }
            }
        }
        unordered_map<size_t, pair<double, Node<S,P> *> > &children = node->get_idx2child();
        for (size_t i = 0; i < children.size(); i++)
        {
            Node<S,P> *child = children[i].second;
            q.push(child);
        }
    }
}

template <class S, class P>
void Node<S,P>::get_cluster_labels(Node<S,P> *root,
                                   const vector<S *> &data,
                                   vector<unsigned int> &cluster_labels)
{
    vector<Node<S,P> *> all_nodes;
    Node<S,P>::breadth_first_traversal(root, all_nodes, false);
    unordered_map<S *, Node<S, P> *> datum2node;
    construct_datum2node(all_nodes, datum2node);

    // determine the classes and clusters
    // assign class label from 0, ..., N = num nodes
    unordered_map<Node<S,P> *, size_t> node2class;
    size_t num_classes = all_nodes.size();
    for (size_t c = 0; c < num_classes; c++) {
        node2class[all_nodes[c]] = c;
    }
    for (size_t i = 0; i < data.size(); i++) {
        Node<S,P> *v = datum2node[data[i]];
        size_t c = node2class[v];
        cluster_labels.push_back(c);
    }
}

template <class S, class P>
void Node<S,P>::construct_datum2node(vector<Node<S,P> *> &all_nodes,
                                     unordered_map<S*, Node<S,P> *> &datum2node)
{
    for (auto it = all_nodes.begin(); it != all_nodes.end(); ++it) {
        for (S *datum : (*it)->get_data())
        {
            datum2node[datum] = (*it);
        }
    }
}

template <class S, class P>
gsl_matrix *Node<S,P>::GetAncestralMatrix(Node<S,P> *root,
                                          const vector<S*> &data,
                                          const unordered_map<S*, Node<S,P> *> &datum2node)
{
    size_t N = data.size();
    vector<unordered_set<S*> > ancestors(N);
    for (size_t i = 0; i < N; i++) {
        S *datum = data.at(i);
        Node<S,P> *v = datum2node.at(datum);
        // trace up to the root and set the row of A
        while (v != root) {
            v = v->get_parent_node();
            unordered_set<S *> d = v->get_data();
            ancestors[i].insert(d.begin(), d.end());
        }
    }
    
    // A_ij = 1 if i occurs before j
    gsl_matrix *A = gsl_matrix_alloc(N, N);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            if (ancestors[j].count(data.at(i))) {
                gsl_matrix_set(A, i, j, 1);
            } else {
                gsl_matrix_set(A, i, j, 0);
            }
        }
    }
    
    return A;
}

#endif /* Node_hpp */
