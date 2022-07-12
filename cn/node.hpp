//
//  node.hpp
//  phylex_cn
//
//  Created by Seong-Hwan Jun on 2021-07-19.
//

#ifndef node_hpp
#define node_hpp

#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

class Node
{
    Node *parent_;
    vector<Node *> children;
    string name_;
public:
    Node(Node *parent, string name) : parent_(parent), name_(name) {
        parent_->AddChild(this);
    }
    inline void AddChild(Node *child_node) {
        children.push_back(child_node);
    }
};

#endif /* node_hpp */
