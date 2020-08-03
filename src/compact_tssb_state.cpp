
#include "compact_tssb_state.hpp"

CompactTSSBState::CompactTSSBState(TSSBState &tssb_state)
{
    newick = write_newick(tssb_state.get_root());
    fill_node_to_param(tssb_state.get_root(), node2param);

    // construct ancestral matrix
    ancestral_matrix = TSSBState::get_ancestral_matrix(tssb_state);
    
    // construct datum2param
    CloneTreeNode *node;
    for (size_t i = 0; i < tssb_state.get_data().size(); i++) {
        node = tssb_state.get_node(tssb_state.get_data().at(i));
        datum2param.push_back(new CloneTreeNodeParam(node->NodeParameter()));
        datum2node.push_back(node->GetName());
    }
    
    // get non-empty nodes
    //    vector<CloneTreeNode *> nodes;
    //    TSSBState::get_all_nodes(true, tssb_state.get_root(), nodes);
    //
    //    for (size_t i = 0; i < nodes.size(); i++) {
    //        unordered_set<S*> data(nodes[i]->get_data()); // copy the elements
    //        clustering.push_back(data);
    //    }
    
    CloneTreeNode::GetClusterLabels(tssb_state.get_root(), tssb_state.get_data(), cluster_labels);
}

CompactTSSBState::~CompactTSSBState()
{
    //    if (ancestral_matrix != 0)
    //        gsl_matrix_free(ancestral_matrix);
    for (size_t i = 0; i < datum2param.size(); i++) {
        // delete the parameter since it is a copy
        // do not delete S*, since it is used globally in the program
        delete datum2param[i];
    }
}

