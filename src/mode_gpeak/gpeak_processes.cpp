#include "gpeak_data.h"

void gpeak_data::process(string file) {
	vector < int > curr_leaves;
	vec_module_ids = vector < string > (vec_nodes.size(), "NA");
	vrb.title("Determine the module IDs for each peak");
	process(tree, curr_leaves);
	vrb.bullet("peak determination finished");

	vrb.title("Write Outcome in [" + file + "]");
	int n_leaf_ymarked = 0, n_leaf_nmarked = 0;
	output_file fd(file);
	for (int n = 0 ; n < vec_module_ids.size() ; n++) {
		if (vec_nodes[n]->leaf() ) {
			if (vec_module_ids[n] != "NA") {
				fd <<vec_nodes[n]->sid << " " <<vec_module_ids[n] << endl;
				n_leaf_ymarked++;
			} else n_leaf_nmarked++;
		}
	}
	fd.close();

}

void gpeak_data::process(gcluster * curr_node, vector < int > & curr_leaves) {
	if (! curr_node->leaf()) {
		vector < int > child1_leaves, child2_leaves;
		process(curr_node->child1, child1_leaves);
		process(curr_node->child2, child2_leaves);
		curr_leaves.insert(curr_leaves.end(), child1_leaves.begin(), child1_leaves.end());
		curr_leaves.insert(curr_leaves.end(), child2_leaves.begin(), child2_leaves.end());
		sort(curr_leaves.begin(), curr_leaves.end());
		set < string > :: iterator it = modules.find(curr_node->sid);
		if (it != modules.end()) for (int n = 0 ; n < curr_leaves.size() ; n++) vec_module_ids[curr_leaves[n]] = curr_node->sid;
	} else curr_leaves.push_back(curr_node->index);
}
