/*Copyright (C) 2015 Olivier Delaneau, Halit Ongen, Emmanouil T. Dermitzakis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include "topo_data.h"

void topo_data::characterizeTopology(tcluster * curr_node, vector < int > & curr_leaves) {
	curr_leaves.clear();
	if (curr_node->leaf()) {
		curr_leaves.push_back(curr_node->index);
	} else {
		vector < int > child1_leaves, child2_leaves;
		characterizeTopology(curr_node->child1, child1_leaves);
		characterizeTopology(curr_node->child2, child2_leaves);

		//Merging children leaves
		curr_leaves.insert(curr_leaves.end(), child1_leaves.begin(), child1_leaves.end());
		curr_leaves.insert(curr_leaves.end(), child2_leaves.begin(), child2_leaves.end());
		sort(curr_leaves.begin(), curr_leaves.end());

		//Number of siblings
		curr_node->child1->n_siblings = child2_leaves.size();
		curr_node->child2->n_siblings = child1_leaves.size();

		//Branching ratio
		curr_node->branching_ratio = child1_leaves.size() * 1.0 / (child1_leaves.size() + child2_leaves.size());

		//Boundary indexes
		curr_node->left_most_index = curr_leaves[0];
		curr_node->right_most_index = curr_leaves.back();

		//Dispersion
		for (int n = 1 ; n < curr_leaves.size() ; n ++) {
			float curr_dispersion = (curr_leaves[n] - curr_leaves[n-1] + 1) * 1.0 / (curr_leaves.back() - curr_node->left_most_index + 1);
			if (curr_node->dispersion < curr_dispersion) curr_node->dispersion = curr_dispersion;
		}

		//Completeness
		curr_node->completeness = curr_node->count * 1.0 / (curr_leaves.back() - curr_node->left_most_index + 1);

		//Rep
		bool hit = false;
		curr_node->n_rep = 0;
		if (curr_leaves[0] == child1_leaves[0] && curr_leaves.back() == child1_leaves.back()) {
			curr_node->n_rep += curr_node->child1->n_rep;
			curr_node->child1->n_rep = 0;
			hit = true;
		}
		if (curr_leaves[0] == child2_leaves[0] && curr_leaves.back() == child2_leaves.back()) {
			curr_node->n_rep += curr_node->child2->n_rep;
			curr_node->child2->n_rep = 0;
			hit = true;
		}
		if (hit) curr_node->n_rep += 1;
	}
}

void topo_data::characterizeTopology() {
	vector < int > root_leaves;
	vrb.title("Characterize topology of the tree");
	characterizeTopology(tree, root_leaves);
}

