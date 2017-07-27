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

#include "call_data.h"

void call_data::tagSignificant(ccluster * curr_node, int curr_mod, float athreshold, float ethreshold, vector < int > & curr_leaves) {
	curr_leaves.clear();
	if (!curr_node->leaf()) {
		bool isModule = false;
		vector < int > child1_leaves, child2_leaves;
		if (curr_node->acorrelation >= athreshold * tree->acorrelation && curr_node->ecorrelation >= ethreshold * tree->ecorrelation) {
			tagSignificant(curr_node->child1, curr_mod+1, athreshold+1, ethreshold+1, child1_leaves);
			tagSignificant(curr_node->child2, curr_mod+1, athreshold+1, ethreshold+1, child2_leaves);
			isModule = true;
			curr_node->module = curr_mod + 1;
		} else {
			tagSignificant(curr_node->child1, curr_mod, athreshold, ethreshold, child1_leaves);
			tagSignificant(curr_node->child2, curr_mod, athreshold, ethreshold, child2_leaves);
			curr_node->module = 0;
		}
		curr_leaves.insert(curr_leaves.end(), child1_leaves.begin(), child1_leaves.end());
		curr_leaves.insert(curr_leaves.end(), child2_leaves.begin(), child2_leaves.end());
		sort(curr_leaves.begin(), curr_leaves.end());

		if (isModule) {
			for (int n = 0 ; n < curr_leaves.size() ; n ++) vec_nodes[curr_leaves[n]]->module = -1;
		}

	} else curr_leaves.push_back(curr_node->index);
}
