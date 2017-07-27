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

#include "score_data.h"

void score_data::computeOverlapScore(scluster * curr_node, set < int > & curr_leaves, bool mode_interogate) {
	curr_leaves.clear();
	if (!curr_node->leaf()) {
		set < int > child1_leaves, child2_leaves;
		computeOverlapScore(curr_node->child1, child1_leaves, mode_interogate);
		computeOverlapScore(curr_node->child2, child2_leaves, mode_interogate);
		curr_leaves.insert(child1_leaves.begin(), child1_leaves.end());
		curr_leaves.insert(child2_leaves.begin(), child2_leaves.end());
		string str_uid = "";
		for (set < int > :: iterator it_str = curr_leaves.begin() ; it_str != curr_leaves.end() ; ++it_str) str_uid += "_" + stb.str(*it_str);
		if (mode_interogate) {
			unordered_map < string, scluster* > :: iterator it_map = map_nodes.find(str_uid);
			if (it_map != map_nodes.end()) it_map->second->bootstrap += 1.0 / tree_samples.size();
		} else map_nodes.insert(pair < string, scluster* > (str_uid, curr_node));
	} else curr_leaves.insert(curr_node->index);
}

void score_data::computeOverlapScore() {
	vrb.title("Compute overlap scores in the Nominal tree");
	set < int > curr_nominal_leaves;
	computeOverlapScore(tree_nominal, curr_nominal_leaves, false);
	for (int b = 0 ; b < tree_samples.size() ; b ++) {
		set < int > curr_sampled_leaves;
		computeOverlapScore(tree_samples[b], curr_sampled_leaves, true);
		vrb.bullet("Processing sampled tree = " + stb.str(b));
	}
}
