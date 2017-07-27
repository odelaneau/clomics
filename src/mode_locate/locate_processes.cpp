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

#include "locate_data.h"

void locate_data::propagateRegulatoryElements(lcluster * curr_node) {
	if (!curr_node->leaf()) {
		propagateRegulatoryElements(curr_node->child1);
		propagateRegulatoryElements(curr_node->child2);

		//Copy regulatory elements from kids
		curr_node->regulatory_elements.clear();
		curr_node->regulatory_elements.insert(curr_node->regulatory_elements.end(), curr_node->child1->regulatory_elements.begin(), curr_node->child1->regulatory_elements.end());
		curr_node->regulatory_elements.insert(curr_node->regulatory_elements.end(), curr_node->child2->regulatory_elements.begin(), curr_node->child2->regulatory_elements.end());
		sort(curr_node->regulatory_elements.begin(), curr_node->regulatory_elements.end());

		//Merge regulatory when necessary
		int target_e1 = 0, target_e2 = 0;
		while (target_e1 >= 0 && target_e2 >= 0) {
			target_e1 = -1;
			target_e2 = -1;
			for (int e1 = 0 ; e1 < curr_node->regulatory_elements.size() && target_e1 < 0 && target_e2 < 0 ; e1 ++) {
				for (int e2 = e1 + 1 ; e2 < curr_node->regulatory_elements.size() && target_e1 < 0 && target_e2 < 0 ; e2 ++) {
					if (curr_node->regulatory_elements[e1].overlap(curr_node->regulatory_elements[e2])) {
						target_e1 = e1;
						target_e2 = e2;
					}
				}
			}
			if (target_e1 >= 0 && target_e2 >= 0) {
				curr_node->regulatory_elements[target_e1] = regulatory_element(curr_node->regulatory_elements[target_e1], curr_node->regulatory_elements[target_e2]);
				curr_node->regulatory_elements.erase(curr_node->regulatory_elements.begin() + target_e2);
			}
		}

		//Update distance, start, end, length;
		curr_node->start = min(curr_node->child1->start, curr_node->child2->start);
		curr_node->end = max(curr_node->child1->end, curr_node->child2->end);
		curr_node->length = curr_node->end - curr_node->start;
		curr_node->distance = 0;
		if (curr_node->child1->end < curr_node->child2->start) curr_node->distance = curr_node->child2->start - curr_node->child1->end;
		if (curr_node->child2->end < curr_node->child1->start) curr_node->distance = curr_node->child1->start - curr_node->child2->end;
	}
}

void locate_data::propagateRegulatoryElements() {
	vrb.title("Propagate genomic locations within the tree");
	propagateRegulatoryElements(tree);
}
