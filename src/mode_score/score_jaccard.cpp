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

int score_data::computeJaccardScore(scluster * curr_sampled_node, vector < bool > & curr_observed_mask, int curr_observed_size, float & best_score) {
	best_score = 0.0;
	if (!curr_sampled_node->leaf()) {
		float child1_score, child2_score;
		int match_child1 = computeJaccardScore(curr_sampled_node->child1, curr_observed_mask, curr_observed_size, child1_score);
		int match_child2 = computeJaccardScore(curr_sampled_node->child2, curr_observed_mask, curr_observed_size, child2_score);

		//Update current score
		float curr_score = (match_child1 + match_child2) * 1.0 / (curr_sampled_node->count + curr_observed_size - match_child1 - match_child2);

		if (curr_score >= child1_score && curr_score >= child2_score) best_score = curr_score;
		if (child1_score >= curr_score && child1_score >= child2_score) best_score = child1_score;
		if (child2_score >= curr_score && child2_score >= child1_score) best_score = child2_score;

		//Update number of matches
		return match_child1 + match_child2;
	} else return (curr_observed_mask[curr_sampled_node->index]?1:0);
}

void score_data::computeJaccardScore(scluster * curr_observed_node, vector < bool > & curr_observed_mask) {
	curr_observed_mask = vector < bool > (tree_nominal->count, false);
	if (!curr_observed_node->leaf()) {
		vector < bool > child1_observed_mask, child2_observed_mask;
		computeJaccardScore(curr_observed_node->child1, child1_observed_mask);
		computeJaccardScore(curr_observed_node->child2, child2_observed_mask);

		for (int n = 0 ; n < curr_observed_mask.size() ; n ++) curr_observed_mask[n] = (child1_observed_mask[n] || child2_observed_mask[n]);

		curr_observed_node->jaccard = 0.0;
		for (int b = 0 ; b < tree_samples.size() ; b++) {
			float best_score = 0.0;
			computeJaccardScore(tree_samples[b], curr_observed_mask, curr_observed_node->count, best_score);
			curr_observed_node->jaccard += best_score;
		}
		curr_observed_node->jaccard /= tree_samples.size();
		vrb.bullet("Processing node = " + stb.str(curr_observed_node->index));
	} else curr_observed_mask[curr_observed_node->index] = true;
}

void score_data::computeJaccardScore() {
	vrb.title("Compute Jaccard scores in the Nominal tree");
	vector < bool > curr_observed_mask;
	computeJaccardScore(tree_nominal, curr_observed_mask);
}
