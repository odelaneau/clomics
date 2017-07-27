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

#include "quantify_data.h"

void quantify_data::normalTransform(vector < float > & V) {
	vector < float > R;
	myranker::rank(V, R);
	double max = 0;
	for (int s = 0 ; s < sample_count ; s ++) {
		R[s] = R[s] - 0.5;
		if (R[s] > max) max = R[s];
	}
	max = max + 0.5;
	for (int s = 0 ; s < sample_count ; s ++) {
		R[s] /= max;
		V[s] = qnorm(R[s], 0.0, 1.0, 1, 0);
	}
}

void quantify_data::quantifyGroups(int method, int arg) {
	vec_quantifications = vector < vector < float > > (vec_groups.size(), vector < float > (sample_count, 0.0));
	vrb.title("Quantify activity of groups");
	for (int g = 0 ; g < vec_groups.size() ; g ++) {
		switch (method) {
		case METH_LOO: quantify_by_loo(vec_groups[g], vec_quantifications[g], arg != 0); break;
		case METH_PCA: quantify_by_pca(vec_groups[g], vec_quantifications[g], arg - 1); break;
		case METH_MEA: quantify_by_mean(vec_groups[g], vec_quantifications[g]); break;
		}
	}
	vrb.bullet("quantification finished");
}

void quantify_data::quantifyTrees(int method, int arg) {
	vector < int > curr_leaves;
	vec_quantifications = vector < vector < float > > (vec_nodes.size(), vector < float > (sample_count, 0.0));
	vrb.title("Quantify activity of internal nodes");
	switch (method) {
	case METH_LOO: vrb.bullet("Stability / " + (arg==0)?"Mean":"Variance"); break;
	case METH_PCA: vrb.bullet("Activity / PC" + stb.str(arg)); break;
	case METH_MEA: vrb.bullet("Activity / Mean"); break;
	}
	quantifyTrees(method, arg, tree, curr_leaves);
	vrb.bullet("quantification finished");
}

void quantify_data::quantifyTrees(int method, int arg, qcluster * curr_node, vector < int > & curr_leaves) {
	if (! curr_node->leaf()) {
		vector < int > child1_leaves, child2_leaves;
		quantifyTrees(method, arg, curr_node->child1, child1_leaves);
		quantifyTrees(method, arg, curr_node->child2, child2_leaves);

		curr_leaves.insert(curr_leaves.end(), child1_leaves.begin(), child1_leaves.end());
		curr_leaves.insert(curr_leaves.end(), child2_leaves.begin(), child2_leaves.end());
		sort(curr_leaves.begin(), curr_leaves.end());

		//
		vec_quantifications[curr_node->index] = vector < float > (sample_count, 0.0);
		if (modules.size() == 0 || modules.count(curr_node->sid)) {
			switch (method) {
			case METH_LOO: quantify_by_loo(curr_leaves, vec_quantifications[curr_node->index], arg != 0); break;
			case METH_PCA: quantify_by_pca(curr_leaves, vec_quantifications[curr_node->index], arg - 1); break;
			case METH_MEA: quantify_by_mean(curr_leaves, vec_quantifications[curr_node->index]); break;
			}
		}
		curr_node->start = min(curr_node->child1->start, curr_node->child2->start);
		curr_node->end = max(curr_node->child1->end, curr_node->child2->end);
		curr_node->first = curr_leaves[0];
		curr_node->last = curr_leaves.back();
	} else {
		curr_leaves.push_back(curr_node->index);
		curr_node->start = phenotype_start[curr_node->index];
		curr_node->end = phenotype_end[curr_node->index];
		vec_quantifications[curr_node->index] = phenotype_val[curr_node->index];
	}
}

void quantify_data::quantify_by_mean(vector < int > & curr_leaves, vector < float > & Q) {
	vrb.bullet("mean on " + stb.str(curr_leaves.size()));
	for (int s = 0 ; s < sample_count ; s++) {
		for (int p = 0 ; p < curr_leaves.size() ; p ++) Q[s] += phenotype_val[curr_leaves[p]][s];
		Q[s] /= curr_leaves.size();
	}
}

void quantify_data::quantify_by_pca(vector < int > & curr_leaves, vector < float > & Q, int ipc) {
	vrb.bullet("pca on " + stb.str(curr_leaves.size()));
	pca P(sample_count, curr_leaves.size());
	P.fill(phenotype_val, curr_leaves);
	P.run(false, false, false);
	P.get(ipc, Q);
}

void quantify_data::quantify_by_loo(vector < int > & curr_leaves, vector < float > & Q, bool variance) {
	vrb.bullet("loo on " + stb.str(curr_leaves.size()));
	for (int s = 0 ; s < sample_count ; s ++) {
		//Initialization
		Q[s] = 0;
		//Compute means and sds
		vector < double > mes = vector < double > (curr_leaves.size(), 0.0);
		vector < double > sds = vector < double > (curr_leaves.size(), 0.0);
		for (int i = 0 ; i < curr_leaves.size() ; i ++) {
			for (int t = 0 ; t < sample_count ; t ++) if (t != s) mes[i] += phenotype_val[curr_leaves[i]][t];
			mes[i] /= sample_count;
			for (int t = 0 ; t < sample_count ; t ++) if (t != s) sds[i] += (phenotype_val[curr_leaves[i]][t] - mes[i]) * (phenotype_val[curr_leaves[i]][t] - mes[i]);
			sds[i] /= (sample_count - 1);
			sds[i] = sqrt(sds[i]);
		}
		//Compute correlations
		vector < double > cors;
		for (int i = 0 ; i < curr_leaves.size() - 1 ; i ++) {
			for (int j = i + 1 ; j < curr_leaves.size() ; j ++) {
				//
				cors.push_back(0.0);
				for (int t = 0 ; t < sample_count ; t ++) if (t != s) cors.back() += ((phenotype_val[curr_leaves[i]][t] - mes[i])/sds[i]) * ((phenotype_val[curr_leaves[j]][t] - mes[j])/sds[j]);
				cors.back() /= (sample_count - 1);
			}
		}
		double mean = 0;
		for (int c = 0 ; c < cors.size() ; c++) mean += cors[c];
		mean /= cors.size();
		if (!variance) Q[s] = mean;
		else {
			for (int c = 0 ; c < cors.size() ; c++) Q[s] += (cors[c] - mean) * (cors[c] - mean);
			Q[s] /= (cors.size() - 1);
		}
	}
}
