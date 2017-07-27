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

#include "merge_data.h"


void merge_data::normalTransform(vector < float > & V) {
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



void merge_data::collapse(int merge_dist) {
	bool merge_done = false;

	vrb.title("Collapsing");
	int n_merge = 0;
	while (!merge_done) {
		int idx = -1;
		for (int p = 1 ; p < phenotype_grp.size() && idx < 0 ; p ++) {
			string chr0 = phenotype_chr[p-1];
			string chr1 = phenotype_chr[p-0];
			int start0 = getStart(p-1);
			int start1 = getStart(p-0);
			int end0 = getEnd(p-1);
			int end1 = getEnd(p-0);
			if (chr0 == chr1 && start0 < (end1 + merge_dist) && start1 < (end0 + merge_dist)) idx = p;
		}
		if (idx >= 0) {
			phenotype_grp[idx-1].insert(phenotype_grp[idx-1].end(),phenotype_grp[idx].begin(), phenotype_grp[idx].end());
			phenotype_grp.erase(phenotype_grp.begin() + idx);
			vrb.bullet("#REs = " + stb.str(phenotype_grp.size()));
			n_merge ++;
		} else merge_done = true;
	}
	vrb.bullet("#merging = " + stb.str(n_merge));
}
