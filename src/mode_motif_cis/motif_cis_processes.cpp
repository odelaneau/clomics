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

#include "motif_cis_data.h"

void motif_cis_data::computeCorrelations(int window) {
	vrb.title("Compute all possible correlations");
	vrb.bullet("n=" + stb.str(phenotype_count * (phenotype_count - 1) / 2));

	for (int i = 0 ; i < phenotype_count - 1 ; i ++) {
		for (int j = i + 1 ; j < phenotype_count ; j ++) {
			int posi = (phenotype_start[i] + phenotype_end[i]) / 2;
			int posj = (phenotype_start[j] + phenotype_end[j]) / 2;
			if (phenotype_start[j] > phenotype_end[i] && abs(posi - posj) <= window) {
				double curr_corr = abs(getCorrelation(phenotype_val[i], phenotype_val[j]));
				for (int im = 0 ; im < phenotype_motives[i].size() ; im ++) {
					for (int jm = 0 ; jm < phenotype_motives[j].size() ; jm ++) {
						int idx0 = phenotype_motives[i][im];
						int idx1 = phenotype_motives[j][jm];
						if (idx1 < idx0) {
							int tmp = idx1;
							idx1 = idx0;
							idx0 = tmp;
						}
						interaction_cor[idx0][idx1] += curr_corr;
						interaction_cnt[idx0][idx1] ++;
					}
				}
			}
			vrb.progress((i+1) * 1.0 / (phenotype_count - 1));
		}
	}
}
