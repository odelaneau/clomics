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

#include "build_data.h"

void build_data::imputePhenotypes() {
	for (int p = 0; p < phenotype_count ; p ++) {
		double mean = 0.0;
		int c_mean= 0;
		for (int s = 0; s < sample_count; s ++) {
			if (phenotype_val[p][s] != bcf_float_missing) {
				mean += phenotype_val [p][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (phenotype_val[p][s] == bcf_float_missing) phenotype_val[p][s] = mean;
	}
}

void build_data::normalizePhenotypes() {
	for (int x = 0 ; x < phenotype_count ; x++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += phenotype_val[x][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			phenotype_val[x][s] -= mean;
			sum += phenotype_val[x][s] * phenotype_val[x][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) phenotype_val[x][s] /= sum;
	}
}

bool build_data::setPhenotypeRegion(string reg) {
	return regionPhenotype.parse(reg);
}

void build_reorder(vector < float > & v, vector < int > & order)  {
	vector < float > v_tmp = v;
	for (int e = 0 ; e < v.size() ; e ++) v_tmp[e] = v[order[e]];
	v = v_tmp;
}

void build_data::bootstrapPhenotypes(float proportion) {
	vrb.title("Sampling with replacement samples [bootstrap]");
	vector < int > indexes = vector < int > (sample_count);
	for (int s = 0 ; s < sample_count ; s++) indexes[s] = rng.getInt(sample_count);
	for (int p = 0 ; p < phenotype_count ; p ++) build_reorder(phenotype_val[p], indexes);
	if (proportion < 1.0) {
		for (int p = 0 ; p < phenotype_count ; p ++) {
			phenotype_val[p].erase(phenotype_val[p].begin() + floor(sample_count * proportion) + 1, phenotype_val[p].end());
		}
		sample_count = phenotype_val[0].size();
	}
	vrb.bullet("Proportion = " + stb.str(proportion, 4) + " #samples = " + stb.str(sample_count));
}

void build_data::jackknifePhenotypes(float proportion) {
	vrb.title("Sampling without replacement individuals [jack-knife]");
	vector < int > indexes = vector < int > (sample_count);
	for (int s = 0 ; s < sample_count ; s++) indexes[s] = s;
	shuffle(indexes.begin(), indexes.end(), rng.getEngine());
	for (int p = 0 ; p < phenotype_count ; p ++) {
		build_reorder(phenotype_val[p], indexes);
		phenotype_val[p].erase(phenotype_val[p].begin() + floor(sample_count * proportion) + 1, phenotype_val[p].end());
	}
	sample_count = phenotype_val[0].size();;
	vrb.bullet("Proportion = " + stb.str(proportion, 4) + " #samples = " + stb.str(sample_count));
}
