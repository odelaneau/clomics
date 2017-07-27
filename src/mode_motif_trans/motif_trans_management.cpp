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

#include "motif_trans_data.h"

void motif_trans_data::imputePhenotypes1() {
	for (int p = 0; p < phenotype_count1 ; p ++) {
		double mean = 0.0;
		int c_mean= 0;
		for (int s = 0; s < sample_count; s ++) {
			if (phenotype_val1[p][s] != bcf_float_missing) {
				mean += phenotype_val1 [p][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (phenotype_val1[p][s] == bcf_float_missing) phenotype_val1[p][s] = mean;
	}
}

void motif_trans_data::normalizePhenotypes1() {
	for (int x = 0 ; x < phenotype_count1 ; x++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += phenotype_val1[x][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			phenotype_val1[x][s] -= mean;
			sum += phenotype_val1[x][s] * phenotype_val1[x][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) phenotype_val1[x][s] /= sum;
	}
}

void motif_trans_data::imputePhenotypes2() {
	for (int p = 0; p < phenotype_count2 ; p ++) {
		double mean = 0.0;
		int c_mean= 0;
		for (int s = 0; s < sample_count; s ++) {
			if (phenotype_val2[p][s] != bcf_float_missing) {
				mean += phenotype_val2 [p][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (phenotype_val2[p][s] == bcf_float_missing) phenotype_val2[p][s] = mean;
	}
}

void motif_trans_data::normalizePhenotypes2() {
	for (int x = 0 ; x < phenotype_count2 ; x++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += phenotype_val2[x][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			phenotype_val2[x][s] -= mean;
			sum += phenotype_val2[x][s] * phenotype_val2[x][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) phenotype_val2[x][s] /= sum;
	}
}

bool motif_trans_data::setPhenotypeRegion1(string reg) {
	return regionPhenotype1.parse(reg);
}

bool motif_trans_data::setPhenotypeRegion2(string reg) {
	return regionPhenotype2.parse(reg);
}
