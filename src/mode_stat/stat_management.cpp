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

#include "stat_data.h"

void stat_data::imputePhenotypes() {
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

bool stat_data::setPhenotypeRegion(string reg) {
	return regionPhenotype.parse(reg);
}
