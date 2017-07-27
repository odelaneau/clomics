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

#include "hic_data.h"

void hic_data::computeInteractions(string fout) {
	vrb.title("Report all possible interactions in [" + fout + "]");
	output_file fd(fout.c_str());
	for (int i = 0 ; i < phenotype_count - 1 ; i ++) {
		for (int j = i + 1 ; j < phenotype_count ; j ++) {
			if (phenotype_hic_cnt[i][j] > 0.0) {
				double corr = getCorrelation(phenotype_val[i], phenotype_val[j]);
				fd << i + 1 << " " << phenotype_chr[i] << " " << phenotype_start[i] << " " << phenotype_end[i] << " " << phenotype_id[i];
				fd << " " << j + 1 << " " << phenotype_chr[j] << " " << phenotype_start[j] << " " << phenotype_end[j] << " " << phenotype_id[j];
				fd << " " << stb.str(phenotype_hic_val[i][j]/phenotype_hic_cnt[i][j], 3) << " " << stb.str(corr, 3)  << " " << getPvalue(corr, sample_count-2) << endl;
			} else if (rng.getDouble() <= thinin) {
				double corr = getCorrelation(phenotype_val[i], phenotype_val[j]);
				fd << i + 1 << " " << phenotype_chr[i] << " " << phenotype_start[i] << " " << phenotype_end[i] << " " << phenotype_id[i];
				fd << " " << j + 1 << " " << phenotype_chr[j] << " " << phenotype_start[j] << " " << phenotype_end[j] << " " << phenotype_id[j];
				fd << " 0 " << stb.str(corr, 3)  << " " << getPvalue(corr, sample_count-2) << endl;
			}
		}
		vrb.progress((i+1) * 1.0 / phenotype_count);
	}
	fd.close();
}
