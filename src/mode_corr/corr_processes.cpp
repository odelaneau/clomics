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

#include "corr_data.h"

void corr_data::computeCorrelations(string fout) {
	vrb.title("Compute all possible correlations and write [" + fout + "]");
	vrb.bullet("n=" + stb.str(phenotype_count * (phenotype_count - 1) / 2));
	output_file fd(fout.c_str());
	for (int i = 0 ; i < phenotype_count - 1 ; i ++) {
		for (int j = i + 1 ; j < phenotype_count ; j ++) {
			double corr = getCorrelation(phenotype_val[i], phenotype_val[j]);
			fd << i + 1 << " " << j + 1;
			fd << " " << (int)(phenotype_start[i]+phenotype_end[i])/2 << " " << (int)(phenotype_start[j]+phenotype_end[j])/2;
			fd << " " << phenotype_id[i] << " " << phenotype_id[j];
			fd << " " << stb.str(corr, 3);
			fd << " " << getPvalue(corr, sample_count - 2);
			if (phenotype_paired.size() > 0) fd << " " << phenotype_paired[i][j];
			fd << endl;
		}
		vrb.progress((i+1) * 1.0 / (phenotype_count - 1));
	}
	fd.close();
}
