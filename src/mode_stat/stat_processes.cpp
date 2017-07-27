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

void stat_data::computeStatistics(string fout) {
	vrb.title("Compute basic statistics and write [" + fout + "]");
	output_file fd(fout.c_str());
	for (int i = 0 ; i < phenotype_count ; i ++) {
		double mean = getMean(phenotype_val[i]);
		double var = getVariance(phenotype_val[i], mean);
		double scaled_mean = getMean(phenotype_val[i], (phenotype_end[i]-phenotype_start[i]));
		double scaled_var = getVariance(phenotype_val[i], scaled_mean, (phenotype_end[i]-phenotype_start[i]));
		fd << i + 1;
		fd << " " << phenotype_chr[i];
		fd << " " << phenotype_start[i];
		fd << " " << phenotype_end[i];
		fd << " " << phenotype_id[i];
		fd << " " << mean;
		fd << " " << sqrt(var);
		fd << " " << scaled_mean;
		fd << " " << sqrt(scaled_var);
		fd << endl;
	}
	fd.close();
}
