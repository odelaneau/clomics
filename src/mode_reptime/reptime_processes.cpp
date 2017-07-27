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

#include "reptime_data.h"

void reptime_data::builtRTtree() {
	vrb.title("Building interval tree for replication timings");
	vector < Interval < float > > interval_vec;
	for (int r = 1 ; r < reptime_pos.size() ; r ++) interval_vec.push_back(Interval < float > (reptime_pos[r-1], reptime_pos[r], (reptime_val[r-1]+reptime_val[r])/2));
	reptime_tree = IntervalTree < float > (interval_vec);
	vrb.bullet("done!");
}

void reptime_data::project(string fout) {
	vrb.title("Report in file [" + fout + "] the extrapolated replication times for all peaks");
	output_file fd(fout.c_str());
	for (int p = 0 ; p < phenotype_count; p ++) {
		vector < Interval < float > > overlaps;
		reptime_tree.findOverlapping(phenotype_start[p], phenotype_end[p], overlaps);
		double meanRT = 0.0;
		for (int o = 0 ; o < overlaps.size() ; o ++) meanRT += overlaps[o].value;
		fd << phenotype_chr[p] << "\t" << phenotype_start[p] << "\t" << phenotype_end[p] << "\t" << phenotype_id[p];
		fd << "\t" << ((overlaps.size() == 0)?"NA":stb.str(meanRT/overlaps.size(), 3));
		fd <<endl;
		vrb.progress((p+1) * 1.0 / phenotype_count);
	}
	fd.close();
}
