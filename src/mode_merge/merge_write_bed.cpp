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

void merge_data::writePhenotypes(string fbed, int n_assays, bool normal) {
	string buffer;
	vector < string > tokens;

	vrb.title("Writing quantifications in [" + fbed + "]");
	output_file fd(fbed);
	vrb.bullet("writing header");
	fd << "#chr\tstart\tend\tid\tdummy\tstrand";
	for (int s = 0 ; s < sample_count ; s ++) fd << "\t" << sample_id[s];
	fd << endl;
	vrb.bullet("writing body");

	vector < float > values = vector < float > (sample_count, 0.0);
	vector < int > assays = vector < int >(n_assays, 0);
	for (int p = 0 ; p < phenotype_grp.size() ; p ++) {
		int start = getStart(p);
		int end = getEnd(p);
		getValues(p, values);
		getAssays(p, assays);
		if (normal) normalTransform(values);
		fd << phenotype_chr[phenotype_grp[p][0]] << "\t" << start << "\t" << end;
		if (phenotype_grp[p].size() > 1) fd << "\tRE" << p << "\t" << assays[0];
		else fd << "\t" << phenotype_id[phenotype_grp[p][0]] << "\t" << assays[0];
		for (int a = 1 ; a < assays.size() ; a ++) fd << "," << assays[a];
		fd << "\t" << (phenotype_neg[phenotype_grp[p][0]]?"-":"+");
		for (int s = 0 ; s < sample_count ; s ++) fd << "\t" << values[s];
		fd << endl;
	}
	fd.close();
}
