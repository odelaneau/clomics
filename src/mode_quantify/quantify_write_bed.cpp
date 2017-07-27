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

#include "quantify_data.h"

void quantify_data::writePhenotypesTrees(string fbed, bool normal) {
	string buffer;
	vector < string > tokens;

	vrb.title("Writing quantifications in [" + fbed + "]");
	output_file fd(fbed);
	vrb.bullet("writing header");
	fd << "#chr\tstart\tend\tid\tdummy\tstrand";
	for (int s = 0 ; s < sample_count ; s ++) fd << "\t" << sample_id[s];
	fd << endl;
	vrb.bullet("writing body");

	for (int n = 0 ; n < vec_nodes.size() ; n ++) {
		if (!vec_nodes[n]->leaf() && (modules.size() == 0 || modules.count(vec_nodes[n]->sid) > 0)) {
			if (normal) normalTransform(vec_quantifications[vec_nodes[n]->index]);
			fd << phenotype_chr[vec_nodes[n]->first] << "\t" << vec_nodes[n]->start << "\t" << vec_nodes[n]->end;
			fd << "\t" << vec_nodes[n]->sid << "\t" << vec_nodes[n]->count << "\t" << (phenotype_neg[vec_nodes[n]->first]?"-":"+");
			for (int s = 0 ; s < sample_count ; s ++) fd << "\t" << vec_quantifications[vec_nodes[n]->index][s];
			fd << endl;
		}
	}
	fd.close();
}

void quantify_data::writePhenotypesGroups(string fbed, bool normal) {
	string buffer;
	vector < string > tokens;

	vrb.title("Writing quantifications in [" + fbed + "]");
	output_file fd(fbed);
	vrb.bullet("writing header");
	fd << "#chr\tstart\tend\tid\tdummy\tstrand";
	for (int s = 0 ; s < sample_count ; s ++) fd << "\t" << sample_id[s];
	fd << endl;
	vrb.bullet("writing body");
	for (int g = 0 ; g < vec_groups.size() ; g ++) {
		if (normal) normalTransform(vec_quantifications[g]);
		string chr = phenotype_chr[vec_groups[g][0]];
		int start = 1000000000;
		int end = 0;
		for (int p = 0 ; p < vec_groups[g].size() ; p ++) if (phenotype_start[vec_groups[g][p]] < start) start = phenotype_start[vec_groups[g][p]];
		for (int p = 0 ; p < vec_groups[g].size() ; p ++) if (phenotype_end[vec_groups[g][p]] > end) end = phenotype_end[vec_groups[g][p]];
		string id = phenotype_id[vec_groups[g][0]];
		for (int p = 1 ; p < vec_groups[g].size() ; p ++) {
			id += "_";
			id += phenotype_id[vec_groups[g][p]];
		}
		fd << chr << "\t" << start << "\t" << end << "\t" << id <<"\t" << vec_groups[g].size() << "\t+";
		for (int s = 0 ; s < sample_count ; s ++) fd << "\t" << vec_quantifications[g][s];
		fd << endl;
	}
	fd.close();
}
