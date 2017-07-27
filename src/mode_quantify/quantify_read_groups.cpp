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

void quantify_data::readGroups(string ftxt) {
	string buffer;
	vector < string > tokens;
	unordered_map < string, unsigned int >::iterator it_map;

	vrb.title("Reading groups in [" + ftxt + "]");
	int n_found = 0, n_total = 0;
	input_file fd(ftxt);
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() < 2) vrb.error("Incorrect number of columns in file, should be >= 2");
		vector < int > curr_grp = vector < int > ();
		for (int t = 0 ; t < tokens.size() ; t ++) {
			it_map = phenotype_map.find(tokens[t]);
			if (it_map != phenotype_map.end()) curr_grp.push_back(it_map->second);
		}
		if (curr_grp.size() > 1) {
			vec_groups.push_back(curr_grp);
			n_found ++;
		} else n_total++;
	}
	vrb.bullet("#found = " + stb.str(n_found));
	vrb.bullet("#total = " + stb.str(n_total));
	fd.close();
}
