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

#include "locate_data.h"

void locate_data::readPhenotypes(vector < string > & fbed) {

	vector < string > tokens;

	vrb.title("Reading phenotype data in " + stb.str(fbed.size()) + " files");
	n_anno = fbed.size();

	for ( int f = 0 ; f < fbed.size() ; f ++) {
		vrb.bullet("File [" + fbed[f] + "]");
		int n_found = 0;
		htsFile *fp = hts_open(fbed[f].c_str(),"r");
		if (!fp) vrb.error("Cannot open file");
		tbx_t *tbx = tbx_index_load(fbed[f].c_str());
		if (!tbx) vrb.error("Cannot open index file");
		kstring_t str = {0,0,0};
		while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
			if (str.l && str.s[0] != tbx->conf.meta_char) {
				stb.split(string(str.s), tokens);
				if (tokens.size() < 6) vrb.error("Incorrect number of columns!");
				int start = atoi(tokens[1].c_str());
				int end = atoi(tokens[2].c_str());
				unordered_map < string, lcluster * >::iterator it_map = map_leaves.find(tokens[3]);
				if (it_map != map_leaves.end()) {
					it_map->second->regulatory_elements.push_back(regulatory_element(start, end, fbed.size()));
					it_map->second->regulatory_elements.back().peak_count[f] ++;
					it_map->second->distance = 0;
					it_map->second->start = start;
					it_map->second->end = end;
					it_map->second->length = end - start + 1;
					n_found ++;
				}
			}
		}
		tbx_destroy(tbx);
		if (hts_close(fp)) vrb.error("Cannot properly close file");
	}

	int n_unfound = 0, n_found = 0;
	for (unordered_map < string, lcluster * >::iterator it_map2 = map_leaves.begin() ; it_map2 != map_leaves.end() ; it_map2 ++) {
		if (it_map2->second->regulatory_elements.size() == 0) n_unfound ++;
		else n_found ++;
	}
	if (n_unfound == 0) vrb.bullet(stb.str(map_leaves.size()) + " leaves positions set up");
	else vrb.error(stb.str(n_unfound) + " leaves positions NOT set up");
}
