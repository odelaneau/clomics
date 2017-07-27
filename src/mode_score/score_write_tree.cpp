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

#include "score_data.h"

void score_data::writeTree(string fout) {
	vrb.title("Write output in [" + fout + "]");
	output_file fd(fout);
	fd << header[0];
	for (int f = 1 ; f < header.size() ; f ++) fd << " " << header[f];
	fd << " BOOT JACC" << endl;
	for (int n = 0 ; n  < vec_nodes.size() ; n ++) {
		fd << vec_nodes[n]->fields[0];
		for (int f = 1 ; f < vec_nodes[n]->fields.size() ; f ++) fd << " " << vec_nodes[n]->fields[f];
		fd << " " << vec_nodes[n]->bootstrap << " " << vec_nodes[n]->jaccard << endl;
	}
	fd.close();
}

