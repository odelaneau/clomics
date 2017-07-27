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

#include "topo_data.h"

void topo_data::writeTree(string fout) {
	vrb.title("Write output tree in [" + fout + "]");
	output_file fd(fout);
	fd << header[0];
	for (int f = 1 ; f < header.size() ; f ++) fd << " " << header[f];
	fd << " BR_RATIO DISPER COMPL LIDX RIDX N_SIB N_REP" << endl;
	for (int n = 0 ; n  < vec_nodes.size() ; n ++) {
		fd << vec_nodes[n]->fields[0];
		for (int f = 1 ; f < vec_nodes[n]->fields.size() ; f ++) fd << " " << vec_nodes[n]->fields[f];
		fd << " " << vec_nodes[n]->branching_ratio << " " << vec_nodes[n]->dispersion << " " << vec_nodes[n]->completeness;
		fd << " " << vec_nodes[n]->left_most_index + 1 << " " << vec_nodes[n]->right_most_index + 1 << " " << vec_nodes[n]->n_siblings << " " << vec_nodes[n]->n_rep;
		fd << endl;
	}
	fd.close();
}
