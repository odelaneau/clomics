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

#include "mbed_data.h"

void mbed_data::writeBED(string fout) {
	vrb.title("Write output BED in [" + fout + "]");
	output_file fd(fout);
	for (int n = 0 ; n  < vec_nodes.size() ; n ++) {
		if (vec_nodes[n]->leaf() && vec_nodes[n]->idx_module >= 0) {
			fd << vec_nodes[n]->start << " " << vec_nodes[n]->end;
			fd << vec_nodes[n]->sid << " " << vec_nodes[vec_nodes[n]->idx_module]->sid << endl;
		}
	}
	fd.close();
}

