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

void locate_data::writeTree(string fout) {
	int n_comb = 1 << n_anno;
	vrb.title("Write TREE in [" + fout + "]");
	output_file fd(fout);
	fd << header[0];
	for (int f = 1 ; f < header.size() ; f ++) fd << " " << header[f];
	fd << " DIST START STOP SIZE NRE";
	for (int a = 0 ; a < n_anno ; a ++) fd << " ANN" << a+1;
	for (int c = 0 ; c < n_comb ; c ++) fd << " CMB" << c+1;
	fd << endl;
	for (int n = 0 ; n  < vec_nodes.size() ; n ++) {
		fd << vec_nodes[n]->fields[0];
		for (int f = 1 ; f < vec_nodes[n]->fields.size() ; f ++) fd << " " << vec_nodes[n]->fields[f];
		fd << " " << vec_nodes[n]->distance << " " << vec_nodes[n]->start << " " << vec_nodes[n]->end << " " << vec_nodes[n]->length << " " << vec_nodes[n]->regulatory_size() << " " << vec_nodes[n]->regulatory_elements.size();
		vector < float > peak_compositionALL = vector < float > (n_anno, 0.0);
		vector < int > peak_compositionIND = vector < int > (n_comb, 0);
		vec_nodes[n]->compositionALL(peak_compositionALL);
		vec_nodes[n]->compositionIND(peak_compositionIND);
		for (int p = 0 ; p < peak_compositionALL.size() ; p ++) fd << " " << peak_compositionALL[p];
		for (int p = 0 ; p < peak_compositionIND.size() ; p ++) fd << " " << peak_compositionIND[p];
		fd << endl;
	}
	fd.close();
}

void locate_data::writeRegulatoryElements(string fout) {
	vrb.title("Write REs in [" + fout + "]");
	output_file fd(fout);
	fd << "NODE";
	for (int a = 0 ; a < n_anno ; a ++) fd << " ANN" << a+1;
	fd << endl;
	for (int n = 0 ; n  < vec_nodes.size() ; n ++) {
		for (int re = 0 ; re < vec_nodes[n]->regulatory_elements.size() ; re ++) {
			fd << vec_nodes[n]->fields[3];
			for (int a = 0 ; a < n_anno ; a ++) fd << " " << vec_nodes[n]->regulatory_elements[re].peak_count[a];
			fd << endl;
		}
	}
	fd.close();
}

