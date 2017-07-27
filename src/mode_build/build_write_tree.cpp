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

#include "build_data.h"

void build_data::writeTree(string fout) {
	vrb.title("Write output in [" + fout + "]");
	output_file fd(fout);
	vector < bcluster * > nodes;
	tree->vectorize(nodes);
	if (!resampling) fd << "IDX CHILD1 CHILD2 UUID COUNT ACORR APVAL BCORR BPVAL CCORR CPVAL DCORR DPVAL ECORR EPVAL PCA1" << endl;
	for (int n = 0 ; n  < nodes.size() ; n ++) {
		fd << nodes[n]->index + 1;
		fd << " " << ((nodes[n]->child1 != NULL)?stb.str(nodes[n]->child1->index + 1):"NA");
		fd << " " << ((nodes[n]->child2 != NULL)?stb.str(nodes[n]->child2->index + 1):"NA");
		fd << " " << nodes[n]->sid << " " << nodes[n]->count;
		if (!resampling) {
			fd << " " << nodes[n]->acorrelation << " " << nodes[n]->apvalue;
			fd << " " << nodes[n]->bcorrelation << " " << nodes[n]->bpvalue;
			fd << " " << nodes[n]->ccorrelation << " " << nodes[n]->cpvalue;
			fd << " " << nodes[n]->dcorrelation << " " << nodes[n]->dpvalue;
			fd << " " << nodes[n]->ecorrelation << " " << nodes[n]->epvalue;
			fd << " " << nodes[n]->pcavariance;
		}
		fd << endl;
	}
	fd.close();
}

