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

void topo_data::readTree(string ftxt) {
	string buffer;
	vector < string > tokens;

	vrb.title("Reading tree in [" + ftxt + "]");

	//Sweep1: reading nodes
	input_file fd1(ftxt);
	getline(fd1, buffer);
	stb.split(buffer, tokens);
	for (int t = 0 ; t < tokens.size() ; t ++) header.push_back(tokens[t]);
	while (getline(fd1, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() < 5) vrb.error("Incorrect number of columns in file, should be >= 5");
		if (!stb.numeric(tokens[0])) vrb.error("Non-numeric node index [" + tokens[0] + "]");
		if (tokens[1] != "NA" && !stb.numeric(tokens[1])) vrb.error("Non-numeric child1 index [" + tokens[1] + "]");
		if (tokens[2] != "NA" && !stb.numeric(tokens[2])) vrb.error("Non-numeric child2 index [" + tokens[2] + "]");
		if (!stb.numeric(tokens[4])) vrb.error("Non-numeric node count [" + tokens[4] + "]");
		int index = atoi(tokens[0].c_str()) - 1;
		int count = atoi(tokens[4].c_str());
		tcluster * new_cluster = new tcluster (tokens[3], index, count);
		for (int t = 0 ; t < tokens.size() ; t ++) new_cluster->fields.push_back(tokens[t]);
		if (index < 0) vrb.error("Incorrect node index [" + tokens[0] + "]");
		else if (index >= vec_nodes.size()) vec_nodes.resize(index + 1, NULL);
		else if (vec_nodes[index] != NULL) vrb.error("Duplicate node index [" + tokens[0] + "]");
		vec_nodes[index] = new_cluster;
	}
	vrb.bullet("#nodes = " + stb.str(vec_nodes.size()));
	fd1.close();

	//Check continuity in indexes
	for (int n = 0 ; n < vec_nodes.size() ; n ++) if (vec_nodes[n] == NULL) vrb.error("Non-contiguous indexes: [" + stb.str(n+1) + "] is missing");

	//Sweep2: connecting nodes
	input_file fd2(ftxt);
	getline(fd2, buffer);
	int n_leaves = 0, n_internals = 0;
	while (getline(fd2, buffer)) {
		stb.split(buffer, tokens);
		int index_pp = atoi(tokens[0].c_str()) - 1;
		if (tokens[1] != "NA" || tokens[2] != "NA") {
			int index_c1 = atoi(tokens[1].c_str()) - 1;
			int index_c2 = atoi(tokens[2].c_str()) - 1;
			if (index_c1 < 0 || index_c1 >= vec_nodes.size()) vrb.error("Incorrect child1 node index [" + tokens[1] + "]");
			if (index_c2 < 0 || index_c2 >= vec_nodes.size()) vrb.error("Incorrect child2 node index [" + tokens[2] + "]");
			vec_nodes[index_pp]->child1 = vec_nodes[index_c1];
			vec_nodes[index_pp]->child2 = vec_nodes[index_c2];
			n_internals ++;
		} else n_leaves ++;
	}
	tree = vec_nodes.back();
	vrb.bullet("#internal_nodes = " + stb.str(n_internals));
	vrb.bullet("#leaf_nodes = " + stb.str(n_leaves));
	if (tree->count != n_leaves) vrb.error("Root node count [" + stb.str(tree->count) + "] and number of leaves [" + stb.str(n_leaves) + "] are different");
	fd2.close();
}
