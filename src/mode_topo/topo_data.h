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

#ifndef _TOPO_DATA_H
#define _TOPO_DATA_H

//INCLUDES
#include "../common/data.h"

class tcluster {
public:
	int index;
	tcluster * child1;
	tcluster * child2;
	string sid;
	int count;
	vector < string > fields;
	float branching_ratio, dispersion, completeness;
	int left_most_index, right_most_index, n_siblings, n_rep;


	tcluster (string _sid, int _index, int _count, tcluster * c1 = NULL, tcluster * c2 = NULL) {
		sid = _sid;
		index = _index;
		left_most_index = _index;
		right_most_index = _index;
		count = _count;
		child1 = c1;
		child2 = c2;
		branching_ratio = 0.0;
		dispersion = 0.0;
		completeness = 0.0;
		n_siblings = 0;
		n_rep = 0;
	};

	bool leaf() {
		return (child1 == NULL || child2 == NULL);
	}
};

class topo_data : public data {
public:

	//TREE
	tcluster * tree;
	vector < tcluster * > vec_nodes;
	vector < string > header;

	//CONSTRUCTOR / DESTRUCTOR
	topo_data() {};
    ~topo_data() {};

    //PROCESS
    void characterizeTopology();
    void characterizeTopology(tcluster *, vector < int > &);

	//READ OR WRITE DATA
    void readTree(string);
    void writeTree(string);

};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void topo_main(vector < string > &);

#endif
