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

#ifndef _MBED_DATA_H
#define _MBED_DATA_H

#define F_EQ	0
#define F_GT	1
#define F_LT	2

//INCLUDES
#include "../common/data.h"

class mcluster {
public:
	int index;
	mcluster * child1;
	mcluster * child2;
	int count;
	int idx_module;
	int start, end;

	string sid;
	vector < string > fields;

	mcluster (string _sid, int _index, int _count, mcluster * c1 = NULL, mcluster * c2 = NULL) {
		sid = _sid;
		index = _index;
		count = _count;
		child1 = c1;
		child2 = c2;
		idx_module = -1;
		start = -1;
		end = -1;
	};

	bool leaf() {
		return (child1 == NULL || child2 == NULL);
	}
};

class mbed_data : public data {
public:
	//TREE
	mcluster * tree;
	vector < mcluster * > vec_nodes;
	vector < string > header;

	set < string > module_set;

	//CONSTRUCTOR / DESTRUCTOR
	mbed_data() {};
    ~mbed_data() {};

    //PROCESS
    void readTree(string);
    void readModules(string);
    void tagModuled(mcluster * , int );
    void writeBED(string fout);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void mbed_main(vector < string > &);

#endif
