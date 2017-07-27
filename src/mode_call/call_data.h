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

#ifndef _CALL_DATA_H
#define _CALL_DATA_H

#define F_EQ	0
#define F_GT	1
#define F_LT	2

//INCLUDES
#include "../common/data.h"

class ccluster {
public:
	int index;
	ccluster * child1;
	ccluster * child2;
	int count;
	double acorrelation, apvalue;
	double bcorrelation, bpvalue;
	double ccorrelation, cpvalue;
	double dcorrelation, dpvalue;
	double ecorrelation, epvalue;
	double pcavariance;
	int module;

	string sid;
	vector < string > fields;

	ccluster (string _sid, int _index, int _count, ccluster * c1 = NULL, ccluster * c2 = NULL) {
		sid = _sid;
		index = _index;
		count = _count;
		child1 = c1;
		child2 = c2;
		acorrelation = 0.0; apvalue = 1.0;
		bcorrelation = 0.0; bpvalue = 1.0;
		ccorrelation = 0.0; cpvalue = 1.0;
		dcorrelation = 0.0; dpvalue = 1.0;
		ecorrelation = 0.0; epvalue = 1.0;
		pcavariance = 0.0;
		module = 0;
	};

	bool leaf() {
		return (child1 == NULL || child2 == NULL);
	}
};

class call_data : public data {
public:
	//TREE
	ccluster * tree;
	vector < ccluster * > vec_nodes;
	vector < string > header;

	//CONSTRUCTOR / DESTRUCTOR
	call_data() {};
    ~call_data() {};

    //PROCESS
    void parseFiltersAndAnnotations(vector < string > &, vector < string > &);

    void tagSignificant(ccluster *, int, float, float, vector < int > &);

    void writeTree(string fout);
    void readTree(string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void call_main(vector < string > &);

#endif
