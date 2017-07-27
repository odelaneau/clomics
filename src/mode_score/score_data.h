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

#ifndef _SCORE_DATA_H
#define _SCORE_DATA_H

//INCLUDES
#include "../common/data.h"

class scluster {
public:
	string sid;
	int index;
	int count;
	vector < string > fields;
	float bootstrap;
	float jaccard;
	scluster * child1;
	scluster * child2;

	scluster (string _sid, int _index, int _count, scluster * c1 = NULL, scluster * c2 = NULL) {
		sid = _sid;
		index = _index;
		count = _count;
		bootstrap = 0;
		jaccard = 0;
		child1 = c1;
		child2 = c2;
	};

	bool leaf() {
		return (child1 == NULL || child2 == NULL);
	}
};

class score_data : public data {
public:

	//TREES
	scluster * tree_nominal;
	vector < scluster * > tree_samples;
	vector < string > header;

	//INTERNAL DATA
	vector < scluster * > vec_nodes;
	unordered_map < string, scluster * > map_nodes;

	//CONSTRUCTOR / DESTRUCTOR
	score_data() {};
    ~score_data() {};

	//READ OR WRITE DATA
    void readNominalTree(string);
    void readSampledTrees(vector < string > &);

    void writeTree(string);

    //SCORE1: OVERLAP ROUTINES
    void computeOverlapScore(scluster *, set < int > &, bool);
    void computeOverlapScore();

    //SCORE2: JACCARD ROUTINES
    int computeJaccardScore(scluster * , vector < bool > & , int, float &);
    void computeJaccardScore(scluster * , vector < bool > & );
    void computeJaccardScore();
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void score_main(vector < string > &);

#endif
