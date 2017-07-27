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

#ifndef _LOCATE_DATA_H
#define _LOCATE_DATA_H

#define DEFAULT_MERGE_DISTANCE 50

//INCLUDES
#include "../common/data.h"

class regulatory_element {
public:
	static int merge_distance;
	int start, end;
	vector < int > peak_count;

	regulatory_element(int _start, int _end, int _size) {
		start = _start;
		end = _end;
		peak_count = vector < int > (_size, 0);
	}

	regulatory_element(regulatory_element & re1, regulatory_element & re2) {
		start = min(re1.start, re2.start);
		end = max(re1.end, re2.end);;
		peak_count = re1.peak_count;
		for (int e = 0 ; e < peak_count.size() ; e ++) peak_count[e] += re2.peak_count[e];
	}

	~regulatory_element() {
		peak_count.clear();
	}

	bool overlap(regulatory_element & re) {
		return (re.start <= (end + merge_distance) && start <= (re.end + merge_distance));
	}

	int encode() {
		int toint = 0;
		for (int a = 0 ; a < peak_count.size() ; a ++) if (peak_count[a]) toint += 1 << a;
		return toint;
	}

	bool operator < (const regulatory_element & d) const {
		if (start < d.start) return true;
		if (start > d.start) return false;
		if (end < d.end) return true;
		else return false;
	}
};

class lcluster {
public:
	string sid;
	int index;
	int count;
	vector < string > fields;
	vector < regulatory_element > regulatory_elements;
	int distance, start, end, length;
	lcluster * child1;
	lcluster * child2;

	lcluster (string _sid, int _index, int _count, lcluster * c1 = NULL, lcluster * c2 = NULL) {
		sid = _sid;
		index = _index;
		count = _count;
		distance = 0;
		start = 0;
		end = 0;
		length = 0;
		child1 = c1;
		child2 = c2;
	};

	bool leaf() {
		return (child1 == NULL || child2 == NULL);
	}

	int regulatory_size () {
		int tmp = 0;
		for (int e = 0 ; e < regulatory_elements.size() ; e++) {
			tmp += regulatory_elements[e].end -regulatory_elements[e].start + 1;
		}
		return tmp;
	}

	void compositionALL(vector < float > & C) {
		float sum = 0.0;
		int n_types = regulatory_elements[0].peak_count.size();
		for (int e = 0 ; e < regulatory_elements.size() ; e ++) for (int p = 0 ; p < n_types ; p++) {
			C[p] += regulatory_elements[e].peak_count[p];
			sum += regulatory_elements[e].peak_count[p];
		}
		for (int p = 0 ; p < n_types ; p++) C[p] /= sum;
	}

	void compositionIND(vector < int > & RE) {
		for (int e = 0 ; e < regulatory_elements.size() ; e ++) RE[regulatory_elements[e].encode()] ++;
	}
};

class locate_data : public data {
public:

	//TREE
	lcluster * tree;
	vector < lcluster * > vec_nodes;
	unordered_map < string, lcluster * > map_leaves;
	vector < string > header;
	int n_anno;

	//CONSTRUCTOR / DESTRUCTOR
	locate_data() {};
    ~locate_data() {};

	//READ OR WRITE DATA
    void readTree(string);
    void writeTree(string);
    void writeRegulatoryElements(string fout);
    void readPhenotypes(vector < string > &);

    //PROPAGATE REGULATORY ELEMENTS WITHIN TREE
    void propagateRegulatoryElements(lcluster *);
    void propagateRegulatoryElements();
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void locate_main(vector < string > &);

#endif

