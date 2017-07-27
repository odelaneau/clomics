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

#ifndef _MERGE_DATA_H
#define _MERGE_DATA_H

//INCLUDES
#include "../common/data.h"

#define DEFAULT_MERGE_DISTANCE 50

class merge_phenotype_data {
public:
	string chr, id;
	unsigned int start, end;

	merge_phenotype_data(string _id, string _chr, unsigned int _start, unsigned int _end) {
		id = _id;
		chr = _chr;
		start = _start;
		end = _end;
	}

	bool operator < (const merge_phenotype_data & d) const {
		if (chr < d.chr) return true;
		if (chr > d.chr) return false;
		if (start < d.start) return true;
		if (start > d.start) return false;
		if (end < d.end) return true;
		else return false;
	}
};

class merge_data : public data {
public:

	//PHENOTYPES STRUCT TO HANDLE MULTIPLE FILES
	unordered_map < string, unsigned int > phenotype_map;
	vector < merge_phenotype_data > phenotype_vec;

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < int > phenotype_assay;						//phenotype assay
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
	vector < bool > phenotype_neg;						//phenotype on the negative strand?

	vector < vector < int > > phenotype_grp;


	//CONSTRUCTOR / DESTRUCTOR
	merge_data() {};
    ~merge_data() {};

	//READ OR GENERATE DATA
    void scanPhenotypes(string);
    void readPhenotypes(string, int);
    void sortPhenotypes();

	//ANALYSIS
	void collapse(int);

	//OUTPUT
	void writePhenotypes(string, int, bool);

	//ROUTINES
	void normalTransform(vector < float > & V);
	int getStart(int p) {
		int min = 1000000000;
		for (int i = 0 ; i < phenotype_grp[p].size() ; i++) {
			if (phenotype_start[phenotype_grp[p][i]] < min) min = phenotype_start[phenotype_grp[p][i]];
			if (phenotype_end[phenotype_grp[p][i]] < min) min = phenotype_end[phenotype_grp[p][i]];
		}
		return min;
	}

	int getEnd(int p) {
		int max = 0;
		for (int i = 0 ; i < phenotype_grp[p].size() ; i++) {
			if (phenotype_start[phenotype_grp[p][i]] > max) max = phenotype_start[phenotype_grp[p][i]];
			if (phenotype_end[phenotype_grp[p][i]] > max) max = phenotype_end[phenotype_grp[p][i]];
		}
		return max;
	}

	void getValues(int p, vector < float > & values) {
		assert(values.size() == sample_count);
		for (int i = 0 ; i < sample_count ; i ++) {
			values[i] = 0.0;
			for (int j = 0 ; j < phenotype_grp[p].size() ; j ++) values[i] += phenotype_val[phenotype_grp[p][j]][i];
		}
	}

	void getAssays(int p, vector < int > & assays) {
		assays = vector < int > (assays.size(), 0);
		for (int i = 0 ; i < phenotype_grp[p].size() ; i ++) assays[phenotype_assay[phenotype_grp[p][i]]] ++;
	}
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void merge_main(vector < string > &);

#endif
