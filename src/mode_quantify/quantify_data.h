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

#ifndef _QUANTIFY_DATA_H
#define _QUANTIFY_DATA_H

//CONSTANTES
#define METH_LOO 0
#define METH_PCA 1
#define METH_MEA 2

//INCLUDES
#include "../common/data.h"

class quantify_phenotype_data {
public:
	string chr, id;
	unsigned int start, end;

	quantify_phenotype_data(string _id, string _chr, unsigned int _start, unsigned int _end) {
		id = _id;
		chr = _chr;
		start = _start;
		end = _end;
	}

	bool operator < (const quantify_phenotype_data & d) const {
		if (chr < d.chr) return true;
		if (chr > d.chr) return false;
		if (start < d.start) return true;
		if (start > d.start) return false;
		if (end < d.end) return true;
		else return false;
	}
};

class qcluster {
public:
	string sid;
	int index;
	int count;
	int start, end;
	int first, last;
	qcluster * child1;
	qcluster * child2;

	qcluster (string _sid, int _index, int _count, qcluster * c1 = NULL, qcluster * c2 = NULL) {
		sid = _sid;
		index = _index;
		first = _index;
		last = _index;
		count = _count;
		start = 0;
		end = 0;
		child1 = c1;
		child2 = c2;
	};

	bool leaf() {
		return (child1 == NULL || child2 == NULL);
	}
};

class quantify_data : public data {
public:
	//REGIONS
	genomic_region regionPhenotype;
	set < string > modules;

	//PHENOTYPES STRUCT TO HANDLE MULTIPLE FILES
	unordered_map < string, unsigned int > phenotype_map;
	vector < quantify_phenotype_data > phenotype_vec;

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
	vector < bool > phenotype_neg;						//phenotype on the negative strand?

	//TREES
	vector < vector < float > > vec_quantifications;
	vector < qcluster * > vec_nodes;
	qcluster * tree;

	//GROUPS
	vector < vector < int > > vec_groups;

	//CONSTRUCTOR / DESTRUCTOR
	quantify_data() {};
    ~quantify_data() {};

	//READ OR GENERATE DATA
    void scanPhenotypes(string);
    void readPhenotypes(string);
    void sortPhenotypes();
    void readTree(string);
    void readModules(string);
    void readGroups(string);
    void normalTransform(vector < float > & V);

	//REGIONS
	bool setPhenotypeRegion(string reg) { return regionPhenotype.parse(reg); }

	//ANALYSIS
	void quantifyGroups(int method, int arg);
	void quantifyTrees(int method, int arg);
	void quantifyTrees(int method, int arg, qcluster * curr_node, vector < int > & curr_leaves);
	void quantify_by_pca(vector < int > & curr_leaves, vector < float > & Q, int ipc);
	void quantify_by_loo(vector < int > & curr_leaves, vector < float > & Q, bool variance);
	void quantify_by_mean(vector < int > & curr_leaves, vector < float > & Q);

	//OUTPUT
	void writePhenotypesTrees(string, bool);
	void writePhenotypesGroups(string, bool);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void quantify_main(vector < string > &);

#endif
