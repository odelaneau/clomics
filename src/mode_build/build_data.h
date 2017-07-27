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

#ifndef _BUILD_DATA_H
#define _BUILD_DATA_H

//INCLUDES
#include "../common/data.h"

class build_phenotype_data {
public:
	string chr, id;
	unsigned int start, end;

	build_phenotype_data(string _id, string _chr, unsigned int _start, unsigned int _end) {
		id = _id;
		chr = _chr;
		start = _start;
		end = _end;
	}

	bool operator < (const build_phenotype_data & d) const {
		if (chr < d.chr) return true;
		if (chr > d.chr) return false;
		if (start < d.start) return true;
		if (start > d.start) return false;
		if (end < d.end) return true;
		else return false;
	}
};

class bcluster {
public:
	string sid;
	int index;
	int count;
	double acorrelation, apvalue;
	double bcorrelation, bpvalue;
	double ccorrelation, cpvalue;
	double dcorrelation, dpvalue;
	double ecorrelation, epvalue;
	double pcavariance;
	bcluster * child1;
	bcluster * child2;

	bcluster (string _sid, int _index, int _count, bcluster * c1 = NULL, bcluster * c2 = NULL) {
		sid = _sid;
		index = _index;
		count = _count;
		acorrelation = 0.0;
		apvalue = 1.0;
		bcorrelation = 0.0;
		bpvalue = 1.0;
		ccorrelation = 0.0;
		cpvalue = 1.0;
		dcorrelation = 0.0;
		dpvalue = 1.0;
		ecorrelation = 0.0;
		epvalue = 1.0;
		pcavariance = 0.0;
		child1 = c1;
		child2 = c2;
	};

	void vectorize(vector < bcluster * > & nodes) {
		if (child1 != NULL) child1->vectorize(nodes);
		if (child2 != NULL) child2->vectorize(nodes);
		if (nodes.size() < index + 1) nodes.resize(index + 1, NULL);
		nodes[index] = this;
	}
};

class build_data : public data {
public:
	//MODE
	bool resampling;

	//REGIONS
	genomic_region regionPhenotype;

	//PHENOTYPES STRUCT TO HANDLE MULTIPLE FILES
	unordered_map < string, unsigned int > phenotype_map;
	vector < build_phenotype_data > phenotype_vec;

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
	vector < bool > phenotype_neg;						//phenotype on the negative strand?

	//TREE
	bcluster * tree;

	//CONSTRUCTOR / DESTRUCTOR
	build_data() { resampling = false; };
    ~build_data() {};

	//READ OR GENERATE DATA
    void scanPhenotypes(string);
    void readPhenotypes(string);
    void sortPhenotypes();
    void bootstrapPhenotypes(float);
    void jackknifePhenotypes(float);

	//PHENOTYPE MANAGEMENT
	void imputePhenotypes();
	void normalizePhenotypes();

	//REGIONS
	bool setPhenotypeRegion(string reg);

	//COMPUTATION METHODS [ALL INLINES FOR SPEED]
	double getCorrelation(vector < float > &, vector < float > &);
	double getPvalue(double, double);

	//ANALYSIS
	void clusterize();

	//OUTPUT
	void writeTree(string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void build_main(vector < string > &);

//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

/*
 * Fast implementation of inner_product optimized for 64 bytes cache lines.
 */
inline double build_data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
	int i = 0;
	int repeat = (sample_count / 4);
	int left = (sample_count % 4);
	double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

	while (repeat --) {
		sum0 += vec1[i] * vec2[i];
		sum1 += vec1[i+1] * vec2[i+1];
		sum2 += vec1[i+2] * vec2[i+2];
		sum3 += vec1[i+3] * vec2[i+3];
		i += 4;
	}

	switch (left) {
	case 3:	sum0 += vec1[i+2] * vec2[i+2];
	case 2:	sum0 += vec1[i+1] * vec2[i+1];
	case 1:	sum0 += vec1[i+0] * vec2[i+0];
	case 0: ;
	}

	return sum0 + sum1 + sum2 + sum3;
}

inline double build_data::getPvalue(double corr, double df) {
	return pf(df * corr * corr / (1 - corr * corr), 1, df, 0, 0);
}

#endif
