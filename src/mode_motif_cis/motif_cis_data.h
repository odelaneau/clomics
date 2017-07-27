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

#ifndef _MOTIF_CIS_DATA_H
#define _MOTIF_CIS_DATA_H

//INCLUDES
#include "../common/data.h"

class motif_cis_phenotype_data {
public:
	string chr, id;
	unsigned int start, end;

	motif_cis_phenotype_data(string _id, string _chr, unsigned int _start, unsigned int _end) {
		id = _id;
		chr = _chr;
		start = _start;
		end = _end;
	}

	bool operator < (const motif_cis_phenotype_data & d) const {
		if (chr < d.chr) return true;
		if (chr > d.chr) return false;
		if (start < d.start) return true;
		if (start > d.start) return false;
		if (end < d.end) return true;
		else return false;
	}
};

class motif_cis_data : public data {
public:
	//REGIONS
	genomic_region regionPhenotype;

	//PHENOTYPES STRUCT TO HANDLE MULTIPLE FILES
	unordered_map < string, unsigned int > phenotype_map;
	vector < motif_cis_phenotype_data > phenotype_vec;

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
	vector < bool > phenotype_neg;						//phenotype on the negative strand?
	IntervalTree < int > phenotype_tree;
	vector < vector < int > > phenotype_motives;		//phenotype to motives mapping

	//MOTIVES
	vector < string > motif_cis_id;
	vector < string > motif_cis_unique_id;
	vector < int > motif_cis_start;
	vector < int > motif_cis_stop;
	vector < int > motif_cis_idx;
	map < string, int > motif_cis_map;

	//INTERACTIONS
	vector < vector < double > > interaction_cor;
	vector < vector < int > > interaction_cnt;

	//CONSTRUCTOR / DESTRUCTOR
	motif_cis_data() {};
    ~motif_cis_data() {};

	//READ OR GENERATE DATA
    void scanPhenotypes(string);
    void readPhenotypes(string);
    void sortPhenotypes();

	//PHENOTYPE MANAGEMENT
	void imputePhenotypes();
	void normalizePhenotypes();

	//MOTIVES
	void readMotives(string, bool);
	void initializeMotives();
	void mappingMotives();

	//REGIONS
	bool setPhenotypeRegion(string reg);

	//COMPUTATION METHODS
	double getCorrelation(vector < float > &, vector < float > &);
	double getCorrelationThreshold(double);

	//ANALYSIS
	void computeCorrelations(int);
	void buildIntervalTree();

	//WRITE
	void writePairwiseData(string fout);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void motif_cis_main(vector < string > &);

//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

/*
 * Fast implementation of inner_product optimized for 64 bytes cache lines.
 */
inline double motif_cis_data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
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

inline double motif_cis_data::getCorrelationThreshold(double pvalue) {
	double p = qf(pvalue, 1, sample_count - 2, 0, 0);
	return sqrt(p / (sample_count - 2 + p));
}

#endif

