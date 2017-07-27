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

#ifndef _MOTIF_TRANS_DATA_H
#define _MOTIF_TRANS_DATA_H

//INCLUDES
#include "../common/data.h"

class motif_trans_phenotype_data {
public:
	string chr, id;
	unsigned int start, end;

	motif_trans_phenotype_data(string _id, string _chr, unsigned int _start, unsigned int _end) {
		id = _id;
		chr = _chr;
		start = _start;
		end = _end;
	}

	bool operator < (const motif_trans_phenotype_data & d) const {
		if (chr < d.chr) return true;
		if (chr > d.chr) return false;
		if (start < d.start) return true;
		if (start > d.start) return false;
		if (end < d.end) return true;
		else return false;
	}
};

class motif_trans_data : public data {
public:
	//REGIONS
	genomic_region regionPhenotype1;
	genomic_region regionPhenotype2;

	//PHENOTYPES STRUCT TO HANDLE MULTIPLE FILES
	unordered_map < string, unsigned int > phenotype_map1;
	vector < motif_trans_phenotype_data > phenotype_vec1;
	unordered_map < string, unsigned int > phenotype_map2;
	vector < motif_trans_phenotype_data > phenotype_vec2;

	//PHENOTYPES 1
	int phenotype_count1;								//phenotype number
	vector < vector < float > > phenotype_val1;			//phenotype values
	vector < string > phenotype_id1;					//phenotype ids
	vector < string > phenotype_chr1;					//phenotype chromosomes
	vector < int > phenotype_start1;					//phenotype start positions
	vector < int > phenotype_end1;						//phenotype end positions
	vector < bool > phenotype_neg1;						//phenotype on the negative strand?
	IntervalTree < int > phenotype_tree1;
	vector < vector < int > > phenotype_motives1;		//phenotype to motives mapping

	//PHENOTYPES 2
	int phenotype_count2;								//phenotype number
	vector < vector < float > > phenotype_val2;			//phenotype values
	vector < string > phenotype_id2;					//phenotype ids
	vector < string > phenotype_chr2;					//phenotype chromosomes
	vector < int > phenotype_start2;					//phenotype start positions
	vector < int > phenotype_end2;						//phenotype end positions
	vector < bool > phenotype_neg2;						//phenotype on the negative strand?
	IntervalTree < int > phenotype_tree2;
	vector < vector < int > > phenotype_motives2;		//phenotype to motives mapping

	//MOTIVES
	vector < string > motif_trans_id1;
	vector < int > motif_trans_start1;
	vector < int > motif_trans_stop1;
	vector < int > motif_trans_idx1;
	vector < string > motif_trans_id2;
	vector < int > motif_trans_start2;
	vector < int > motif_trans_stop2;
	vector < int > motif_trans_idx2;
	map < string, int > motif_trans_map;
	vector < string > motif_trans_unique_id;

	//INTERACTIONS
	vector < vector < double > > interaction_cor;
	vector < vector < int > > interaction_cnt;

	//CONSTRUCTOR / DESTRUCTOR
	motif_trans_data() {};
    ~motif_trans_data() {};

	//READ OR GENERATE DATA
    void scanPhenotypes1(string);
    void readPhenotypes1(string);
    void sortPhenotypes1();
    void scanPhenotypes2(string);
    void readPhenotypes2(string);
    void sortPhenotypes2();

	//PHENOTYPE MANAGEMENT
	void imputePhenotypes1();
	void normalizePhenotypes1();
	void imputePhenotypes2();
	void normalizePhenotypes2();

	//MOTIVES
	void readMotives1(string, bool);
	void readMotives2(string, bool);
	void initializeMotives();
	void mappingMotives1();
	void mappingMotives2();

	//REGIONS
	bool setPhenotypeRegion1(string reg);
	bool setPhenotypeRegion2(string reg);

	//COMPUTATION METHODS
	double getCorrelation(vector < float > &, vector < float > &);
	double getCorrelationThreshold(double);

	//ANALYSIS
	void computeCorrelations();
	void buildIntervalTree1();
	void buildIntervalTree2();

	//WRITE
	void writePairwiseData(string fout);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void motif_trans_main(vector < string > &);

//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

/*
 * Fast implementation of inner_product optimized for 64 bytes cache lines.
 */
inline double motif_trans_data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
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

inline double motif_trans_data::getCorrelationThreshold(double pvalue) {
	double p = qf(pvalue, 1, sample_count - 2, 0, 0);
	return sqrt(p / (sample_count - 2 + p));
}

#endif

