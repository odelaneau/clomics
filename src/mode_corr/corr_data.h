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

#ifndef _CORR_DATA_H
#define _CORR_DATA_H

//INCLUDES
#include "../common/data.h"

class corr_phenotype_data {
public:
	string chr, id;
	unsigned int start, end;

	corr_phenotype_data(string _id, string _chr, unsigned int _start, unsigned int _end) {
		id = _id;
		chr = _chr;
		start = _start;
		end = _end;
	}

	bool operator < (const corr_phenotype_data & d) const {
		if (chr < d.chr) return true;
		if (chr > d.chr) return false;
		if (start < d.start) return true;
		if (start > d.start) return false;
		if (end < d.end) return true;
		else return false;
	}
};

class corr_data : public data {
public:
	//REGIONS
	genomic_region regionPhenotype;

	//PHENOTYPES STRUCT TO HANDLE MULTIPLE FILES
	unordered_map < string, unsigned int > phenotype_map;
	vector < corr_phenotype_data > phenotype_vec;

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
	vector < bool > phenotype_neg;						//phenotype on the negative strand?
	vector < vector < int > > phenotype_paired;

	//CONSTRUCTOR / DESTRUCTOR
	corr_data() {};
    ~corr_data() {};

	//READ OR GENERATE DATA
    void scanPhenotypes(string);
    void readPhenotypes(string);
    void sortPhenotypes();

	//PHENOTYPE MANAGEMENT
	void imputePhenotypes();
	void normalizePhenotypes();
	void pairingPhenotypes(string);

	//REGIONS
	bool setPhenotypeRegion(string reg);

	//COMPUTATION METHODS [ALL INLINES FOR SPEED]
	double getCorrelation(vector < float > &, vector < float > &);
	double getPvalue(double corr, double df);

	//ANALYSIS
	void computeCorrelations(string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void corr_main(vector < string > &);

//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

/*
 * Fast implementation of inner_product optimized for 64 bytes cache lines.
 */
inline double corr_data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
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

inline double corr_data::getPvalue(double corr, double df) {
	double pval = pf(df * corr * corr / (1 - corr * corr), 1, df, 0, 0);
	if (pval <= std::numeric_limits<double>::min()) pval =std::numeric_limits<double>::min();
	return pval;
}


#endif
