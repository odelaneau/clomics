#ifndef _SPAN_DATA_H
#define _SPAN_DATA_H

//INCLUDES
#include "../common/data.h"

class span_phenotype_data {
public:
	string chr, id;
	unsigned int start, end;

	span_phenotype_data(string _id, string _chr, unsigned int _start, unsigned int _end) {
		id = _id;
		chr = _chr;
		start = _start;
		end = _end;
	}

	bool operator < (const span_phenotype_data & d) const {
		if (chr < d.chr) return true;
		if (chr > d.chr) return false;
		if (start < d.start) return true;
		if (start > d.start) return false;
		if (end < d.end) return true;
		else return false;
	}
};

class span_data : public data {
public:
	//REGIONS
	genomic_region regionPhenotype;

	//PHENOTYPES
	int phenotype_count;
	unordered_map < string, unsigned int > phenotype_map;
	vector < span_phenotype_data > phenotype_vec;

	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions

	//ANNOTATIONS
	IntervalTree < int > annotation_tree;
	map < string, int > annotation_id_map;
	vector < string > annotation_id_vec;

	//CONSTRUCTOR / DESTRUCTOR
	span_data() {};
    ~span_data() {};

	//READ OR GENERATE DATA
    void scanPhenotypes(string);
    void sortPhenotypes();
    void readPhenotypes(string);
    void readAnnotations(string);

	//DATA MANAGEMENT
    void normalizePhenotypes();
    double getCorrelation(vector < float > & vec1, vector < float > & vec2);
    double getPvalue(double corr, double df);

	//REGIONS
	bool setPhenotypeRegion(string);

	//ANALYSIS
	void computeInteractionsSingle(string, int, bool);
	void computeInteractionsPaired(string, int, bool);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void span_main(vector < string > &);

inline double span_data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
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

inline double span_data::getPvalue(double corr, double df) {
	double pval = pf(df * corr * corr / (1 - corr * corr), 1, df, 0, 0);
	if (pval <= std::numeric_limits<double>::min()) pval = std::numeric_limits<double>::min();
	return pval;
}


#endif
