#ifndef _TRANS_DATA_H
#define _TRANS_DATA_H

//INCLUDES
#include "../common/data.h"

class trans_phenotype_data {
public:
	string chr, id;
	unsigned int start, end;

	trans_phenotype_data(string _id, string _chr, unsigned int _start, unsigned int _end) {
		id = _id;
		chr = _chr;
		start = _start;
		end = _end;
	}

	bool operator < (const trans_phenotype_data & d) const {
		if (chr < d.chr) return true;
		if (chr > d.chr) return false;
		if (start < d.start) return true;
		if (start > d.start) return false;
		if (end < d.end) return true;
		else return false;
	}
};

class trans_data : public data {
public:
	//PHENOTYPES CHUNK1
	int phenotype_min_pos1;
	int phenotype_max_pos1;
	genomic_region regionPhenotype1;
	int phenotype_count1;
	vector < string > phenotype_chr1;
	vector < string > phenotype_id1;
	vector < int > phenotype_start1;
	vector < double > phenotype_rpos1;
	vector < int > phenotype_end1;
	IntervalTree < int > phenotype_tree1;
	vector < vector < float > > phenotype_val1;
	vector < trans_phenotype_data > phenotype_vec1;
	unordered_map < string, unsigned int > phenotype_map1;

	//PHENOTYPES CHUNK2
	int phenotype_min_pos2;
	int phenotype_max_pos2;
	genomic_region regionPhenotype2;
	int phenotype_count2;
	vector < string > phenotype_chr2;
	vector < string > phenotype_id2;
	vector < int > phenotype_start2;
	vector < int > phenotype_end2;
	vector < double > phenotype_rpos2;
	IntervalTree < int > phenotype_tree2;
	vector < vector < float > > phenotype_val2;
	vector < trans_phenotype_data > phenotype_vec2;
	unordered_map < string, unsigned int > phenotype_map2;

	//HI-C CHUNK1
	int hic_bin1;
	vector < float > hic_norm1;
	vector < vector < int > > hic2pheno1;

	//HI-C CHUNK2
	int hic_bin2;
	vector < float > hic_norm2;
	vector < vector < int > > hic2pheno2;

	//SQUARED DATA
	double thinin;
	//vector < vector < float > > corr_val;
	vector < vector < float > > hic_val;
	vector < vector < short > > hic_cnt;

	//OTHER DATA
	int centromere_position;

	//CONSTRUCTOR / DESTRUCTOR
	trans_data() {};
    ~trans_data() {};

	//READ OR GENERATE DATA
    void scanPhenotypes1(string);
    void sortPhenotypes1();
    void readPhenotypes1(string);
    void scanPhenotypes2(string);
    void sortPhenotypes2();
	void readPhenotypes2(string);
    void readHiCnorms1(string);
    void readHiCnorms2(string);

	//REGIONS
	bool setPhenotypeRegion1(string);
	bool setPhenotypeRegion2(string);
	void allocate();

	//COMPUTATION METHODS [ALL INLINES FOR SPEED]
	void normalizePhenotypes1();
	void normalizePhenotypes2();
	double getCorrelation(vector < float > &, vector < float > &);
	double getPvalue(double corr, double df);
	void computeRelativePositions(int, int);

	//ANALYSIS
	void readHiCvalues(string);
	void computeHicInteractions(string);
	void computeAllInteractions(string);
	void computeSigInteractions(string, double);
	void computeAllInteractionsWithDesc(string fout);
	void binRelativePositions(int n_bins, string fout);
	void computeBinnedInteractions(string fout, double bin_size);

};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void trans_main(vector < string > &);

//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//
inline double trans_data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
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

inline double trans_data::getPvalue(double corr, double df) {
	double pval = pf(df * corr * corr / (1 - corr * corr), 1, df, 0, 0);
	if (pval <= std::numeric_limits<double>::min()) pval =std::numeric_limits<double>::min();
	return pval;
}

#endif
