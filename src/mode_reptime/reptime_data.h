#ifndef _REPTIME_DATA_H
#define _REPTIME_DATA_H

//INCLUDES
#include "../common/data.h"

class reptime_phenotype_data {
public:
	string chr, id;
	unsigned int start, end;

	reptime_phenotype_data(string _id, string _chr, unsigned int _start, unsigned int _end) {
		id = _id;
		chr = _chr;
		start = _start;
		end = _end;
	}

	bool operator < (const reptime_phenotype_data & d) const {
		if (chr < d.chr) return true;
		if (chr > d.chr) return false;
		if (start < d.start) return true;
		if (start > d.start) return false;
		if (end < d.end) return true;
		else return false;
	}
};

class reptime_data : public data {
public:
	//REGIONS
	genomic_region regionPhenotype;

	//PHENOTYPES
	int phenotype_count;										//phenotype number
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_chr;						//phenotype chr
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
	vector < reptime_phenotype_data > phenotype_vec;	//useful for merging

	//Replication times
	int reptime_count;
	vector < int > reptime_pos;							//position at which rep time has been measured
	vector < float > reptime_val;						//rep time value
	IntervalTree < float > reptime_tree;					//rep time tree representation

	//CONSTRUCTOR / DESTRUCTOR
	reptime_data() {};
    ~reptime_data() {};

	//READ OR GENERATE DATA
    void scanPhenotypes(string);
    void readPhenotypes(string);
    void sortPhenotypes();
    void readReptimes(string);

	//DATA MANAGEMENT
    void builtRTtree();

	//REGIONS
	bool setPhenotypeRegion(string);

	//ANALYSIS
	void project(string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void reptime_main(vector < string > &);

#endif
