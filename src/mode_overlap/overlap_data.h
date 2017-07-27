#ifndef _OVERLAP_DATA_H
#define _OVERLAP_DATA_H

//INCLUDES
#include "../common/data.h"

class overlap_data : public data {
public:
	//REGIONS
	genomic_region regionPhenotype;

	//References
	vector < string > annotation_unique_type;			//annotation unique types
	vector < string > annotation_unique_chr;			//annotation unique chromosomes

	//Annotations
	int annotation_count;
	vector < int > annotation_chr;						//start position of the annotation
	vector < int > annotation_start;					//start position of the annotation
	vector < int > annotation_end;						//end position of the annotation
	vector < int > annotation_type;						//annotation type
	vector < IntervalTree < int > > annotation_tree;	//annotation tree

	//CONSTRUCTOR / DESTRUCTOR
	overlap_data() {};
    ~overlap_data() {};

	//READ OR GENERATE DATA
    void readAnnotations(string);

	//DATA MANAGEMENT
    void builtAnnotationtrees();

	//ANALYSIS
	void processingPhenotypes(string, string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void overlap_main(vector < string > &);

#endif
