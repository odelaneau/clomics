#ifndef _TRANS_FILTER_DATA_H
#define _TRANS_FILTER_DATA_H

//INCLUDES
#include "../common/data.h"

class trans_filter_data : public data {
public:
	//INTERACTIONS
	vector < pair < string, string > > inter_chr;
	vector < pair < int, int > > inter_start;
	vector < pair < int, int > > inter_end;

	//DATA
	int bin;
	double pva;
	double hic;

	//CONSTRUCTOR / DESTRUCTOR
	trans_filter_data() {};
    ~trans_filter_data() {};

	//READ OR GENERATE DATA
    void readInteractions(string);
    void writeInteractions(string file);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void trans_filter_main(vector < string > &);

#endif
