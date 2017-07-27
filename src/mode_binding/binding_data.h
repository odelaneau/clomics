#ifndef _BINDING_DATA_H
#define _BINDING_DATA_H

//INCLUDES
#include "../common/data.h"

class binding_data : public data {
public:
	//REFERENCE SEQUENCES
	map < string, string > refseq;

	//VARIANTS
	vector < string > chr;
	vector < int > pos;
	vector < string > ref;
	vector < string > alt;

	//CONSTRUCTOR / DESTRUCTOR
	binding_data() {};
    ~binding_data() {};

    //READ OR GENERATE DATA
    void readReference(string);
    void readPosition(string);

    //PROCESS
    void process();

    //WRITE OUTPUT
    void writeFasta(string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void binding_main(vector < string > &);

#endif
