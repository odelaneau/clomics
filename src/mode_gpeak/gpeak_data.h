#ifndef _GPEAK_DATA_H
#define _GPEAK_DATA_H

//INCLUDES
#include "../common/data.h"

class gcluster {
public:
	string sid;
	int index;
	gcluster * child1;
	gcluster * child2;

	gcluster (string _sid, int _index, gcluster * c1 = NULL, gcluster * c2 = NULL) {
		sid = _sid;
		index = _index;
		child1 = c1;
		child2 = c2;
	};

	bool leaf() {
		return (child1 == NULL || child2 == NULL);
	}
};

class gpeak_data : public data {
public:
	//MODULES
	set < string > modules;

	//TREE
	vector < gcluster * > vec_nodes;
	vector < string > vec_module_ids;
	gcluster * tree;

	//CONSTRUCTOR / DESTRUCTOR
	gpeak_data() {};
    ~gpeak_data() {};

	//READ OR GENERATE DATA
    void readTree(string);
    void readModules(string);

	//OUTPUT
	void process(string);
	void process(gcluster * , vector < int > & );
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void gpeak_main(vector < string > &);

#endif
