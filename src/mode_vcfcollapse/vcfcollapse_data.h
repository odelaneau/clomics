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

#ifndef _VCFCOLLAPSE_DATA_H
#define _VCFCOLLAPSE_DATA_H

//INCLUDES
#include "../common/data.h"

class vcfcollapse_data : public data {
public :
	//PARAMETERS
	float param_max_maf;
	float param_min_maf;

	//GENOTYPES
	int genotype_count;									//variant site number
	vector < vector < float > > genotype_val;			//variant site genotype dosages
	vector < string > genotype_chr;						//variant site chromosome
	vector < string > genotype_id;						//variant site IDs
	vector < int > genotype_start;						//variant site start positions
	vector < int > genotype_end;						//variant site end positions
	vector < string > genotype_bed;

	//ANNOTATIONS
	int ann_count;										//Annotation number
	vector < int > ann_start;							//Annotation start position
	vector < int > ann_end;								//Annotation end position
	vector < string > ann_chr;							//Annotation chromosome
	vector < string > ann_id;							//Annotation chromosome
	map  < string, int > ann_chr_map;
	vector <  IntervalTree < int > > ann_int;			//Annotation intervals

	//RESULTS
	int burden_count;									//variant site number
	vector < vector < float > > burden_val;				//variant site genotype dosages
	vector < string > burden_chr;						//variant site chromosome
	vector < string > burden_id;						//variant site IDs
	vector < int > burden_start;						//variant site start positions
	vector < int > burden_end;							//variant site end positions

	//CONSTRUCTOR/DESTRUCTOR
	vcfcollapse_data() {
	}

	~vcfcollapse_data() {
		genotype_count = 0;
		genotype_val.clear();
		genotype_chr.clear();
		genotype_id.clear();
		genotype_start.clear();
		genotype_end.clear();
	}

	//
	void normalTransform(vector < float > & );
	void readGenotypes(string);
	void readAnnotations(string);
	void buildTrees();
	void collapseGenotypes(bool);
	void writeGenotypes(string);
};

void vcfcollapse_main(vector < string > & );

#endif
