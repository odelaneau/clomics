#include "trans_data.h"

bool trans_data::setPhenotypeRegion1(string reg) {
	return regionPhenotype1.parse(reg);
}

bool trans_data::setPhenotypeRegion2(string reg) {
	return regionPhenotype2.parse(reg);
}

void trans_data::allocate() {
	//corr_val = vector < vector < float > > (phenotype_count1, vector < float >(phenotype_count2, 0));
	hic_val = vector < vector < float > > (phenotype_count1, vector < float >(phenotype_count2, 0));
	hic_cnt = vector < vector < short > > (phenotype_count1, vector < short >(phenotype_count2, 0));
}

void trans_data::normalizePhenotypes1() {
	for (int x = 0 ; x < phenotype_count1 ; x++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += phenotype_val1[x][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			phenotype_val1[x][s] -= mean;
			sum += phenotype_val1[x][s] * phenotype_val1[x][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) phenotype_val1[x][s] /= sum;
	}
}

void trans_data::normalizePhenotypes2() {
	for (int x = 0 ; x < phenotype_count2 ; x++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += phenotype_val2[x][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			phenotype_val2[x][s] -= mean;
			sum += phenotype_val2[x][s] * phenotype_val2[x][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) phenotype_val2[x][s] /= sum;
	}
}
