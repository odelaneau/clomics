#include "span_data.h"

void span_data::normalizePhenotypes() {
	for (int x = 0 ; x < phenotype_count ; x++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += phenotype_val[x][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			phenotype_val[x][s] -= mean;
			sum += phenotype_val[x][s] * phenotype_val[x][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) phenotype_val[x][s] /= sum;
	}
}

bool span_data::setPhenotypeRegion(string reg) {
	return regionPhenotype.parse(reg);
}
