#include "motif_trans_data.h"

void motif_trans_data::computeCorrelations() {
	vrb.title("Compute all possible correlations");
	vrb.bullet("n=" + stb.str(phenotype_count1 * phenotype_count2));

	for (int p1 = 0 ; p1 < phenotype_count1 ; p1 ++) {
		for (int p2 = 0 ; p2 < phenotype_count2 ; p2 ++) {
			double curr_corr = abs(getCorrelation(phenotype_val1[p1], phenotype_val2[p2]));
			for (int m1 = 0 ; m1 < phenotype_motives1[p1].size() ; m1 ++) {
				for (int m2 = 0 ; m2 < phenotype_motives2[p2].size() ; m2 ++) {
					int idx0 = phenotype_motives1[p1][m1];
					int idx1 = phenotype_motives2[p2][m2];
					if (idx1 < idx0) {
						int tmp = idx1;
						idx1 = idx0;
						idx0 = tmp;
					}
					interaction_cor[idx0][idx1] += curr_corr;
					interaction_cnt[idx0][idx1] ++;
				}
			}
		}
		vrb.progress((p1+1) * 1.0 / (phenotype_count1 - 1));
	}
}
