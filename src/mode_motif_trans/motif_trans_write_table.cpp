#include "motif_trans_data.h"

void motif_trans_data::writePairwiseData(string fout) {
	vrb.title("Write output in [" + fout + "]");
	output_file fd(fout);
	fd << "TF1 TF2 MEAN" << endl;
	for (int m0=0 ; m0 < motif_trans_unique_id.size() ; m0 ++) {
		for (int m1=m0 ; m1 < motif_trans_unique_id.size() ; m1 ++) {
			fd << motif_trans_unique_id[m0] << " " << motif_trans_unique_id[m1] << " " << interaction_cor[m0][m1]  << " " << interaction_cnt[m0][m1] << endl;
		}
	}
	fd.close();
}

