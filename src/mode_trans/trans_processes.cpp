#include "trans_data.h"

void trans_data::computeBinnedInteractions(string fout, double bin_size) {
	vrb.title("Report interactions within BINS in [" + fout + "]");

	int n_I = (int)ceil(phenotype_end1.back() / bin_size) + 1;
	int n_J = (int)ceil(phenotype_end2.back() / bin_size) + 1;

	vector < vector < double > > CORR = vector < vector < double > > (n_I, vector < double > (n_J, 0.0));
	vector < vector < int > > COUNT = vector < vector < int > > (n_I, vector < int > (n_J, 0));

	for (int i = 0 ; i < phenotype_count1 ; i ++) {
		for (int j = 0 ; j < phenotype_count2 ; j ++) {
			int idx_I = (int)round((phenotype_start1[i] + phenotype_end1[i])/(2 * bin_size));
			int idx_J = (int)round((phenotype_start2[j] + phenotype_end2[j])/(2 * bin_size));
			double corr = abs(getCorrelation(phenotype_val1[i], phenotype_val2[j]));
			if (CORR[idx_I][idx_J] < corr) CORR[idx_I][idx_J] = corr;
			COUNT[idx_I][idx_J] = 1;
		}
	}

	output_file fd(fout.c_str());
	for (int I = 0; I < n_I ; I ++) for (int J = 0; J < n_J ; J ++) fd << I << " " << I * bin_size << " " << J << " " << J * bin_size << " " << COUNT[I][J] << " " << CORR[I][J] << endl;
	fd.close();
}


void trans_data::computeHicInteractions(string fout) {
	vrb.title("Report all possible interactions with Hi-C in [" + fout + "]");
	output_file fd(fout.c_str());
	for (int i = 0 ; i < phenotype_count1 ; i ++) {
		for (int j = 0 ; j < phenotype_count2 ; j ++) {
			if (hic_cnt[i][j] > 0) {
				double corr = getCorrelation(phenotype_val1[i], phenotype_val2[j]);
				fd << i + 1 << " " << phenotype_chr1[i] << " " << phenotype_start1[i] << " " << phenotype_end1[i] << " " << phenotype_id1[i];
				if (phenotype_rpos1.size() > 0) fd << " " << phenotype_rpos1[i];
				fd << " " << j + 1 << " " << phenotype_chr2[j] << " " << phenotype_start2[j] << " " << phenotype_end2[j] << " " << phenotype_id2[j];
				if (phenotype_rpos2.size() > 0) fd << " " << phenotype_rpos2[j];
				fd << " " << stb.str(hic_val[i][j]/hic_cnt[i][j], 3) << " " << stb.str(corr, 3)  << " " << getPvalue(corr, sample_count-2) << endl;
			} else if (rng.getDouble() <= thinin) {
				double corr = getCorrelation(phenotype_val1[i], phenotype_val2[j]);
				fd << i + 1 << " " << phenotype_chr1[i] << " " << phenotype_start1[i] << " " << phenotype_end1[i] << " " << phenotype_id1[i];
				if (phenotype_rpos1.size() > 0) fd << " " << phenotype_rpos1[i];
				fd << " " << j + 1 << " " << phenotype_chr2[j] << " " << phenotype_start2[j] << " " << phenotype_end2[j] << " " << phenotype_id2[j];
				if (phenotype_rpos2.size() > 0) fd << " " << phenotype_rpos2[j];
				fd << " 0 " << stb.str(corr, 3)  << " " << getPvalue(corr, sample_count-2) << endl;
			}
		}
		vrb.progress((i+1) * 1.0 / phenotype_count1);
	}
	fd.close();
}

void trans_data::computeAllInteractions(string fout) {
	vrb.title("Report all possible interactions in [" + fout + "]");
	output_file fd(fout.c_str());
	for (int i = 0 ; i < phenotype_count1 ; i ++) {
		for (int j = 0 ; j < phenotype_count2 ; j ++) {
			double corr = getCorrelation(phenotype_val1[i], phenotype_val2[j]);
			fd << i + 1 << " " << j + 1 <<  " " << stb.str(corr, 3) << " " << getPvalue(corr, sample_count-2) << endl;
		}
		vrb.progress((i+1) * 1.0 / phenotype_count1);
	}
	fd.close();
}

void trans_data::computeSigInteractions(string fout, double threshold) {
	vrb.title("Report significant interactions in [" + fout + "] with pvalue<" + stb.str(threshold));
	output_file fd(fout.c_str());
	for (int i = 0 ; i < phenotype_count1 ; i ++) {
		for (int j = 0 ; j < phenotype_count2 ; j ++) {
			double corr = getCorrelation(phenotype_val1[i], phenotype_val2[j]);
			double pval = getPvalue(corr, sample_count-2);
			if (pval < threshold) {
				fd << i + 1 << " " << (phenotype_start1[i]+phenotype_end1[i])/2;
				fd << " " << j + 1 << " " << (phenotype_start2[j]+phenotype_end2[j])/2;
				fd << " " << corr << " " << pval;
				if (phenotype_rpos1.size() == 0) fd << endl;
				else fd << " " << phenotype_rpos1[i] << " " << phenotype_rpos2[j] << endl;
			}
		}
		vrb.progress((i+1) * 1.0 / phenotype_count1);
	}
	fd.close();
}


void trans_data::computeAllInteractionsWithDesc(string fout) {
	vrb.title("Report all possible interactions with full output in [" + fout + "]");
	output_file fd(fout.c_str());
	for (int i = 0 ; i < phenotype_count1 ; i ++) {
		for (int j = 0 ; j < phenotype_count2 ; j ++) {
			double corr = getCorrelation(phenotype_val1[i], phenotype_val2[j]);
			fd << i + 1 << " " << phenotype_chr1[i] <<  	" " << (phenotype_start1[i]+phenotype_end1[i])/2 << " " << phenotype_id1[i];;
			fd << " " << j + 1 << " " << phenotype_chr2[j] << " " << (phenotype_start2[j]+phenotype_end2[j])/2 << " " << phenotype_id2[j];;
			fd <<  " " << stb.str(corr, 3)  << " " << getPvalue(corr, sample_count-2) << endl;
		}
		vrb.progress((i+1) * 1.0 / phenotype_count1);
	}
	fd.close();
}

