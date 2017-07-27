#include "trans_data.h"

void trans_data::readHiCnorms1(string fnorm) {
	int line = 0, start = 0, end = 0, n_kept = 0, n_remo = 0;
	string buffer;

	vrb.title("Reading Hi-C normalization data for chunk 1 in [" + fnorm + "]");
	input_file fd (fnorm);
	while (getline(fd, buffer)) {
		start = line * hic_bin1;
		end = (line + 1) * hic_bin1 - 1;
		vector < Interval < int > > phenotype_overlap;
		phenotype_tree1.findOverlapping(start, end, phenotype_overlap);
		hic2pheno1.push_back(vector < int >());
		hic_norm1.push_back(-1.0);
		if (phenotype_overlap.size() > 0 && buffer != "NaN") {
			hic_norm1.back() = atof(buffer.c_str());
			for (int o = 0 ; o < phenotype_overlap.size() ; o ++) hic2pheno1.back().push_back(phenotype_overlap[o].value);
			n_kept ++;
		} else n_remo ++;
		line ++;
	}
	vrb.bullet(stb.str(n_kept) + " Hi-C bins retained");
	vrb.bullet(stb.str(n_remo) + " Hi-C bins removed");
}

void trans_data::readHiCnorms2(string fnorm) {
	int line = 0, start = 0, end = 0, n_kept = 0, n_remo = 0;
	string buffer;

	vrb.title("Reading Hi-C normalization data for chunk 2 in [" + fnorm + "]");
	input_file fd (fnorm);
	while (getline(fd, buffer)) {
		start = line * hic_bin2;
		end = (line + 1) * hic_bin2 - 1;
		vector < Interval < int > > phenotype_overlap;
		phenotype_tree2.findOverlapping(start, end, phenotype_overlap);
		hic2pheno2.push_back(vector < int >());
		hic_norm2.push_back(-1.0);
		if (phenotype_overlap.size() > 0 && buffer != "NaN") {
			hic_norm2.back() = atof(buffer.c_str());
			for (int o = 0 ; o < phenotype_overlap.size() ; o ++) hic2pheno2.back().push_back(phenotype_overlap[o].value);
			n_kept ++;
		} else n_remo ++;
		line ++;
	}
	vrb.bullet(stb.str(n_kept) + " Hi-C bins retained");
	vrb.bullet(stb.str(n_remo) + " Hi-C bins removed");
}

void trans_data::readHiCvalues(string fvalues) {
	long n_set = 0, n_unset = 0, n_update = 0;
	string buffer;
	vector < string > tokens;

	vrb.title("Reading Hi-C interaction data in [" + fvalues + "]");
	input_file fd (fvalues);
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		assert(tokens.size() == 3);

		int idx1 = atoi(tokens[0].c_str()) / hic_bin1;
		int idx2 = atoi(tokens[1].c_str()) / hic_bin2;

		if (hic_norm1[idx1] >= 0.0 && hic_norm2[idx2] >= 0.0) {
			for (int p1 = 0 ; p1 < hic2pheno1[idx1].size() ; p1 ++) for (int p2 = 0 ; p2 < hic2pheno2[idx2].size() ; p2 ++) {
				hic_val[hic2pheno1[idx1][p1]][hic2pheno2[idx2][p2]] += atof(tokens[2].c_str()) / (hic_norm1[idx1] * hic_norm2[idx2]);
				hic_cnt[hic2pheno1[idx1][p1]][hic2pheno2[idx2][p2]] ++;
				n_update ++;
			}
			n_set ++;
		} else n_unset++;
	}
	vrb.bullet(stb.str(n_update) + " matrix updates");
	vrb.bullet(stb.str(n_set) + " lines used");
	vrb.bullet(stb.str(n_unset) + " lines unused");
}

