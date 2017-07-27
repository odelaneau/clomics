#include "hic_data.h"

void hic_data::readHiCnorms(string fnorm) {
	int line = 0, start = 0, end = 0, n_kept = 0, n_remo = 0;
	string buffer;

	vrb.title("Reading Hi-C normalization data in [" + fnorm + "]");
	input_file fd (fnorm);
	while (getline(fd, buffer)) {

		//
		start = line * hic_bin;
		end = (line + 1) * hic_bin - 1;

		//
		vector < Interval < int > > phenotype_overlap;
		phenotype_tree.findOverlapping(start, end, phenotype_overlap);

		//
		hic2pheno.push_back(vector < int >());
		hic_norm.push_back(-1.0);
		if (phenotype_overlap.size() > 0 && buffer != "NaN") {
			hic_norm.back() = atof(buffer.c_str());
			for (int o = 0 ; o < phenotype_overlap.size() ; o ++) hic2pheno.back().push_back(phenotype_overlap[o].value);
			n_kept ++;
		} else n_remo ++;

		//
		line ++;
	}
	vrb.bullet(stb.str(n_kept) + " Hi-C bins retained");
	vrb.bullet(stb.str(n_remo) + " Hi-C bins removed");

}

void hic_data::readHiCvalues(string fvalues) {
	long n_set = 0, n_unset = 0, n_update = 0;
	string buffer;
	vector < string > tokens;

	vrb.title("Reading Hi-C interaction data in [" + fvalues + "]");
	input_file fd (fvalues);
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		assert(tokens.size() == 3);

		//
		int idx0 = atoi(tokens[0].c_str()) / hic_bin;
		int idx1 = atoi(tokens[1].c_str()) / hic_bin;

		//
		if (hic_norm[idx0] >= 0.0 && hic_norm[idx1] >= 0.0) {
			for (int p0 = 0 ; p0 < hic2pheno[idx0].size() ; p0 ++) for (int p1 = 0 ; p1 < hic2pheno[idx1].size() ; p1 ++) {
				phenotype_hic_val[hic2pheno[idx0][p0]][hic2pheno[idx1][p1]] += atof(tokens[2].c_str()) / (hic_norm[idx0] * hic_norm[idx1]);
				phenotype_hic_cnt[hic2pheno[idx0][p0]][hic2pheno[idx1][p1]] ++;
				phenotype_hic_val[hic2pheno[idx1][p1]][hic2pheno[idx0][p0]] += atof(tokens[2].c_str()) / (hic_norm[idx0] * hic_norm[idx1]);
				phenotype_hic_cnt[hic2pheno[idx1][p1]][hic2pheno[idx0][p0]] ++;
				n_update ++;
			}
			n_set ++;
		} else n_unset++;
	}
	vrb.bullet(stb.str(n_update) + " matrix updates");
	vrb.bullet(stb.str(n_set) + " lines used");
	vrb.bullet(stb.str(n_unset) + " lines unused");
}

