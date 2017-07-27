#include "overlap_data.h"

void overlap_data::builtAnnotationtrees() {
	vrb.title("Building interval trees for replication timings");
	vector < Interval < int > > interval_vec;
	annotation_tree = vector < IntervalTree < int > > (annotation_unique_chr.size());
	for (int c = 0 ; c < annotation_unique_chr.size() ; c++) {
		interval_vec.clear();
		for (int a = 0 ; a < annotation_count ; a ++) if (annotation_chr[a] == c) interval_vec.push_back(Interval < int > (annotation_start[a], annotation_end[a], a));
		annotation_tree[c] = IntervalTree < int > (interval_vec);
		vrb.bullet("Chromosome [" + annotation_unique_chr[c] + "] filled with " + stb.str(interval_vec.size()) + " intervals");
	}
}

void overlap_data::processingPhenotypes(string fbed, string fout) {
	string buffer;
	vector < string > tokens;

	//Open BED file
	int n_overlaps = 0, n_total = 0;
	vrb.title("Processing phenotype data in [" + fbed + "] and writing overlaps in [" + fout + "]");
	input_file fdi (fbed);
	output_file fdo (fout);
	fdo << "chr\tstart\tend\tid";
	for (int a = 0 ; a < annotation_unique_type.size() ; a ++) fdo << "\t" << annotation_unique_type[a];
	fdo << endl;
	while (getline(fdi, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() < 4) vrb.error("Incorrect number of column somewhere in the file!");

		int idx_chr = -1;
		for (int e = 0 ; e < annotation_unique_chr.size() && idx_chr < 0 ; e++) if (annotation_unique_chr[e] == tokens[0]) idx_chr = e;
		if (idx_chr >= 0) {
			int start = atoi(tokens[1].c_str());
			int end = atoi(tokens[2].c_str());
			vector < bool > mask = vector < bool > (annotation_unique_type.size(), false);
			vector < Interval < int > > overlaps;
			annotation_tree[idx_chr].findOverlapping(start, end, overlaps);
			for (int o = 0 ; o < overlaps.size() ; o ++) mask[annotation_type[overlaps[o].value]] = true;
			fdo << tokens[0] << "\t" << tokens[1] << "\t" << tokens[2] << "\t" << tokens[3];
			for (int m = 0 ; m < mask.size() ; m ++) fdo << "\t" << mask[m];
			fdo << endl;
			n_overlaps++;
		}
		n_total++;
	}
	vrb.bullet(stb.str(n_overlaps) + " / " + stb.str(n_total) + " phenotypes overlapped by at least one annotation");
	fdo.close();
	fdi.close();
}
