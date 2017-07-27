#include "trans_filter_data.h"

void trans_filter_data::readInteractions(string file) {
	string buffer;
	long n_total = 0, n_inter = 0, n_merged = 0;
	vector < string > tokens;
	map < string, int >::iterator it1, it2;

	vrb.title("Reading interactions [" + file + "]");
	input_file fd (file.c_str());
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		assert(tokens.size() == 9);
		double chic = atof(tokens[6].c_str());
		double cpva = atof(tokens[8].c_str());
		if (chic >= hic && cpva <= pva) {
			//cout << buffer << endl;
			string chr1 = tokens[1];
			string chr2 = tokens[4];
			int start1 = atoi(tokens[2].c_str()) - bin;
			if (start1 < 0) start1 = 0;
			int end1 = atoi(tokens[2].c_str()) + bin;
			int start2 = atoi(tokens[5].c_str()) - bin;
			if (start2 < 0) start2 = 0;
			int end2 = atoi(tokens[5].c_str()) + bin;

			int idx = -1;
			for (int e = 0; e < inter_start.size() && idx < 0; e ++) {
				if (inter_chr[e].first == chr1 && inter_chr[e].second == chr2 && inter_start[e].first < end1 && start1 < inter_end[e].first && inter_start[e].second < end2 && start2 < inter_end[e].second) idx = e;
			}

			if (idx >= 0) {
				if (start1 < inter_start[idx].first) inter_start[idx].first = start1;
				if (end1 > inter_end[idx].first) inter_end[idx].first = end1;
				if (start2 < inter_start[idx].second) inter_start[idx].second = start2;
				if (end2 > inter_end[idx].second) inter_end[idx].second = end2;
			} else {
				inter_chr.push_back(pair < string, string > (chr1, chr2));
				inter_start.push_back(pair < int, int > (start1, start2));
				inter_end.push_back(pair < int, int > (end1, end2));
				n_merged ++;
			}
			n_inter ++;
		}
		n_total ++;
	}
	fd.close();
	vrb.bullet("Number of interactions read = " + stb.str(n_total));
	vrb.bullet("Number of interactions filtered = " + stb.str(n_inter));
	vrb.bullet("Number of interactions merged = " + stb.str(n_merged));
}
