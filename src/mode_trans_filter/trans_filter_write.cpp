#include "trans_filter_data.h"

void trans_filter_data::writeInteractions(string file) {
	vrb.title("Writing merged interactions [" + file + "]");
	output_file fd (file.c_str());

	for (int e = 0 ; e < inter_chr.size() ; e ++)
		fd << inter_chr[e].first << " " << (int)((inter_start[e].first + inter_end[e].first) / 2) << " " << inter_chr[e].second << " " << (int)((inter_start[e].second + inter_end[e].second) / 2) << " " << inter_end[e].first - inter_start[e].first << " " << inter_end[e].second - inter_start[e].second << endl;
	fd.close();
}
