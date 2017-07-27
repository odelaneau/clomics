#include "span_data.h"

void span_data::readAnnotations(string fbed) {
	vector < string > tokens;

	//Open BED file
	vrb.title("Reading annotations in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};

	//Read phenotypes
	hts_itr_t *itr = tbx_itr_querys(tbx, regionPhenotype.get().c_str());
	//vrb.bullet("target region [" + regionPhenotype.get() + "]");
	if (!itr) vrb.error("Cannot jump to region!");

	//Read data
	map < string, int >::iterator it;
	vector < Interval < int > > interval_vec;
	while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
		stb.split(string(str.s), tokens);
		if (tokens.size() < 4) vrb.error("Incorrect number of columns!");
		it = annotation_id_map.find(tokens[3]);
		if (it == annotation_id_map.end()) {
			interval_vec.push_back(Interval < int > (atoi(tokens[1].c_str()), atoi(tokens[2].c_str()), annotation_id_vec.size()));
			annotation_id_map.insert(pair < string, int > (tokens[3], annotation_id_vec.size()));
			annotation_id_vec.push_back(tokens[3]);
		} else interval_vec.push_back(Interval < int > (atoi(tokens[1].c_str()), atoi(tokens[2].c_str()), it->second));
	}
	annotation_tree = IntervalTree < int > (interval_vec);
	tbx_itr_destroy(itr);

	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(interval_vec.size()) + " annotations read");
}
