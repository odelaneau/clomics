#include "motif_cis_data.h"

void motif_cis_data::readMotives(string fbed, bool perm) {
	vector < string > tokens;

	//Open BED file
	vrb.title("Reading TF binding motives in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};

	//Read phenotypes
	hts_itr_t *itr = tbx_itr_querys(tbx, regionPhenotype.get().c_str());
	vrb.bullet("target region [" + regionPhenotype.get() + "]");
	if (!itr) vrb.error("Cannot jump to region!");
	//Read data
	while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
		stb.split(string(str.s), tokens);
		if (tokens.size() < 4) vrb.error("Incorrect number of columns!");
		motif_cis_id.push_back(tokens[3]);
		motif_cis_start.push_back(atoi(tokens[1].c_str()));
		motif_cis_stop.push_back(atoi(tokens[2].c_str()));
	}
	tbx_itr_destroy(itr);

	//
	if (perm) shuffle(motif_cis_id.begin(), motif_cis_id.end(), rng.getEngine());

	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(motif_cis_id.size()) + " motives read");
}

void motif_cis_data::initializeMotives() {
	vrb.title("Initialize required structures for TF binding data");
	//
	map < string , int >::iterator motif_cis_it;
	for (int t = 0 ; t < motif_cis_id.size() ; t ++) {
		motif_cis_it = motif_cis_map.find(motif_cis_id[t]);
		if (motif_cis_it == motif_cis_map.end()) motif_cis_map.insert(pair<string, int> (motif_cis_id[t], motif_cis_map.size()));
	}
	vrb.bullet(stb.str(motif_cis_map.size()) + " distinct motives found");

	//
	motif_cis_idx = vector < int >(motif_cis_id.size(), -1);
	for (int t = 0 ; t < motif_cis_id.size() ; t ++) {
		motif_cis_it = motif_cis_map.find(motif_cis_id[t]);
		motif_cis_idx[t] = motif_cis_it->second;
	}
	//
	motif_cis_unique_id = vector < string >(motif_cis_map.size());
	for (map < string, int > :: iterator it = motif_cis_map.begin(); it != motif_cis_map.end() ; ++ it) motif_cis_unique_id[it->second] = it->first;
	vrb.bullet("indexing motives is done");
	//
	interaction_cnt = vector < vector < int > > (motif_cis_map.size(), vector < int > (motif_cis_map.size(), 0));
	interaction_cor = vector < vector < double > > (motif_cis_map.size(), vector < double > (motif_cis_map.size(), 0));
}

void motif_cis_data::mappingMotives() {
	vrb.title("Mapping TF motives to phenotypes");
	//
	phenotype_motives = vector < vector < int > > (phenotype_count);
	//
	int n_mapped = 0, n_unmapped = 0;
	for (int t = 0 ; t < motif_cis_id.size() ; t ++) {
		vector < Interval < int > > phenotype_overlap;
		phenotype_tree.findOverlapping(motif_cis_start[t], motif_cis_stop[t], phenotype_overlap);
		for (int o = 0 ; o < phenotype_overlap.size() ; o ++) {
			phenotype_motives[phenotype_overlap[o].value].push_back(motif_cis_idx[t]);
		}
		if (phenotype_overlap.size() > 0) n_mapped ++;
		else n_unmapped++;
	}
	vrb.bullet(stb.str(n_mapped) + " mapped / " + stb.str(n_unmapped) + " unmapped");
}
