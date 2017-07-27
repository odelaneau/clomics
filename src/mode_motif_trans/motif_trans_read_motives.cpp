#include "motif_trans_data.h"

void motif_trans_data::readMotives1(string fbed, bool perm) {
	vector < string > tokens;

	//Open BED file
	vrb.title("Reading TF binding motives 1 in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};

	//Read phenotypes
	hts_itr_t *itr = tbx_itr_querys(tbx, regionPhenotype1.get().c_str());
	vrb.bullet("target region [" + regionPhenotype1.get() + "]");
	if (!itr) vrb.error("Cannot jump to region!");
	//Read data
	while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
		stb.split(string(str.s), tokens);
		if (tokens.size() < 4) vrb.error("Incorrect number of columns!");
		motif_trans_id1.push_back(tokens[3]);
		motif_trans_start1.push_back(atoi(tokens[1].c_str()));
		motif_trans_stop1.push_back(atoi(tokens[2].c_str()));
	}
	tbx_itr_destroy(itr);

	//
	if (perm) shuffle(motif_trans_id1.begin(), motif_trans_id1.end(), rng.getEngine());

	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(motif_trans_id1.size()) + " motives read");
}

void motif_trans_data::readMotives2(string fbed, bool perm) {
	vector < string > tokens;

	//Open BED file
	vrb.title("Reading TF binding motives 2 in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};

	//Read phenotypes
	hts_itr_t *itr = tbx_itr_querys(tbx, regionPhenotype2.get().c_str());
	vrb.bullet("target region [" + regionPhenotype2.get() + "]");
	if (!itr) vrb.error("Cannot jump to region!");
	//Read data
	while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
		stb.split(string(str.s), tokens);
		if (tokens.size() < 4) vrb.error("Incorrect number of columns!");
		motif_trans_id2.push_back(tokens[3]);
		motif_trans_start2.push_back(atoi(tokens[1].c_str()));
		motif_trans_stop2.push_back(atoi(tokens[2].c_str()));
	}
	tbx_itr_destroy(itr);

	//
	if (perm) shuffle(motif_trans_id2.begin(), motif_trans_id2.end(), rng.getEngine());

	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(motif_trans_id2.size()) + " motives read");
}

void motif_trans_data::initializeMotives() {
	vrb.title("Initialize required structures for TF binding data");
	//
	map < string , int >::iterator motif_trans_it;
	for (int t = 0 ; t < motif_trans_id1.size() ; t ++) {
		motif_trans_it = motif_trans_map.find(motif_trans_id1[t]);
		if (motif_trans_it == motif_trans_map.end()) motif_trans_map.insert(pair<string, int> (motif_trans_id1[t], motif_trans_map.size()));
	}
	for (int t = 0 ; t < motif_trans_id2.size() ; t ++) {
		motif_trans_it = motif_trans_map.find(motif_trans_id2[t]);
		if (motif_trans_it == motif_trans_map.end()) motif_trans_map.insert(pair<string, int> (motif_trans_id2[t], motif_trans_map.size()));
	}
	vrb.bullet(stb.str(motif_trans_map.size()) + " distinct motives found");

	//
	motif_trans_idx1 = vector < int >(motif_trans_id1.size(), -1);
	motif_trans_idx2 = vector < int >(motif_trans_id2.size(), -1);
	for (int t = 0 ; t < motif_trans_id1.size() ; t ++) {
		motif_trans_it = motif_trans_map.find(motif_trans_id1[t]);
		motif_trans_idx1[t] = motif_trans_it->second;
	}
	for (int t = 0 ; t < motif_trans_id2.size() ; t ++) {
		motif_trans_it = motif_trans_map.find(motif_trans_id2[t]);
		motif_trans_idx2[t] = motif_trans_it->second;
	}
	//
	motif_trans_unique_id = vector < string >(motif_trans_map.size());
	for (map < string, int > :: iterator it = motif_trans_map.begin(); it != motif_trans_map.end() ; ++ it) motif_trans_unique_id[it->second] = it->first;
	vrb.bullet("indexing motives is done");
	//
	interaction_cnt = vector < vector < int > > (motif_trans_map.size(), vector < int > (motif_trans_map.size(), 0));
	interaction_cor = vector < vector < double > > (motif_trans_map.size(), vector < double > (motif_trans_map.size(), 0));
}

void motif_trans_data::mappingMotives1() {
	vrb.title("Mapping TF motives to phenotypes 1");
	//
	phenotype_motives1 = vector < vector < int > > (phenotype_count1);
	//
	int n_mapped = 0, n_unmapped = 0;
	for (int t = 0 ; t < motif_trans_id1.size() ; t ++) {
		vector < Interval < int > > phenotype_overlap;
		phenotype_tree1.findOverlapping(motif_trans_start1[t], motif_trans_stop1[t], phenotype_overlap);
		for (int o = 0 ; o < phenotype_overlap.size() ; o ++) {
			phenotype_motives1[phenotype_overlap[o].value].push_back(motif_trans_idx1[t]);
		}
		if (phenotype_overlap.size() > 0) n_mapped ++;
		else n_unmapped++;
	}
	vrb.bullet(stb.str(n_mapped) + " mapped / " + stb.str(n_unmapped) + " unmapped");
}

void motif_trans_data::mappingMotives2() {
	vrb.title("Mapping TF motives to phenotypes 2");
	//
	phenotype_motives2 = vector < vector < int > > (phenotype_count2);
	//
	int n_mapped = 0, n_unmapped = 0;
	for (int t = 0 ; t < motif_trans_id2.size() ; t ++) {
		vector < Interval < int > > phenotype_overlap;
		phenotype_tree2.findOverlapping(motif_trans_start2[t], motif_trans_stop2[t], phenotype_overlap);
		for (int o = 0 ; o < phenotype_overlap.size() ; o ++) {
			phenotype_motives2[phenotype_overlap[o].value].push_back(motif_trans_idx2[t]);
		}
		if (phenotype_overlap.size() > 0) n_mapped ++;
		else n_unmapped++;
	}
	vrb.bullet(stb.str(n_mapped) + " mapped / " + stb.str(n_unmapped) + " unmapped");
}
