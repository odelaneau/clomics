/*Copyright (C) 2015 Olivier Delaneau, Halit Ongen, Emmanouil T. Dermitzakis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include "motif_trans_data.h"

void motif_trans_data::scanPhenotypes1(string fbed) {
	int n_includedP = 0;
	int n_excludedP = 0;
	vector < string > tokens;

	//Open BED file
	vrb.title("Scanning phenotype data 1 in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line!");

	//Read phenotypes
	if (regionPhenotype1.chr != "NA"){
		hts_itr_t *itr = tbx_itr_querys(tbx, regionPhenotype1.get().c_str());
		vrb.bullet("target region [" + regionPhenotype1.get() + "]");
		if (!itr) vrb.error("Cannot jump to region!");
		//Read data
		while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
			stb.split(string(str.s), tokens);
			if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
			if (filter_phenotype.check(tokens[3])) {
				phenotype_vec1.push_back(motif_trans_phenotype_data(tokens[3], tokens[0], atoi(tokens[1].c_str()) + 1, atoi(tokens[2].c_str())));
				n_includedP++;
			} else n_excludedP ++;
		}
		tbx_itr_destroy(itr);
	} else {
		while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
			stb.split(string(str.s), tokens);
			if (str.l && str.s[0] != tbx->conf.meta_char) {
				if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
				if (filter_phenotype.check(tokens[3])) {
					phenotype_vec1.push_back(motif_trans_phenotype_data(tokens[3], tokens[0], atoi(tokens[1].c_str()) + 1, atoi(tokens[2].c_str())));
					n_includedP++;
				} else n_excludedP ++;
			}
		}
	}

	//Finalize & verbose
	tbx_destroy(tbx);
	phenotype_count1 = n_includedP;
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded by user");
}

void motif_trans_data::scanPhenotypes2(string fbed) {
	int n_includedP = 0;
	int n_excludedP = 0;
	vector < string > tokens;

	//Open BED file
	vrb.title("Scanning phenotype data 2 in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line!");

	//Read phenotypes
	if (regionPhenotype2.chr != "NA"){
		hts_itr_t *itr = tbx_itr_querys(tbx, regionPhenotype2.get().c_str());
		vrb.bullet("target region [" + regionPhenotype2.get() + "]");
		if (!itr) vrb.error("Cannot jump to region!");
		//Read data
		while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
			stb.split(string(str.s), tokens);
			if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
			if (filter_phenotype.check(tokens[3])) {
				phenotype_vec2.push_back(motif_trans_phenotype_data(tokens[3], tokens[0], atoi(tokens[1].c_str()) + 1, atoi(tokens[2].c_str())));
				n_includedP++;
			} else n_excludedP ++;
		}
		tbx_itr_destroy(itr);
	} else {
		while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
			stb.split(string(str.s), tokens);
			if (str.l && str.s[0] != tbx->conf.meta_char) {
				if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
				if (filter_phenotype.check(tokens[3])) {
					phenotype_vec2.push_back(motif_trans_phenotype_data(tokens[3], tokens[0], atoi(tokens[1].c_str()) + 1, atoi(tokens[2].c_str())));
					n_includedP++;
				} else n_excludedP ++;
			}
		}
	}

	//Finalize & verbose
	tbx_destroy(tbx);
	phenotype_count2 = n_includedP;
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded by user");
}

void motif_trans_data::sortPhenotypes1() {
	vrb.title("Sorting phenotype data 1 and allocating required memory");
	sort(phenotype_vec1.begin(), phenotype_vec1.end());
	for (int p = 0 ; p < phenotype_vec1.size() ; p ++) phenotype_map1.insert(pair < string , unsigned int > (phenotype_vec1[p].id, p));
	phenotype_count1 = phenotype_vec1.size();
	phenotype_val1 = vector < vector < float > > (phenotype_count1, vector < float > (sample_count, 0.0));
	phenotype_id1 = vector < string > (phenotype_count1);
	phenotype_chr1 = vector < string > (phenotype_count1);
	phenotype_start1 = vector < int > (phenotype_count1, -1);
	phenotype_end1 = vector < int > (phenotype_count1, -1);
	vrb.bullet("#samples = " + stb.str(sample_count) + " #phenotypes = " + stb.str(phenotype_count1));
}

void motif_trans_data::sortPhenotypes2() {
	vrb.title("Sorting phenotype data 2 and allocating required memory");
	sort(phenotype_vec2.begin(), phenotype_vec2.end());
	for (int p = 0 ; p < phenotype_vec2.size() ; p ++) phenotype_map2.insert(pair < string , unsigned int > (phenotype_vec2[p].id, p));
	phenotype_count2 = phenotype_vec2.size();
	phenotype_val2 = vector < vector < float > > (phenotype_count2, vector < float > (sample_count, 0.0));
	phenotype_id2 = vector < string > (phenotype_count2);
	phenotype_chr2 = vector < string > (phenotype_count2);
	phenotype_start2 = vector < int > (phenotype_count2, -1);
	phenotype_end2 = vector < int > (phenotype_count2, -1);
	vrb.bullet("#samples = " + stb.str(sample_count) + " #phenotypes = " + stb.str(phenotype_count2));
}

void motif_trans_data::readPhenotypes1(string fbed) {
	int n_includedS = 0;
	int n_includedP = 0;
	int n_excludedP = 0;
	int n_negativeStrd = 0;
	vector < int > mappingS;
	unordered_map < string, unsigned int >::iterator it_map;

	//Open BED file
	vrb.title("Reading phenotype data 1 in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line!");

	//Process sample names
	vector < string > tokens;
	stb.split(string(str.s), tokens);
	if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
	for (int t = 6 ; t < tokens.size() ; t ++) {
		mappingS.push_back(findSample(tokens[t]));
		if (mappingS.back() >= 0) n_includedS++;
	}

    //Read phenotypes
    if (regionPhenotype1.chr != "NA"){
        hts_itr_t *itr = tbx_itr_querys(tbx, regionPhenotype1.get().c_str());
        vrb.bullet("target region [" + regionPhenotype1.get() + "]");
        if (!itr) vrb.error("Cannot jump to region!");
        //Read data
        while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
            stb.split(string(str.s), tokens);
            if (tokens.size() < 5) vrb.error("Incorrect number of columns!");
            it_map = phenotype_map1.find(tokens[3]);
            if (it_map != phenotype_map1.end()) {
            	phenotype_id1[it_map->second] = tokens[3];
                phenotype_chr1[it_map->second] = tokens[0];
                phenotype_start1[it_map->second] = atoi(tokens[1].c_str()) + 1;
                phenotype_end1[it_map->second] = atoi(tokens[2].c_str());
                phenotype_neg1.push_back(tokens[5] == "-");
                if (phenotype_neg1.back()) n_negativeStrd ++;
                for (int t = 6 ; t < tokens.size() ; t ++) {
                    if (mappingS[t-6] >= 0) {
                        if (tokens[t] == "NA") phenotype_val1[it_map->second][mappingS[t-6]] = bcf_float_missing;
                        else phenotype_val1[it_map->second][mappingS[t-6]] = stof(tokens[t]);
                    }
                }
                n_includedP++;
            } else n_excludedP ++;
        }
        tbx_itr_destroy(itr);
    }else{
        while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
            stb.split(string(str.s), tokens);
            if (str.l && str.s[0] != tbx->conf.meta_char) {
                if (tokens.size() < 5) vrb.error("Incorrect number of columns!");
                it_map = phenotype_map1.find(tokens[3]);
                if (it_map != phenotype_map1.end()) {
                	phenotype_id1[it_map->second] = tokens[3];
                    phenotype_chr1[it_map->second] = tokens[0];
                    phenotype_start1[it_map->second] = atoi(tokens[1].c_str()) + 1;
                    phenotype_end1[it_map->second] = atoi(tokens[2].c_str());
                    phenotype_neg1.push_back(tokens[5] == "-");
                    if (phenotype_neg1.back()) n_negativeStrd ++;
                    for (int t = 6 ; t < tokens.size() ; t ++) {
                        if (mappingS[t-6] >= 0) {
                            if (tokens[t] == "NA") phenotype_val1[it_map->second][mappingS[t-6]] = bcf_float_missing;
                            else phenotype_val1[it_map->second][mappingS[t-6]] = stof(tokens[t]);
                        }
                    }
                    n_includedP++;
                } else n_excludedP ++;
            }
        }
    }

	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded");
	if (n_negativeStrd > 0) vrb.bullet(stb.str(n_negativeStrd) + " phenotypes on the negative strand");
    if (phenotype_count1 == 0) vrb.leave("Cannot find phenotypes in target region!");
}

void motif_trans_data::readPhenotypes2(string fbed) {
	int n_includedS = 0;
	int n_includedP = 0;
	int n_excludedP = 0;
	int n_negativeStrd = 0;
	vector < int > mappingS;
	unordered_map < string, unsigned int >::iterator it_map;

	//Open BED file
	vrb.title("Reading phenotype data 2 in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line!");

	//Process sample names
	vector < string > tokens;
	stb.split(string(str.s), tokens);
	if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
	for (int t = 6 ; t < tokens.size() ; t ++) {
		mappingS.push_back(findSample(tokens[t]));
		if (mappingS.back() >= 0) n_includedS++;
	}

    //Read phenotypes
    if (regionPhenotype2.chr != "NA"){
        hts_itr_t *itr = tbx_itr_querys(tbx, regionPhenotype2.get().c_str());
        vrb.bullet("target region [" + regionPhenotype2.get() + "]");
        if (!itr) vrb.error("Cannot jump to region!");
        //Read data
        while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
            stb.split(string(str.s), tokens);
            if (tokens.size() < 5) vrb.error("Incorrect number of columns!");
            it_map = phenotype_map2.find(tokens[3]);
            if (it_map != phenotype_map2.end()) {
            	phenotype_id2[it_map->second] = tokens[3];
                phenotype_chr2[it_map->second] = tokens[0];
                phenotype_start2[it_map->second] = atoi(tokens[1].c_str()) + 1;
                phenotype_end2[it_map->second] = atoi(tokens[2].c_str());
                phenotype_neg2.push_back(tokens[5] == "-");
                if (phenotype_neg2.back()) n_negativeStrd ++;
                for (int t = 6 ; t < tokens.size() ; t ++) {
                    if (mappingS[t-6] >= 0) {
                        if (tokens[t] == "NA") phenotype_val2[it_map->second][mappingS[t-6]] = bcf_float_missing;
                        else phenotype_val2[it_map->second][mappingS[t-6]] = stof(tokens[t]);
                    }
                }
                n_includedP++;
            } else n_excludedP ++;
        }
        tbx_itr_destroy(itr);
    }else{
        while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
            stb.split(string(str.s), tokens);
            if (str.l && str.s[0] != tbx->conf.meta_char) {
                if (tokens.size() < 5) vrb.error("Incorrect number of columns!");
                it_map = phenotype_map2.find(tokens[3]);
                if (it_map != phenotype_map2.end()) {
                	phenotype_id2[it_map->second] = tokens[3];
                    phenotype_chr2[it_map->second] = tokens[0];
                    phenotype_start2[it_map->second] = atoi(tokens[1].c_str()) + 1;
                    phenotype_end2[it_map->second] = atoi(tokens[2].c_str());
                    phenotype_neg2.push_back(tokens[5] == "-");
                    if (phenotype_neg2.back()) n_negativeStrd ++;
                    for (int t = 6 ; t < tokens.size() ; t ++) {
                        if (mappingS[t-6] >= 0) {
                            if (tokens[t] == "NA") phenotype_val2[it_map->second][mappingS[t-6]] = bcf_float_missing;
                            else phenotype_val2[it_map->second][mappingS[t-6]] = stof(tokens[t]);
                        }
                    }
                    n_includedP++;
                } else n_excludedP ++;
            }
        }
    }

	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded");
	if (n_negativeStrd > 0) vrb.bullet(stb.str(n_negativeStrd) + " phenotypes on the negative strand");
    if (phenotype_count2 == 0) vrb.leave("Cannot find phenotypes in target region!");
}


void motif_trans_data::buildIntervalTree1() {
	vrb.title("Building interval tree for phenotype data 1");
	vector < Interval < int > > interval_vec;
	for (int p = 0 ; p < phenotype_start1.size() ; p ++) interval_vec.push_back(Interval < int > (phenotype_start1[p], phenotype_end1[p], p));
	phenotype_tree1 = IntervalTree < int > (interval_vec);
	vrb.bullet("done!");
}

void motif_trans_data::buildIntervalTree2() {
	vrb.title("Building interval tree for phenotype data 2");
	vector < Interval < int > > interval_vec;
	for (int p = 0 ; p < phenotype_start2.size() ; p ++) interval_vec.push_back(Interval < int > (phenotype_start2[p], phenotype_end2[p], p));
	phenotype_tree2 = IntervalTree < int > (interval_vec);
	vrb.bullet("done!");
}
