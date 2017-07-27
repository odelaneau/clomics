#include "span_data.h"

void span_data::scanPhenotypes(string fbed) {
	int n_includedP = 0;
	int n_excludedP = 0;
	vector < string > tokens;

	//Open BED file
	vrb.title("Scanning phenotype data in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line!");

	//Read phenotypes
	if (regionPhenotype.chr != "NA"){
		hts_itr_t *itr = tbx_itr_querys(tbx, regionPhenotype.get().c_str());
		vrb.bullet("target region [" + regionPhenotype.get() + "]");
		if (!itr) vrb.error("Cannot jump to region!");
		//Read data
		while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
			stb.split(string(str.s), tokens);
			if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
			if (filter_phenotype.check(tokens[3])) {
				phenotype_vec.push_back(span_phenotype_data(tokens[3], tokens[0], atoi(tokens[1].c_str()) + 1, atoi(tokens[2].c_str())));
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
					phenotype_vec.push_back(span_phenotype_data(tokens[3], tokens[0], atoi(tokens[1].c_str()) + 1, atoi(tokens[2].c_str())));
					n_includedP++;
				} else n_excludedP ++;
			}
		}
	}

	//Finalize & verbose
	tbx_destroy(tbx);
	phenotype_count = n_includedP;
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded by user");
}

void span_data::sortPhenotypes() {
	vrb.title("Sorting phenotype data and allocating required memory");
	sort(phenotype_vec.begin(), phenotype_vec.end());
	for (int p = 0 ; p < phenotype_vec.size() ; p ++) phenotype_map.insert(pair < string , unsigned int > (phenotype_vec[p].id, p));
	phenotype_count = phenotype_vec.size();
	phenotype_val = vector < vector < float > > (phenotype_count, vector < float > (sample_count, 0.0));
	phenotype_id = vector < string > (phenotype_count);
	phenotype_chr = vector < string > (phenotype_count);
	for (int p = 0 ; p < phenotype_vec.size() ; p ++) {
		phenotype_start.push_back(phenotype_vec[p].start);
		phenotype_end.push_back(phenotype_vec[p].end);
	}
	vrb.bullet("#samples = " + stb.str(sample_count) + " #phenotypes = " + stb.str(phenotype_count));
}


void span_data::readPhenotypes(string fbed) {
	int n_includedS = 0;
	int n_includedP = 0;
	int n_excludedP = 0;
	int n_negativeStrd = 0;
	vector < int > mappingS;
	unordered_map < string, unsigned int >::iterator it_map;

	//Open BED file
	vrb.title("Reading phenotype data in [" + fbed + "]");
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
    if (regionPhenotype.chr != "NA"){
        hts_itr_t *itr = tbx_itr_querys(tbx, regionPhenotype.get().c_str());
        vrb.bullet("target region [" + regionPhenotype.get() + "]");
        if (!itr) vrb.error("Cannot jump to region!");
        //Read data
        while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
            stb.split(string(str.s), tokens);
            if (tokens.size() < 5) vrb.error("Incorrect number of columns!");
            it_map = phenotype_map.find(tokens[3]);
            if (it_map != phenotype_map.end()) {
            	phenotype_id[it_map->second] = tokens[3];
                phenotype_chr[it_map->second] = tokens[0];
                phenotype_start[it_map->second] = atoi(tokens[1].c_str()) + 1;
                phenotype_end[it_map->second] = atoi(tokens[2].c_str());
                for (int t = 6 ; t < tokens.size() ; t ++) {
                    if (mappingS[t-6] >= 0) {
                        if (tokens[t] == "NA") phenotype_val[it_map->second][mappingS[t-6]] = bcf_float_missing;
                        else phenotype_val[it_map->second][mappingS[t-6]] = stof(tokens[t]);
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
                it_map = phenotype_map.find(tokens[3]);
                if (it_map != phenotype_map.end()) {
                	phenotype_id[it_map->second] = tokens[3];
                    phenotype_chr[it_map->second] = tokens[0];
                    phenotype_start[it_map->second] = atoi(tokens[1].c_str()) + 1;
                    phenotype_end[it_map->second] = atoi(tokens[2].c_str());
                    for (int t = 6 ; t < tokens.size() ; t ++) {
                        if (mappingS[t-6] >= 0) {
                            if (tokens[t] == "NA") phenotype_val[it_map->second][mappingS[t-6]] = bcf_float_missing;
                            else phenotype_val[it_map->second][mappingS[t-6]] = stof(tokens[t]);
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
    if (phenotype_count == 0) vrb.leave("Cannot find phenotypes in target region!");
}

