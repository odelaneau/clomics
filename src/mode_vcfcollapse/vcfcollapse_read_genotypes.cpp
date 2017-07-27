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

#include "vcfcollapse_data.h"

void vcfcollapse_data::readGenotypes(string fvcf) {
	int n_includedG = 0;
	int n_excludedG_mult = 0;
	int n_excludedG_void = 0;
	int n_excludedG_user = 0;
	int n_excludedG_chr = 0;
	int n_excludedG_ann = 0;
	int n_includedS = 0;

	vector < int > mappingS;

	vrb.title("Reading genotype in [" + fvcf + "]");

	//Opening files
	bcf_srs_t * sr =  bcf_sr_init();
	if(!(bcf_sr_add_reader (sr, fvcf.c_str()))) {
		switch (sr->errnum) {
		case not_bgzf: vrb.error("File not compressed with bgzip!"); break;
		case idx_load_failed: vrb.error("Impossible to load index file!"); break;
		case file_type_error: vrb.error("File format not detected by htslib!"); break;
		default : vrb.error("Unknown error!");
		}
	}

	//Sample processing
	int n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i0 = 0 ; i0 < n_samples ; i0 ++) {
		mappingS.push_back(findSample(string(sr->readers[0].header->samples[i0])));
		if (mappingS.back() >= 0) n_includedS++;
	}

	//Read genotype data
	unsigned int linecount=0;
	int ngt, ngt_arr = 0, nds, nds_arr = 0, * gt_arr = NULL, nsl, nsl_arr = 0, * sl_arr = NULL;
	float * ds_arr = NULL;
	bcf1_t * line;
	while(bcf_sr_next_line (sr)) {
        linecount ++;
        if (linecount % 100000 == 0) vrb.bullet("Read " + stb.str(linecount) + " lines");
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele == 2) {
			bcf_unpack(line, BCF_UN_STR);
			string sid = string(line->d.id);
			string chr = string(bcf_hdr_id2name(sr->readers[0].header, line->rid));
			int pos = line->pos + 1;

			map < string, int >::iterator itC = ann_chr_map.find(chr);
			if (itC != ann_chr_map.end()) {

				vector < Interval < int > > ann_in_cis;
				ann_int[itC->second].findOverlapping(pos, ann_in_cis);

				if (ann_in_cis.size() > 0) {
					ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);
					nds = bcf_get_format_float(sr->readers[0].header, line,"DS", &ds_arr, &nds_arr);
					if (nds == n_samples || ngt == 2*n_samples) {
						if (filter_genotype.check(sid)) {
							genotype_id.push_back(sid);
							genotype_chr.push_back(chr);
							string genotype_ref = string(line->d.allele[0]);
							genotype_start.push_back(pos);
							nsl = bcf_get_info_int32(sr->readers[0].header, line, "END", &sl_arr, &nsl_arr);
							if (nsl >= 0 && nsl_arr == 1) genotype_end.push_back(sl_arr[0]);
							else genotype_end.push_back(genotype_start.back() + genotype_ref.size() - 1);
							genotype_bed.push_back(ann_id[ann_in_cis[0].value]);
							genotype_val.push_back(vector < float > (sample_count, 0.0));

							for(int i = 0 ; i < n_samples ; i ++) {
								if (mappingS[i] >= 0) {
									if (nds > 0) genotype_val.back()[mappingS[i]] = ds_arr[i];
									else {
										if (gt_arr[2*i+0] == bcf_gt_missing || gt_arr[2*i+1] == bcf_gt_missing) genotype_val.back()[mappingS[i]] = bcf_float_missing;
										else genotype_val.back()[mappingS[i]] = bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]);
									}
								}
							}
							n_includedG++;
						} else n_excludedG_user ++;
					} else n_excludedG_void ++;
				} else n_excludedG_ann ++;
			} else n_excludedG_chr ++;
		} else n_excludedG_mult ++;
	}
	//Finalize
	free(gt_arr);
	free(ds_arr);
	bcf_sr_destroy(sr);
	genotype_count = n_includedG;
	vrb.bullet(stb.str(n_includedG) + " variants included");
	if (n_excludedG_user > 0) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
	if (n_excludedG_mult > 0) vrb.bullet(stb.str(n_excludedG_mult) + " multi-allelic variants excluded");
	if (n_excludedG_void > 0) vrb.bullet(stb.str(n_excludedG_void) + " uninformative variants excluded [no GT/DS]");
	if (n_excludedG_chr > 0) vrb.bullet(stb.str(n_excludedG_chr) + " variants excluded because chr is unfound");
	if (n_excludedG_ann > 0) vrb.bullet(stb.str(n_excludedG_ann) + " non-annotated variants excluded");
	if (genotype_count == 0) vrb.leave("Cannot find genotypes in target region!");

}
