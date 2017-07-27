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

void vcfcollapse_data::writeGenotypes(string fvcf) {
	vrb.title("Writing VCF file [" + fvcf + "]");

    //VCF OPEN
	htsFile * bcf_fd = NULL;
	bool compressed_vcf = fvcf.substr(fvcf.find_last_of(".") + 1) == "gz";
	if (compressed_vcf) bcf_fd = bcf_open(fvcf.c_str(), "wz");
	else bcf_fd = bcf_open(fvcf.c_str(), "wu");
	if (bcf_fd == NULL) vrb.error("Impossible to create VCF file");
	else if (compressed_vcf) vrb.bullet("BGZIP for VCF compression is ON");
	else vrb.bullet("BGZIP compression for VCF is OFF (add .gz to filename to activate it)");
	bcf_hdr_t * bcf_hdr = bcf_hdr_init("w");
	kstring_t str = {0,0,0};
	ksprintf(&str, "##QTLtools vcfcollapse Version=1\n");
	bcf_hdr_append(bcf_hdr, str.s);
    free(str.s);

	//VCF INFO
	vrb.bullet("Writing VCF header [INFO, CONTIG, FORMAT, SAMPLES]");
	bcf_hdr_append(bcf_hdr,"##INFO=<ID=END,Number=1,Type=Integer,Description=\"END POSITION\">");

	//VCF CONTIG
	vrb.bullet("Writing " + stb.str(ann_chr_map.size()) + " CONTIG fields");
	for (map < string , int > :: iterator it_c = ann_chr_map.begin() ; it_c != ann_chr_map.end() ; ++ it_c) {
		string tmp_str = "##contig=<ID=" + it_c->first + ">";
		bcf_hdr_append(bcf_hdr, tmp_str.c_str());
	}

	//VCF FORMAT
    bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"RARE VARIANT DOSAGE\">");

    //VCF SAMPLES
    for (int i = 0 ; i < sample_count ; i ++) bcf_hdr_add_sample(bcf_hdr, sample_id[i].c_str());
    bcf_hdr_add_sample(bcf_hdr, NULL);

    //VCF HEADER
    bcf_hdr_write(bcf_fd, bcf_hdr);

	//Loop across BURDEN
    int n_rm_maf = 0;
    bcf1_t * bcf_rec = bcf_init1();
    float * DS = (float*) malloc(sample_count* 1 * sizeof(float));
    for (int b = 0 ; b < burden_count ; b ++) {

    	double maf = 0.0;
    	for (int i = 0 ; i < sample_count ; i ++) if (burden_val[b][i] > 0) maf += 1.0;
    	maf /= sample_count;

    	if (maf >= param_min_maf) {
			bcf_clear1(bcf_rec);
			bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, burden_chr[b].c_str());
			bcf_rec->pos = burden_start[b];
			bcf_update_id(bcf_hdr, bcf_rec, burden_id[b].c_str());
			string str_alleles = "R,A";
			bcf_update_alleles_str(bcf_hdr, bcf_rec, str_alleles.c_str());
			bcf_rec->qual = 100;
			int tmpi = bcf_hdr_id2int(bcf_hdr, BCF_DT_ID, "PASS");
			bcf_update_filter(bcf_hdr, bcf_rec, &tmpi, 1);
			tmpi = burden_end[b];
			bcf_update_info_int32(bcf_hdr, bcf_rec, "END", &tmpi, 1);
			for (int i = 0; i < sample_count ; i ++) DS[i] = burden_val[b][i];
			bcf_update_format_float(bcf_hdr, bcf_rec, "DS", DS, sample_count);
			bcf_write1(bcf_fd, bcf_hdr, bcf_rec);
    	} else n_rm_maf++;
    }
    vrb.bullet("#remove because low maf = " + stb.str(n_rm_maf));
    free(DS);
    hts_close(bcf_fd);
    bcf_destroy1(bcf_rec);
    bcf_hdr_destroy(bcf_hdr);
}
