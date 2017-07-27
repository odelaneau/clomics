#include "span_data.h"

void span_data::computeInteractionsSingle(string fout, int max_dist, bool mode_peak) {
	vrb.title("Report 1 all possible interactions in [" + fout + "]");
	output_file fd(fout.c_str());

	//
	vector < pair < string, int > > mapping;
	for ( int a = 0 ; a < annotation_id_vec.size() ; a ++) mapping.push_back(pair < string, int > (annotation_id_vec[a], a));
	sort(mapping.begin(), mapping.end());

	//
	fd << "IDX1 CHR1 POS1 ID1 IDX2 CHR2 POS2 ID2 CORR PVAL";
	for ( int a = 0 ; a < annotation_id_vec.size() ; a ++) fd << " " << mapping[a].first;
	fd << endl;

	//
	for (int i = 0 ; i < phenotype_count - 1 ; i ++) {
		for (int j = i + 1 ; j < phenotype_count ; j ++) {
			int posi =  (phenotype_start[i]+phenotype_end[i])/2;
			int posj =  (phenotype_start[j]+phenotype_end[j])/2;
			int distance = abs(posj - posi);
			if (mode_peak) distance = j - i;
			if (distance <= max_dist && phenotype_end[i] < phenotype_start[j]) {
				double corr = getCorrelation(phenotype_val[i], phenotype_val[j]);
				double pval = getPvalue(corr, sample_count-2);

				vector < Interval < int > > overlaps;
				vector < bool > mask = vector < bool > (annotation_id_vec.size(), false);
				annotation_tree.findOverlapping(phenotype_end[i], phenotype_start[j], overlaps);
				for (int o = 0 ; o < overlaps.size() ; o ++) if (overlaps[o].start > phenotype_end[i] && overlaps[o].stop < phenotype_start[j]) mask[overlaps[o].value] = true;

				fd << i + 1 << " " << phenotype_chr[i] << " " << posi << " " << phenotype_id[i];
				fd << " " << j + 1 << " " << phenotype_chr[j] << " " << posj << " " << phenotype_id[j];
				fd << " " << stb.str(corr, 3) << " " << pval;
				for (int a = 0 ; a  < mask.size() ; a ++) fd << " " << mask[mapping[a].second];
				fd << endl;
			}
		}
		vrb.progress((i+1) * 1.0 / phenotype_count);
	}
	fd.close();
}

/*
void span_data::computeInteractionsPaired(string fout, int max_dist, bool mode_peak) {
	vrb.title("Report 2 all possible interactions in [" + fout + "]");
	output_file fd(fout.c_str());

	//
	vector < pair < string, int > > mapping;
	for ( int a = 0 ; a < annotation_id_vec.size() ; a ++) mapping.push_back(pair < string, int > (annotation_id_vec[a], a));
	sort(mapping.begin(), mapping.end());

	//
	fd << "IDX1 CHR1 POS1 ID1 IDX2 CHR2 POS2 ID2 CORR PVAL";
	for (int a0 = 0 ; a0 < annotation_id_vec.size() ; a0 ++) for (int a1 = a0 ; a1 < annotation_id_vec.size() ; a1 ++) fd << " " << mapping[a0].first << "_" << mapping[a1].first;
	fd << endl;

	//
	for (int i = 0 ; i < phenotype_count - 1 ; i ++) {
		for (int j = i + 1 ; j < phenotype_count ; j ++) {
			int posi =  (phenotype_start[i]+phenotype_end[i])/2;
			int posj =  (phenotype_start[j]+phenotype_end[j])/2;
			int distance = abs(posj - posi);
			if (mode_peak) distance = j - i;

			vector < Interval < int > > overlaps_i;
			vector < bool > mask_i = vector < bool > (annotation_id_vec.size(), false);
			annotation_tree.findOverlapping(phenotype_start[i], phenotype_end[i], overlaps_i);
			for (int o = 0 ; o < overlaps_i.size() ; o ++) mask_i[overlaps_i[o].value] = true;

			if (distance <= max_dist && phenotype_end[i] < phenotype_start[j] && overlaps_i.size() > 0) {
				double corr = getCorrelation(phenotype_val[i], phenotype_val[j]);
				double pval = getPvalue(corr, sample_count-2);

				vector < Interval < int > > overlaps_j;
				vector < bool > mask_j = vector < bool > (annotation_id_vec.size(), false);
				annotation_tree.findOverlapping(phenotype_start[j], phenotype_end[j], overlaps_j);
				for (int o = 0 ; o < overlaps_j.size() ; o ++) mask_j[overlaps_j[o].value] = true;

				if (overlaps_j.size() > 0) {
					fd << i + 1 << " " << phenotype_chr[i] << " " << posi << " " << phenotype_id[i];
					fd << " " << j + 1 << " " << phenotype_chr[j] << " " << posj << " " << phenotype_id[j];
					fd << " " << stb.str(corr, 3) << " " << pval;

					vector < vector < bool > > mask_ij = vector < vector < bool > > (annotation_id_vec.size(), vector < bool > (annotation_id_vec.size(), false));
					for (int ai = 0 ; ai  < mask_i.size() ; ai ++) for (int aj = 0 ; aj  < mask_j.size() ; aj ++) mask_ij[ai][aj] = mask_i[ai] && mask_j[aj];
					for (int a0 = 0 ; a0 < annotation_id_vec.size() ; a0 ++) for (int a1 = a0 ; a1 < annotation_id_vec.size() ; a1 ++) fd << " " << (mask_ij[mapping[a0].second][mapping[a1].second] || mask_ij[mapping[a1].second][mapping[a0].second]);
					fd << endl;
				}
			}
		}
		vrb.progress((i+1) * 1.0 / phenotype_count);
	}
	fd.close();
}
*/

void span_data::computeInteractionsPaired(string fout, int max_dist, bool mode_peak) {
	vrb.title("Compute means:");
	vector < double > me_obs = vector < double > (annotation_id_vec.size() * annotation_id_vec.size(), 0.0);
	vector < double > me_exp = vector < double > (annotation_id_vec.size() * annotation_id_vec.size(), 0.0);
	vector < double > sd_obs = vector < double > (annotation_id_vec.size() * annotation_id_vec.size(), 0.0);
	vector < double > sd_exp = vector < double > (annotation_id_vec.size() * annotation_id_vec.size(), 0.0);
	vector < int > nu_obs = vector < int > (annotation_id_vec.size() * annotation_id_vec.size(), 0);
	vector < int > nu_exp = vector < int > (annotation_id_vec.size() * annotation_id_vec.size(), 0);
	for (int i = 0 ; i < phenotype_count - 1 ; i ++) {
		for (int j = i + 1 ; j < phenotype_count ; j ++) {
			int posi =  (phenotype_start[i]+phenotype_end[i])/2;
			int posj =  (phenotype_start[j]+phenotype_end[j])/2;
			int distance = abs(posj - posi);
			if (mode_peak) distance = j - i;
			vector < Interval < int > > overlaps_i;
			vector < bool > mask_i = vector < bool > (annotation_id_vec.size(), false);
			annotation_tree.findOverlapping(phenotype_start[i], phenotype_end[i], overlaps_i);
			for (int o = 0 ; o < overlaps_i.size() ; o ++) mask_i[overlaps_i[o].value] = true;
			if (distance <= max_dist && phenotype_end[i] < phenotype_start[j] && overlaps_i.size() > 0) {
				double corr = getCorrelation(phenotype_val[i], phenotype_val[j]);
				vector < Interval < int > > overlaps_j;
				vector < bool > mask_j = vector < bool > (annotation_id_vec.size(), false);
				annotation_tree.findOverlapping(phenotype_start[j], phenotype_end[j], overlaps_j);
				for (int o = 0 ; o < overlaps_j.size() ; o ++) mask_j[overlaps_j[o].value] = true;
				if (overlaps_j.size() > 0) {
					for (int oi = 0 ; oi < mask_i.size() ; oi ++) {
						for (int oj = 0 ; oj < mask_j.size() ; oj ++) {
							int idx = oi * mask_i.size() + oj;
							if (mask_i[oi] && mask_j[oj]) {
								me_obs[idx] += corr;
								nu_obs[idx] ++;
							} else {
								me_exp[idx] += corr;
								nu_exp[idx] ++;
							}
						}
					}
				}
			}
		}
		vrb.progress((i+1) * 1.0 / phenotype_count);
	}

	vrb.title("Compute sds:");
	for (int i = 0 ; i < phenotype_count - 1 ; i ++) {
		for (int j = i + 1 ; j < phenotype_count ; j ++) {
			int posi =  (phenotype_start[i]+phenotype_end[i])/2;
			int posj =  (phenotype_start[j]+phenotype_end[j])/2;
			int distance = abs(posj - posi);
			if (mode_peak) distance = j - i;
			vector < Interval < int > > overlaps_i;
			vector < bool > mask_i = vector < bool > (annotation_id_vec.size(), false);
			annotation_tree.findOverlapping(phenotype_start[i], phenotype_end[i], overlaps_i);
			for (int o = 0 ; o < overlaps_i.size() ; o ++) mask_i[overlaps_i[o].value] = true;
			if (distance <= max_dist && phenotype_end[i] < phenotype_start[j] && overlaps_i.size() > 0) {
				double corr = getCorrelation(phenotype_val[i], phenotype_val[j]);
				vector < Interval < int > > overlaps_j;
				vector < bool > mask_j = vector < bool > (annotation_id_vec.size(), false);
				annotation_tree.findOverlapping(phenotype_start[j], phenotype_end[j], overlaps_j);
				for (int o = 0 ; o < overlaps_j.size() ; o ++) mask_j[overlaps_j[o].value] = true;
				if (overlaps_j.size() > 0) {
					for (int oi = 0 ; oi < mask_i.size() ; oi ++) {
						for (int oj = 0 ; oj < mask_j.size() ; oj ++) {
							int idx = oi * mask_i.size() + oj;
							if (mask_i[oi] && mask_j[oj]) {
								sd_obs[idx] += (corr - me_obs[idx]/nu_obs[idx]) * (corr - me_obs[idx]/nu_obs[idx]);
							} else {
								sd_exp[idx] += (corr - me_exp[idx]/nu_exp[idx]) * (corr - me_exp[idx]/nu_exp[idx]);
							}
						}
					}
				}
			}
		}
		vrb.progress((i+1) * 1.0 / phenotype_count);
	}

	//
	vector < pair < string, int > > mapping;
	for ( int a = 0 ; a < annotation_id_vec.size() ; a ++) mapping.push_back(pair < string, int > (annotation_id_vec[a], a));
	sort(mapping.begin(), mapping.end());

	//
	output_file fd(fout.c_str());
	for (int a0 = 0 ; a0 < annotation_id_vec.size() ; a0 ++) for (int a1 = 0 ; a1 < annotation_id_vec.size() ; a1 ++) {
		int idx = mapping[a0].second * annotation_id_vec.size() + mapping[a1].second;
		fd << mapping[a0].first << " " << mapping[a1].first << " " << me_obs[idx] << " " << me_exp[idx] << " " << sd_obs[idx] << " " << sd_exp[idx] << " " << nu_obs[idx] << " " << nu_exp[idx] << endl;
	}
	fd.close();
}
