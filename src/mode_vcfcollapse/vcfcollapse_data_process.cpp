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

void vcfcollapse_data::normalTransform(vector < float > & V) {
	vector < float > R;
	myranker::rank(V, R);
	double max = 0;
	for (int s = 0 ; s < sample_count ; s ++) {
		R[s] = R[s] - 0.5;
		if (R[s] > max) max = R[s];
	}
	max = max + 0.5;
	for (int s = 0 ; s < sample_count ; s ++) {
		R[s] /= max;
		V[s] = qnorm(R[s], 0.0, 1.0, 1, 0);
	}
}


void vcfcollapse_data::buildTrees() {
	vrb.title("Building interval trees");
	vector < vector < Interval < int > > > Tvec;
	for (int a = 0 ; a < ann_count ; a ++) {
		map < string, int >::iterator itC = ann_chr_map.find(ann_chr[a]);
		if (itC == ann_chr_map.end()) {
			ann_chr_map.insert(pair < string, int > (ann_chr[a], ann_chr_map.size()));
			Tvec.push_back(vector < Interval < int > > (1, Interval < int > (ann_start[a], ann_end[a], a)));
		} else Tvec[itC->second].push_back(Interval < int > (ann_start[a], ann_end[a], a));
	}
	ann_int = vector <  IntervalTree < int > > (ann_chr_map.size(), IntervalTree < int > ());
	for (int c = 0 ; c < ann_chr_map.size() ; c ++) ann_int[c] = IntervalTree < int > (Tvec[c]);
	vrb.bullet("#detected chromosomes in annotation data = " + stb.str(ann_chr_map.size()));
}


void vcfcollapse_data::collapseGenotypes(bool normalize) {
	vrb.title("Collapsing genotypes");

	int n_kept = 0;
	vector < bool > flag_maf = vector < bool > (genotype_count, false);
	vector < bool > flag_flip = vector < bool > (genotype_count, false);
	for (int g = 0 ; g < genotype_count ; g ++) {
		float c = 0.0; float nm = 0.0;
		for (int i = 0; i < sample_count ; i ++) if (genotype_val[g][i] != bcf_float_missing) {
			c += genotype_val[g][i];
			nm ++;
		}
		float freq = c / (2 * nm);
		if (freq > 0.5) {
			freq = 1- freq;
			flag_flip[g] = true;
		}
		if (freq > 0 && freq < param_max_maf) {
			flag_maf[g] = true;
			n_kept ++;
		}
	}
	vrb.bullet("#variants passing MAF filter = " + stb.str(n_kept));

	map < string, int > bed2idx;
	for (int g = 0 ; g < genotype_count ; g++) {
		if (flag_maf[g]) {
			map < string, int > :: iterator itB = bed2idx.find(genotype_bed[g]);
			if (itB == bed2idx.end()) bed2idx.insert(pair < string, int > (genotype_bed[g], bed2idx.size()));
		}
	}
	vrb.bullet("#collapsed variants = " + stb.str(bed2idx.size()));

	burden_count = bed2idx.size();
	burden_val = vector < vector < float > > (burden_count, vector < float > (sample_count, 0));
	burden_chr = vector < string > (burden_count);
	burden_id = vector < string > (burden_count);
	burden_start = vector < int > (burden_count, 1000000000);
	burden_end = vector < int > (burden_count, 0);


	int n_flip = 0;
	for (int g = 0 ; g < genotype_count ; g++) {
		if (flag_maf[g]) {
			map < string, int > :: iterator itB = bed2idx.find(genotype_bed[g]);
			burden_chr[itB->second] = genotype_chr[g];
			if (burden_start[itB->second] > genotype_start[g]) burden_start[itB->second] = genotype_start[g];
			if (burden_end[itB->second] < genotype_end[g]) burden_end[itB->second] = genotype_end[g];
			burden_id[itB->second] = itB->first;

			if (flag_flip[g]) n_flip ++;

			for (int i = 0 ; i < sample_count ; i ++) {
				if (flag_flip[g]) burden_val[itB->second][i] += 2 - genotype_val[g][i];
				else burden_val[itB->second][i] += genotype_val[g][i];
			}

			if (normalize) normalTransform(burden_val[itB->second]);
		}
	}
	vrb.bullet("#flipped variants = " + stb.str(n_flip));
}
