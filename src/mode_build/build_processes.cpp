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

#include "build_data.h"

void build_data::clusterize() {

	//STEP1: compute all pairwise correlations vidx = i * phenotype_count - (i + 1) * i / 2 + j - i - 1
	vrb.title("Compute all n=" + stb.str(phenotype_count * (phenotype_count - 1) / 2) + " possible correlations");
	vector < vector < float > > Kdist = vector <  vector < float > > (phenotype_count, vector < float > (phenotype_count, 0.0));
	for (int i = 0 ; i < phenotype_count - 1 ; i ++) {
		for (int j = i + 1 ; j < phenotype_count ; j ++) Kdist[i][j] = getCorrelation(phenotype_val[i], phenotype_val[j]);
		vrb.progress((i+1) * 1.0 / (phenotype_count - 1));
	}

	//STEP2: initialize vector of nodes
	int k = 0;
	vector < int > Kidx = vector < int > (phenotype_count, 0);
	vector < bcluster * > Knode = vector < bcluster * > (phenotype_count, NULL);
	vector < vector < int > > Kcontent = vector < vector < int > > (phenotype_count, vector < int > (1, -1));
	for ( ; k < phenotype_count ; k ++) {
		Knode[k] = new bcluster(phenotype_id[k], k, 1);
		Kcontent[k][0] = k;
		Kidx[k] = k;
	}

	//STEP3:
	vrb.title("Clusterize phenotypes by aggregation");
	unsigned int n_iteration = 0;
	vector < float > median_vector = vector < float > (500);
	while (Knode.size() > 1) {

		//Find max
		int max_i = -1, max_j = -1;
		float max_val = 0.0;
		for (int i = 0 ; i < Kidx.size() - 1 ; i ++) for (int j = i + 1 ; j < i + 200 && j < Kidx.size() ; j ++) {
			if (abs(Kdist[Kidx[i]][Kidx[j]]) >= max_val) {
				max_val = abs(Kdist[Kidx[i]][Kidx[j]]);
				max_i = i;
				max_j = j;
			}
		}

		//Update correlation matrix for 0 <= i < i_max
		int n_node_max_i = Knode[max_i]->count;
		int n_node_max_j = Knode[max_j]->count;
		for (int i = 0 ; i < Kidx.size() ; i ++) {
			if (i != max_i && i != max_j) {
				if (i < max_i && i < max_j) Kdist[Kidx[i]][Kidx[max_i]] = (Kdist[Kidx[i]][Kidx[max_i]] * n_node_max_i + Kdist[Kidx[i]][Kidx[max_j]] * n_node_max_j) / (n_node_max_i + n_node_max_j);
				if (i > max_i && i < max_j) Kdist[Kidx[max_i]][Kidx[i]] = (Kdist[Kidx[max_i]][Kidx[i]] * n_node_max_i + Kdist[Kidx[i]][Kidx[max_j]] * n_node_max_j) / (n_node_max_i + n_node_max_j);
				if (i > max_i && i > max_j) Kdist[Kidx[max_i]][Kidx[i]] = (Kdist[Kidx[max_i]][Kidx[i]] * n_node_max_i + Kdist[Kidx[max_j]][Kidx[i]] * n_node_max_j) / (n_node_max_i + n_node_max_j);
			}
		}

		//Merge nearest clusters
		Knode[max_i] = new bcluster(string(phenotype_chr[0] + "_internal_" + stb.str(k+1)), k, Knode[max_i]->count + Knode[max_j]->count, Knode[max_i], Knode[max_j]);
		Knode.erase(Knode.begin() + max_j);
		Kcontent[max_i].insert(Kcontent[max_i].end(), Kcontent[max_j].begin(), Kcontent[max_j].end());
		Kcontent.erase(Kcontent.begin() + max_j);
		Kidx.erase(Kidx.begin() + max_j);
		k++;


		if (!resampling) {

			//PCA
			pca P(sample_count, Kcontent[max_i].size());
			P.fill(phenotype_val, Kcontent[max_i]);
			P.run(false, false, false);
			Knode[max_i]->pcavariance = P.getVariance(0);

			//APVALUE (aggregate)
			int n_acells = 0;
			Knode[max_i]->acorrelation = 0.0;
			for (int i = 0 ; i < Kcontent[max_i].size() - 1 ; i ++)
				for (int j = i+1 ; j < Kcontent[max_i].size(); j ++) {
					Knode[max_i]->acorrelation += abs(Kdist[Kcontent[max_i][i]][Kcontent[max_i][j]]);
					n_acells ++;
				}
			Knode[max_i]->acorrelation /= n_acells;
			Knode[max_i]->apvalue = getPvalue(Knode[max_i]->acorrelation, sample_count - 2);

			//BPVALUE (branching)
			Knode[max_i]->bcorrelation = max_val;
			Knode[max_i]->bpvalue = getPvalue(max_val, sample_count - 2);

			//CPVALUE & DPVALUE (cumulative & TOP corner)
			int smallest_idx = 1000000, largest_idx = 0, n_ccells = 0, n_dcells = 0;
			for (int i = 0 ; i < Kcontent[max_i].size(); i ++) {
				if (Kcontent[max_i][i] > largest_idx) largest_idx = Kcontent[max_i][i];
				if (Kcontent[max_i][i] < smallest_idx) smallest_idx = Kcontent[max_i][i];
			}
			Knode[max_i]->ccorrelation = 0.0;
			Knode[max_i]->dcorrelation = 0.0;
			for (int i = smallest_idx ; i < largest_idx ; i ++) {
				for (int j = i+1 ; j <= largest_idx ; j ++) {
					if (i < (smallest_idx + 5) && j > (largest_idx - 5)) {
						Knode[max_i]->dcorrelation += abs(Kdist[i][j]);
						n_dcells ++;
					}
					Knode[max_i]->ccorrelation += abs(Kdist[i][j]);
					n_ccells ++;
				}
			}
			Knode[max_i]->ccorrelation /= n_ccells;
			Knode[max_i]->dcorrelation /= n_dcells;
			Knode[max_i]->cpvalue = getPvalue(Knode[max_i]->ccorrelation, sample_count - 2);
			Knode[max_i]->dpvalue = getPvalue(Knode[max_i]->dcorrelation, sample_count - 2);

			//E CORRELATION
			int n_ecells = 0;
			Knode[max_i]->ecorrelation = 0.0;
			for (int i = smallest_idx+1 ; i <= largest_idx ; i ++) { Knode[max_i]->ecorrelation += abs(Kdist[smallest_idx][i]); n_ecells ++; }
			for (int i = smallest_idx+1 ; i < largest_idx ; i ++) { Knode[max_i]->ecorrelation += abs(Kdist[i][largest_idx]); n_ecells ++; }
			Knode[max_i]->ecorrelation /= n_ecells;
			Knode[max_i]->epvalue = getPvalue(Knode[max_i]->ecorrelation, sample_count - 2);
			/*
			int n_ecells = 0;
			double prob = 500.0 / n_ccells;
			for (int i = smallest_idx ; i < largest_idx && n_ecells < 500; i ++) for (int j = i+1 ; j <= largest_idx  && n_ecells < 500; j ++) {
				if (rng.getDouble() <= prob) {
					median_vector[n_ecells] = abs(Kdist[i][j]);
					n_ecells ++;
				}
			}
			sort(median_vector.begin(), median_vector.begin() + n_ecells);
			Knode[max_i]->ecorrelation = median_vector[n_ecells/2];
			Knode[max_i]->epvalue = getPvalue(Knode[max_i]->ecorrelation, sample_count - 2);
			*/
		}

		//Update iterator
		n_iteration ++;

		//Verbose
		vrb.progress((phenotype_count - Knode.size()) * 1.0 / (phenotype_count - 1));
		//vrb.bullet("Iteration i=" + stb.str(n_iteration) + "\tn1=" + stb.str(max_i) + " (s=" + stb.str(n_node_max_i) + ")\tn2=" + stb.str(max_j) + " (s=" + stb.str(n_node_max_j) + ")\t corr=" + stb.str(max_val, 4));
	}

	tree = Knode[0];
}
