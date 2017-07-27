#include "corr_data.h"

void corr_data::pairingPhenotypes(string file) {
	string buffer;
	vector < string > tokens;

	//
	vrb.title("Reading phenotype pairing data in [" + file + "]");
	map < string, map < string, bool > > pairs;
	int n_line = 0;
	input_file fd(file);
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		map < string, map < string, bool > > :: iterator it = pairs.find (tokens[1]);
		if (it == pairs.end()) {
			map < string, bool > s;
			s.insert(pair < string, bool > (tokens[0], tokens[2] == "1"));
			pairs.insert(pair < string , map < string, bool > > (tokens[1], s));
		} else it->second.insert(pair < string, bool > (tokens[0], tokens[2] == "1"));
		n_line ++;
	}
	vrb.bullet("#pairs=" + stb.str(n_line));

	//
	vrb.title("Processing pairs");
	int n_pair_set = 0;
	phenotype_paired = vector < vector < int > > (phenotype_count, vector < int > (phenotype_count, -1));
	for (map < string, map < string, bool > > :: iterator itM = pairs.begin() ; itM != pairs.end() ; itM ++) {
		for (map < string, bool > :: iterator itS1 = itM->second.begin() ; itS1 != itM->second.end() ; itS1 ++) {
			for (map < string, bool > :: iterator itS2 = itM->second.begin() ; itS2 != itM->second.end() ; itS2 ++) {
				unordered_map < string, unsigned int >::iterator itP1 = phenotype_map.find(itS1->first);
				unordered_map < string, unsigned int >::iterator itP2 = phenotype_map.find(itS2->first);
				if (itP1 != phenotype_map.end() && itP2 != phenotype_map.end()) {
					phenotype_paired[itP1->second][itP2->second] = (itS1->second == itS2->second);
					n_pair_set ++;
				}
			}
		}
	}
	vrb.bullet("#pairs_set=" + stb.str(n_pair_set));
}

