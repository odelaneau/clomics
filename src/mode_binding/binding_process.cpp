#include "binding_data.h"

void binding_data::process() {
	vrb.title("Check concordance between alleles and reference genome");
	for (int p = 0 ; p < chr.size() ; p ++) {
		map < string, string > :: iterator it = refseq.find(chr[p]);
		string refseq_content = it->second.substr(pos[p] - 1, ref[p].size());
		cout << ref[p] << " " << refseq_content << endl;
	}
}

