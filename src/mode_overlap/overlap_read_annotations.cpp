#include "overlap_data.h"

void overlap_data::readAnnotations(string fbed) {
	string buffer;
	vector < string > tokens;

	//Open BED file
	vrb.title("Reading annotation data in [" + fbed + "]");
	annotation_count = 0;
	input_file fd (fbed);
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() < 4) vrb.error("Incorrect number of column somewhere in the file!");

		int idx_chr = -1;
		for (int e = 0 ; e < annotation_unique_chr.size() && idx_chr < 0 ; e++) if (annotation_unique_chr[e] == tokens[0]) idx_chr = e;
		if (idx_chr == -1) {
			annotation_chr.push_back(annotation_unique_chr.size());
			annotation_unique_chr.push_back(tokens[0]);
		} else annotation_chr.push_back(idx_chr);

		int idx_type = -1;
		for (int e = 0 ; e < annotation_unique_type.size() && idx_type < 0 ; e++) if (annotation_unique_type[e] == tokens[3]) idx_type = e;
		if (idx_type == -1) {
			annotation_type.push_back(annotation_unique_type.size());
			annotation_unique_type.push_back(tokens[3]);
		} else annotation_type.push_back(idx_type);

		annotation_start.push_back(atoi(tokens[1].c_str()));
		annotation_end.push_back(atoi(tokens[2].c_str()));
		annotation_count++;
	}
	vrb.bullet(stb.str(annotation_count) + " annotations read");

	vrb.title("List of chromosomes discovered:");
	for (int e = 0 ; e < annotation_unique_chr.size() ; e++) vrb.bullet(annotation_unique_chr[e]);

}
