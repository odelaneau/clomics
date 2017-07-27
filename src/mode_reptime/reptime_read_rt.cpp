#include "reptime_data.h"

void reptime_data::readReptimes(string frt) {
	string buffer;
	vector < string > tokens;

	vrb.title("Reading Replication times in [" + frt + "]");
	input_file fd (frt);
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);

		if (tokens.size() != 3) vrb.error("Incorrect number of column somewhere in the file!");

		if (tokens[0] == phenotype_chr[0]) {
			reptime_pos.push_back(atoi(tokens[1].c_str()));
			reptime_val.push_back(atof(tokens[2].c_str()));
		}
	}

	reptime_count = reptime_val.size();
	vrb.bullet(stb.str(reptime_count) + " times read");
}
