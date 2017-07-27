#include "binding_data.h"

void binding_data::readReference(string file) {
	string buffer, idseq = "", valseq = "";
	input_file fd(file);

	vrb.title("Read reference genome in [" + file + "]");
	while(getline(fd, buffer)) {
		if (buffer[0] == '>') {
			if (idseq != "") refseq.insert(pair < string, string > (idseq, valseq));
			idseq = buffer.substr(1);
			valseq = "";
		} else valseq += buffer;
	}
	refseq.insert(pair < string, string > (idseq, valseq));
	fd.close();
	vrb.bullet("#seq = " + stb.str(refseq.size()));
}

void binding_data::readPosition(string file) {
	string buffer;
	vector < string > tokens;
	input_file fd(file);

	vrb.title("Read variant positions in [" + file + "]");
	while(getline(fd, buffer)) {
		stb.split(buffer, tokens);
		chr.push_back(tokens[0]);
		pos.push_back(atoi(tokens[1].c_str()));
		ref.push_back(tokens[2]);
		alt.push_back(tokens[3]);
	}
	fd.close();
	vrb.bullet("#pos = " + stb.str(chr.size()));

}
