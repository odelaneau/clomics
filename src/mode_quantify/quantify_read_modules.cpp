#include "quantify_data.h"

void quantify_data::readModules(string ftxt) {
	string buffer;
	vector < string > tokens;

	vrb.title("Reading modules in [" + ftxt + "]");

	//Sweep1: reading nodes
	input_file fd(ftxt);
	while (getline(fd, buffer)) {
		modules.insert(buffer);
	}
	vrb.bullet("#modules = " + stb.str(modules.size()));
	fd.close();

}
