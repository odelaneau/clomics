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

#include "score_data.h"

void score_main(vector < string > & argv) {
	score_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("spl", boost::program_options::value< vector < string > >()->multitoken(), "Bootstrapped/Jack-knife trees")
		("nom", boost::program_options::value< string >(), "Nominal tree")
		("out", boost::program_options::value< string >(), "Nominal tree populated with quality scores.");

	D.option_descriptions.add(opt_files);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [score] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("SCORE CLUSTERIZATION OF MOLECULAR PHENOTYPE DATA");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("nom")) vrb.error("Please provide a nominal tree with --nom [file.txt]");
	if (!D.options.count("spl")) vrb.error("Please provide resampled trees with --spl [file1.txt ... fileN.txt]");
	if (!D.options.count("out")) vrb.error("Please provide an output tree with --out [file.txt]");

	//---------------
	// 9. READ FILES
	//---------------
	D.processBasicOptions();
	vector < string > tree_list = D.options["spl"].as < vector < string > > ();
	D.readNominalTree(D.options["nom"].as < string > ());
	D.readSampledTrees(tree_list);


	//-----------------
	// 11. RUN ANALYSIS
	//-----------------
	D.computeOverlapScore();
	D.computeJaccardScore();

	//-----------------
	// 12. WRITE OUTPUT
	//-----------------
	D.writeTree(D.options["out"].as < string > ());
}
