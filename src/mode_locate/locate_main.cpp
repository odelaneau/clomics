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

#include "locate_data.h"

int regulatory_element::merge_distance = DEFAULT_MERGE_DISTANCE;

void locate_main(vector < string > & argv) {
	locate_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("bed", boost::program_options::value< vector < string > >()->multitoken(), "List of phenotype files in BED format.")
		("tree", boost::program_options::value< string >(), "Tree file.")
		("out", boost::program_options::value< string >(), "Output file for tree.");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("distance", boost::program_options::value< int >()->default_value(DEFAULT_MERGE_DISTANCE), "Two peaks with this distance in between are merged!");

	D.option_descriptions.add(opt_files).add(opt_parameters);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [locate] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("SET GENOMIC POSITIONS OF TREE NODES");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}
	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("bed")) vrb.error("Specify phenotypes with --bed [file.bed]");
	if (!D.options.count("out")) vrb.error("Specify output tree with --out [file.out]");
	if (!D.options.count("tree")) vrb.error("Specify a tree with --tree [file.txt]");

	vrb.bullet("Merging distance between peaks = " + stb.str(D.options["distance"].as < int > ()));
	regulatory_element::merge_distance = D.options["distance"].as < int > ();

	//---------------
	// 6. READ FILES
	//---------------
	D.processBasicOptions();
	D.readTree(D.options["tree"].as < string > ());
	vector < string > bed_list = D.options["bed"].as < vector < string > > ();
	D.readPhenotypes(bed_list);

	//----------------
	// 8. RUN ANALYSIS
	//----------------
	D.propagateRegulatoryElements();

	//----------------
	// 9. WRITE OUTPUT
	//----------------
	D.writeTree(D.options["out"].as < string > ());
}
