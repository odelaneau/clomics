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

#include "call_data.h"

void call_main(vector < string > & argv) {
	call_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("tree", boost::program_options::value< string >(), "Input tree file")
		("out", boost::program_options::value< string >(), "Output tree file");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("threshold", boost::program_options::value< float >()->default_value(3.0), "Pvalue threshold");

	D.option_descriptions.add(opt_files).add(opt_parameters);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [call] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("CALL MODULES");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}
	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("tree")) vrb.error("Specify a tree with --tree [file.txt]");
	if (!D.options.count("out")) vrb.error("Specify output file with --out [file.out]");

	//---------------
	// 6. READ FILES
	//---------------
	vector < int > dummy;
	D.processBasicOptions();
	D.readTree(D.options["tree"].as < string > ());
    D.tagSignificant(D.tree, 0, D.options["threshold"].as < float > (), D.options["threshold"].as < float > (), dummy);

    //-----------------
	// 10. WRITE OUTPUT
	//-----------------
	D.writeTree(D.options["out"].as < string > ());
}
