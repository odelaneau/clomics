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

#include "corr_data.h"

void corr_main(vector < string > & argv) {
	corr_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("bed", boost::program_options::value< vector < string > >()->multitoken(), "List of phenotype files in BED format.")
		("out", boost::program_options::value< string >(), "Output file for constructed tree.")
		("pair", boost::program_options::value< string >(), "Phenotype pairing");

	boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
	opt_parallel.add_options()
		("region", boost::program_options::value< string >(), "Region of interest.");

	D.option_descriptions.add(opt_files).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [corr] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("CLUSTERIZE MOLECULAR PHENOTYPE DATA");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}
	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("bed")) vrb.error("Specify phenotypes with --bed [file.bed]");
	if (!D.options.count("out")) vrb.error("Specify output tree with --out [file.out]");

	//--------------
	// 5. SET REGION
	//--------------
	if (D.options.count("region") && !D.setPhenotypeRegion(D.options["region"].as < string > ()))
		vrb.error("Impossible to interpret region [" + D.options["region"].as < string > () + "]");

	//---------------
	// 6. READ FILES
	//---------------
	vector < string > bed_list = D.options["bed"].as < vector < string > > ();
	D.processBasicOptions();
	for (int b = 0 ; b < bed_list.size() ; b ++) D.readSampleFromBED(bed_list[b]);
	D.mergeSampleLists();
	for (int b = 0 ; b < bed_list.size() ; b ++) D.scanPhenotypes(bed_list[b]);
	D.sortPhenotypes();
	for (int b = 0 ; b < bed_list.size() ; b ++) D.readPhenotypes(bed_list[b]);

	//-----------------------
	// 7. INITIALIZE ANALYSIS
	//-----------------------
	D.imputePhenotypes();
	D.normalizePhenotypes();
	if (D.options.count("pair")) D.pairingPhenotypes(D.options["pair"].as < string > ());

	//----------------
	// 8. RUN ANALYSIS
	//----------------
	D.computeCorrelations(D.options["out"].as < string > ());
}
