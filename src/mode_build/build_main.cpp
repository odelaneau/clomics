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

#define _DECLARE_TOOLBOX_HERE
#include "build_data.h"

void build_main(vector < string > & argv) {
	build_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("bed", boost::program_options::value< vector < string > >()->multitoken(), "List of phenotype files in BED format.")
		("out", boost::program_options::value< string >(), "Output file for constructed tree.");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("bootstrap", boost::program_options::value< float >(), "Bootstrap re-sampling X% of data")
		("jackknife", boost::program_options::value< float >(), "Jackknife re-sampling X% of data");

	boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
	opt_parallel.add_options()
		("region", boost::program_options::value< string >(), "Region of interest.");

	D.option_descriptions.add(opt_files).add(opt_parameters).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [build] command line :" << string(e.what()) << endl;
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
	if (D.options.count("bootstrap") && D.options.count("jackknife")) vrb.error("Specify only one of --bootstrap and --jackknife");
	if (D.options.count("jackknife")) {
		float prop = D.options["jackknife"].as < float > ();
		if (prop < 0.1) vrb.error("When using --jackknife X, X needs to be greater than 0.1");
		if (prop >= 1.0) vrb.error("When using --jackknife X, X needs to be lower than 1.0");
		vrb.bullet("Jackknife re-sample " + stb.str(prop, 1) + "% of the phenotype data");
	}
	if (D.options.count("bootstrap")) {
		float prop = D.options["bootstrap"].as < float > ();
		if (prop < 0.1) vrb.error("When using --bootstrap X, X needs to be greater than 0.1");
		if (prop > 1.0) vrb.error("When using --bootstrap X, X needs to be lower or equal than 1.0");
		vrb.bullet("Bootstrap re-sample " + stb.str(prop, 1) + "% of the phenotype data");
	}

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

	if (D.options.count("bootstrap")) {
		D.bootstrapPhenotypes(D.options["bootstrap"].as < float > ());
		D.resampling = true;
	}
	if (D.options.count("jackknife")) {
		D.jackknifePhenotypes(D.options["jackknife"].as < float > ());
		D.resampling = true;
	}

	//-----------------------
	// 7. INITIALIZE ANALYSIS
	//-----------------------
	D.imputePhenotypes();
	D.normalizePhenotypes();

	//----------------
	// 8. RUN ANALYSIS
	//----------------
	D.clusterize();

	//----------------
	// 9. WRITE OUTPUT
	//----------------
	D.writeTree(D.options["out"].as < string > ());
}
