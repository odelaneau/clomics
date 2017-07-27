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

#include "motif_trans_data.h"

void motif_trans_main(vector < string > & argv) {
	motif_trans_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("bed", boost::program_options::value< vector < string > >()->multitoken(), "List of phenotype files in BED format")
		("ann", boost::program_options::value< string >(), "Annotations to look at")
		("out", boost::program_options::value< string >(), "Output file");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("permute", "Permute data to draw from null");

	boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
	opt_parallel.add_options()
		("regions", boost::program_options::value< vector < string > >()->multitoken(), "Region of interest.");

	D.option_descriptions.add(opt_files).add(opt_parameters).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [motif-trans] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("CHECK FOR PAIRWISE ENRICHMENT OF ANNOTATION DATA INTO INTERACTION DATA IN TRANS");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}
	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("bed")) vrb.error("Specify phenotypes with --bed [file.bed]");
	if (!D.options.count("ann")) vrb.error("Specify motives with --ann [file.bed]");
	if (!D.options.count("out")) vrb.error("Specify output tree with --out [file.out]");
	if (!D.options.count("regions")) vrb.error("Specify two regions to be processed with --regions chr12 chr14]");

	//--------------
	// 5. SET REGION
	//--------------
	vector < string > region_list = D.options["regions"].as < vector < string > > ();
	if (!D.setPhenotypeRegion1(region_list[0])) vrb.error("Impossible to interpret region1 [" + region_list[0] + "]");
	if (!D.setPhenotypeRegion2(region_list[1])) vrb.error("Impossible to interpret region2 [" + region_list[1] + "]");

	//---------------
	// 6. READ FILES
	//---------------
	vector < string > bed_list = D.options["bed"].as < vector < string > > ();
	D.processBasicOptions();
	for (int b = 0 ; b < bed_list.size() ; b ++) D.readSampleFromBED(bed_list[b]);
	D.mergeSampleLists();
	for (int b = 0 ; b < bed_list.size() ; b ++) {
		D.scanPhenotypes1(bed_list[b]);
		D.scanPhenotypes2(bed_list[b]);
	}
	D.sortPhenotypes1();
	D.sortPhenotypes2();
	for (int b = 0 ; b < bed_list.size() ; b ++) {
		D.readPhenotypes1(bed_list[b]);
		D.readPhenotypes2(bed_list[b]);
	}
	D.buildIntervalTree1();
	D.buildIntervalTree2();

	//-----------------------
	// 7. INITIALIZE ANALYSIS
	//-----------------------
	D.imputePhenotypes1();
	D.imputePhenotypes2();
	D.normalizePhenotypes1();
	D.normalizePhenotypes2();

	//-----------
	// 9. MOTIVES
	//-----------
	D.readMotives1(D.options["ann"].as < string > (), D.options.count("permute"));
	D.readMotives2(D.options["ann"].as < string > (), D.options.count("permute"));
	D.initializeMotives();
	D.mappingMotives1();
	D.mappingMotives2();

	//-------------
	// 10. RUN ANAL
	//-------------
	D.computeCorrelations();
	D.writePairwiseData(D.options["out"].as < string > ());
}
