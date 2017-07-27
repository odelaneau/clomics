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

#include "hic_data.h"

void hic_main(vector < string > & argv) {
	hic_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("bed", boost::program_options::value< vector < string > >()->multitoken(), "List of phenotype files in BED format.")
		("hic", boost::program_options::value< vector < string > >()->multitoken(), "Hi-C interaction and normalization files.")
		("out", boost::program_options::value< string >(), "Output file for interactions.");

	boost::program_options::options_description opt_param ("\x1B[32mParameters\33[0m");
	opt_param.add_options()
		("null", boost::program_options::value< double >()->default_value(0.01), "Percentage null to be reported in the output");

	boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
	opt_parallel.add_options()
		("region", boost::program_options::value< string >(), "Region of interest.");

	D.option_descriptions.add(opt_files).add(opt_param).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [hic] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("PROJECT HIC DATA ONTO DEFINED INTERVALS");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("bed")) vrb.error("Specify phenotypes with --bed [file.bed]");
	if (!D.options.count("hic")) vrb.error("Specify hic data with --hic [file.mat file1.norm file2.norm bin_size]");
	if (!D.options.count("region")) vrb.error("Specify region with --region [chr1]");
	if (!D.options.count("out")) vrb.error("Specify output tree with --out [file.out]");

	vector < string > trans_list = D.options["hic"].as < vector < string > > ();
	if (trans_list.size() != 3) vrb.error("--hic take 3 arguments: interation_matrix normalization_vector bin_size");
	D.hic_bin = atoi(trans_list[2].c_str());
	vrb.bullet("Bin size of Hi-C data = " + stb.str(D.hic_bin));

	D.thinin = D.options["null"].as < double > ();
	vrb.bullet("Percentage of null interactions to be reported = " + stb.str(D.thinin*100, 1));

	//--------------
	// 5. SET REGION
	//--------------
	if (!D.setPhenotypeRegion(D.options["region"].as < string > ())) vrb.error("Impossible to interpret region [" + D.options["region"].as < string > () + "]");

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
	D.normalizePhenotypes();

	D.readHiCnorms(trans_list[1]);
	D.readHiCvalues(trans_list[0]);

	//----------------
	// 8. RUN ANALYSIS
	//----------------
	D.computeInteractions(D.options["out"].as < string > ());
}
