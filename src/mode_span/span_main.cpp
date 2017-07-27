#include "span_data.h"

void span_main(vector < string > & argv) {
	span_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("bed", boost::program_options::value< vector < string > >()->multitoken(), "List of phenotype files in BED format.")
		("ann", boost::program_options::value< string >(), "Annotation file in BED format.")
		("out", boost::program_options::value< string >(), "Output file for interactions.");

	boost::program_options::options_description opt_param ("\x1B[32mParameters\33[0m");
	opt_param.add_options()
		("window", boost::program_options::value< int >()->default_value(1000000), "Cis-window to be explored")
		("distance-by-peaks", "Cis-window size applies to distance in number of peaks and not physical distance")
		("test-pairs", "Cis-window size applies to distance in number of peaks and not physical distance");

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
		cerr << "Error parsing [span] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("OVERLAP INTERACTIONS WITH ANNOTATIONS");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("bed")) vrb.error("Specify phenotypes with --bed [files.bed]");
	if (!D.options.count("ann")) vrb.error("Specify annotation data with --ann [files.bed]");
	if (!D.options.count("region")) vrb.error("Specify region with --region [chr1]");
	if (!D.options.count("out")) vrb.error("Specify output tree with --out [file.out]");

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

	D.readAnnotations(D.options["ann"].as < string > ());

	//----------------
	// 8. RUN ANALYSIS
	//----------------
	if (D.options.count("test-pairs")) D.computeInteractionsPaired(D.options["out"].as < string > (), D.options["window"].as < int > (), D.options.count("distance-by-peaks"));
	else D.computeInteractionsSingle(D.options["out"].as < string > (), D.options["window"].as < int > (), D.options.count("distance-by-peaks"));
}
