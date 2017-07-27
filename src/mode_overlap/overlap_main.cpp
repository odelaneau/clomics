#include "overlap_data.h"

void overlap_main(vector < string > & argv) {
	overlap_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("bed", boost::program_options::value< string >(), "List of phenotype files in BED format.")
		("ann", boost::program_options::value< string >(), "List of annotation file in BED format.")
		("out", boost::program_options::value< string >(), "Output file overlap proportions.");

	boost::program_options::options_description opt_param ("\x1B[32mParameters\33[0m");
	opt_param.add_options()
		("permute", "Permute annotation types");

	D.option_descriptions.add(opt_files).add(opt_param);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [overlap] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("COUNT NUMBER OF OVERLAPS BETWEEN TWO BED FILES");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("bed")) vrb.error("Specify phenotypes with --bed [file.bed]");
	if (!D.options.count("ann")) vrb.error("Specify annotations with --ann [file.bed]");
	if (!D.options.count("out")) vrb.error("Specify output file with --out [file.out]");

	//---------------
	// 6. READ FILES
	//---------------
	D.readAnnotations(D.options["ann"].as < string > ());

	if (D.options.count("permute")) shuffle(D.annotation_type.begin(), D.annotation_type.end(), rng.getEngine());

	//----------------
	// 8. RUN ANALYSIS
	//----------------
	D.builtAnnotationtrees();
	D.processingPhenotypes(D.options["bed"].as < string > (), D.options["out"].as < string > ());
}

