#include "trans_filter_data.h"

void trans_filter_main(vector < string > & argv) {
	trans_filter_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("inter", boost::program_options::value< string >(), "Interaction files")
		("out", boost::program_options::value< string >(), "Output file");

	boost::program_options::options_description opt_param ("\x1B[32mParameters\33[0m");
	opt_param.add_options()
		("bin", boost::program_options::value< int >()->default_value(5000), "Merge together positions close by this size")
		("hic", boost::program_options::value< double >()->default_value(2), "Minimal Hi-C value")
		("pva", boost::program_options::value< double >()->default_value(0.001), "Minimal correlation p-value");

	D.option_descriptions.add(opt_files).add(opt_param);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [trans_filter] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("FILTER TRANS INTERACTION FILE");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}
	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("inter")) vrb.error("Specify interaction files with --inter [file.bed]");
	if (!D.options.count("out")) vrb.error("Specify output file with --out [file.bed]");

	D.bin = D.options["bin"].as < int > ();
	D.hic = D.options["hic"].as < double > ();
	D.pva = D.options["pva"].as < double > ();
	vrb.bullet("Merging bin size = " + stb.str(D.bin));
	vrb.bullet("Hi-C value threshold = " + stb.str(D.hic));
	vrb.bullet("P-value = " + stb.str(D.pva));

	//---------------
	// 6. READ FILES
	//---------------
	D.readInteractions(D.options["inter"].as < string > ());
	D.writeInteractions(D.options["out"].as < string > ());

}
