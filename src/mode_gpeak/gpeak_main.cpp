#include "gpeak_data.h"

void gpeak_main(vector < string > & argv) {
	gpeak_data D;
	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();
	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("tree", boost::program_options::value< string >(), "Input tree file.")
		("mod", boost::program_options::value< string >(), "Input module file.")
		("out", boost::program_options::value< string >(), "Output file for peaks per module");

	D.option_descriptions.add(opt_files);
	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [gpeak] command line :" << string(e.what()) << endl;
		exit(0);
	}
	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("DETERMINE PEAKS PER MODULE");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}
	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("tree"))  vrb.error("Specify one tree file using --tree");
	if (!D.options.count("mod"))  vrb.error("Specify module IDs using --mod");
	if (!D.options.count("out")) vrb.error("Specify output bed file with --out [file.bed]");

	//---------------
	// 6. READ FILES
	//---------------
	D.readTree(D.options["tree"].as < string > ());
	D.readModules(D.options["mod"].as < string > ());

	D.process(D.options["out"].as < string > ());
}
