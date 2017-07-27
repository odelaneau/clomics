#include "binding_data.h"

void binding_main(vector < string > & argv) {
	binding_data D;
	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();
	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("ref", boost::program_options::value< string >(), "Input Reference")
		("pos", boost::program_options::value< string >(), "Input Position")
		("out", boost::program_options::value< string >(), "Output sequences");

	D.option_descriptions.add(opt_files);
	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [binding] command line :" << string(e.what()) << endl;
		exit(0);
	}
	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("GENERATE SEQUENCES WITH ALTERNATIVE ALLELES");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}
	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("ref"))  vrb.error("Specify the reference genome with --ref");
	if (!D.options.count("pos"))  vrb.error("Specify the variable position with --pos");
	if (!D.options.count("out")) vrb.error("Specify output file with --out");

	//---------------
	// 6. READ FILES
	//---------------
	D.readReference(D.options["ref"].as < string > ());
	D.readPosition(D.options["pos"].as < string > ());
	D.process();
	D.writeFasta(D.options["out"].as < string > ());
}
