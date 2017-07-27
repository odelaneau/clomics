#include "quantify_data.h"

void quantify_main(vector < string > & argv) {
	quantify_data D;
	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();
	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("bed", boost::program_options::value< vector < string > >()->multitoken(), "List of phenotype files in BED format.")
		("tree", boost::program_options::value< vector < string > >()->multitoken(), "Input tree file.")
		("grp", boost::program_options::value< string >(), "Input pair file.")
		("out", boost::program_options::value< string >(), "Output file for quantified nodes.");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("normal", "Normal transform quantifications")
		("mean", "Quantify module by mean quantification")
		("loo", boost::program_options::value< int >(), "Quantify stability and not activity using leave-one-out strategy (0: mean / 1:variance)")
		("pca", boost::program_options::value< int >(), "Quantify activity and not stability using Principal Component Analysis (arg = PC to be used)");

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
		cerr << "Error parsing [quantify] command line :" << string(e.what()) << endl;
		exit(0);
	}
	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("QUANTIFY MODULES");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}
	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if ((D.options.count("loo") + D.options.count("pca") + D.options.count("mean")) != 1)  vrb.error("Specify one quantification method using either --pca, --loo or --mean");
	if (!D.options.count("bed")) vrb.error("Specify phenotypes with --bed [file.bed]");
	if ((D.options.count("tree") + D.options.count("grp")) != 1)  vrb.error("Specify one reference file using either --tree or --grp");
	if (!D.options.count("out")) vrb.error("Specify output bed file with --out [file.bed]");

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

	if (D.options.count("tree")) {
		vector < string > file_list = D.options["tree"].as < vector < string > > ();
		if (file_list.size() != 2) vrb.error("Option --tree requires two arguments");
		D.readTree(file_list[0]);
		D.readModules(file_list[1]);
		if (D.options.count("mean")) D.quantifyTrees(METH_MEA, 1);
		if (D.options.count("loo")) D.quantifyTrees(METH_LOO, D.options["loo"].as < int > ());
		if (D.options.count("pca")) D.quantifyTrees(METH_PCA, D.options["pca"].as < int > ());
		D.writePhenotypesTrees(D.options["out"].as < string > (), D.options.count("normal"));
	} else if (D.options.count("grp")) {
		D.readGroups(D.options["grp"].as < string > ());
		if (D.options.count("mean")) D.quantifyGroups(METH_MEA, 1);
		if (D.options.count("loo")) D.quantifyGroups(METH_LOO, D.options["loo"].as < int > ());
		if (D.options.count("pca")) D.quantifyGroups(METH_PCA, D.options["pca"].as < int > ());
		D.writePhenotypesGroups(D.options["out"].as < string > (), D.options.count("normal"));
	}
}
