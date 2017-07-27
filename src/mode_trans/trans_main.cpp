#include "trans_data.h"

void trans_main(vector < string > & argv) {
	trans_data D;

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
		("threshold", boost::program_options::value< double >()->default_value(1.0), "More output for significant hits")
		("centromere", boost::program_options::value< vector < int > >()->multitoken(), "Centromere positions to calculate relative positions")
		("bincentromere", boost::program_options::value< int >(), "Count number of tests per relative positions bin")
		("null", boost::program_options::value< double >()->default_value(0.01), "Percentage null to be reported in the output")
		("bin", boost::program_options::value< double >(), "Mean correlation per bin")
		("full", "Full report (BIG)");

	boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
	opt_parallel.add_options()
		("regions", boost::program_options::value< vector < string > >()->multitoken(), "Region of interest.");

	D.option_descriptions.add(opt_files).add(opt_param).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [trans] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("TRANS INTERACTION WITH HI-C SUPPORT");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}
	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("bed")) vrb.error("Specify phenotypes with --bed [file.bed]");
	//if (!D.options.count("hic")) vrb.error("Specify hic data with --hic [file.mat file1.norm file2.norm bin_size]");
	if (!D.options.count("regions")) vrb.error("Specify two regions with --regions [chr1 chr2]");
	if (!D.options.count("out")) vrb.error("Specify output tree with --out [file.out]");

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
	D.allocate();
	D.normalizePhenotypes1();
	D.normalizePhenotypes2();

	if (D.options.count("full")) {
		D.computeAllInteractionsWithDesc(D.options["out"].as < string > ());
	} else if (D.options.count("bin")) {
		D.computeBinnedInteractions(D.options["out"].as < string > (), D.options["bin"].as < double > ());
	} else if (D.options.count("hic")) {
		vector < string > trans_list = D.options["hic"].as < vector < string > > ();
		if (trans_list.size() != 4) vrb.error("--hic take 4 arguments: interation_matrix normalization_vector1 normalization_vector2 bin_size");
		D.hic_bin1 = atoi(trans_list[3].c_str());
		D.hic_bin2 = atoi(trans_list[3].c_str());
		vrb.bullet("Bin size of Hi-C data = " + stb.str(D.hic_bin1));
		D.thinin = D.options["null"].as < double > ();
		vrb.bullet("Percentage of null interactions to be reported = " + stb.str(D.thinin*100, 1));
		D.readHiCnorms1(trans_list[1]);
		D.readHiCnorms2(trans_list[2]);
		D.readHiCvalues(trans_list[0]);
		if (D.options.count("centromere")) {
				vector < int > centromere_list = D.options["centromere"].as < vector < int > > ();
				assert(centromere_list.size() == 2);
				D.computeRelativePositions(centromere_list[0], centromere_list[1]);
		}
		D.computeHicInteractions(D.options["out"].as < string > ());
	} else if (D.options.count("bincentromere")) {
		vector < int > centromere_list = D.options["centromere"].as < vector < int > > ();
		assert(centromere_list.size() == 2);
		D.computeRelativePositions(centromere_list[0], centromere_list[1]);
		D.binRelativePositions(D.options["bincentromere"].as < int > (), D.options["out"].as < string > ());
	} else if (D.options["threshold"].defaulted()) {
		D.computeAllInteractions(D.options["out"].as < string > ());
	} else {
		if (D.options.count("centromere")) {
			vector < int > centromere_list = D.options["centromere"].as < vector < int > > ();
			assert(centromere_list.size() == 2);
			D.computeRelativePositions(centromere_list[0], centromere_list[1]);
		}
		D.computeSigInteractions(D.options["out"].as < string > (), D.options["threshold"].as < double > ());
	}
}
