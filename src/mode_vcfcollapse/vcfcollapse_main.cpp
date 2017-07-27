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

#include "vcfcollapse_data.h"

void vcfcollapse_main(vector < string > & argv) {
	vcfcollapse_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();
	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("vcf", boost::program_options::value< string >(), "Genotypes in VCF/BCF format.")
		("bed", boost::program_options::value< string >(), "Annotations in BED format.")
		("out", boost::program_options::value< string >(), "Output VCF file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mFilters\33[0m");
	opt_parameters.add_options()
		("normal", "Normal transform quantifications")
		("filter-max-maf", boost::program_options::value< float >()->default_value(0.05), "Maximal MAF for input variant to be considered.")
		("filter-min-maf", boost::program_options::value< float >()->default_value(0.05), "Minimum MAF for output variant to be considered.");

	D.option_descriptions.add(opt_files).add(opt_parameters);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	boost::program_options::variables_map options;
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [vcfcollapse] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("COLLPASE RARE VARIANTS ACCORDING TO ANNOTATIONS");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("vcf")) vrb.error("Genotype data needs to be specified with --vcf [file.vcf]");
	if (!D.options.count("bed")) vrb.error("Annotation data needs to be specified with --bed [file.bed]");
	if (!D.options.count("out")) vrb.error("Output needs to be specified with --out [file.out]");

	//TO DO CHECK PARAMETER VALUES
	D.param_max_maf = D.options["filter-max-maf"].as < float > ();
	D.param_min_maf = D.options["filter-min-maf"].as < float > ();
	vrb.bullet("input MAF <= " + stb.str(D.param_min_maf));
	vrb.bullet("output MAF >= " + stb.str(D.param_max_maf));

	//------------------------------------------
	// 5. READ FILES / INITIALIZE / RUN ANALYSIS
	//------------------------------------------
	D.processBasicOptions();
	D.readSampleFromVCF(D.options["vcf"].as < string > ());
	D.mergeSampleLists();
	D.readAnnotations(D.options["bed"].as < string > ());
	D.buildTrees();
	D.readGenotypes(D.options["vcf"].as < string > ());
	D.collapseGenotypes(D.options.count("normal"));
	D.writeGenotypes(D.options["out"].as < string > ());
}
