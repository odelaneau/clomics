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

#include "mode_build/build_data.h"
#include "mode_score/score_data.h"
#include "mode_locate/locate_data.h"
#include "mode_topo/topo_data.h"
#include "mode_corr/corr_data.h"
#include "mode_quantify/quantify_data.h"
#include "mode_hic/hic_data.h"
#include "mode_call/call_data.h"
#include "mode_trans/trans_data.h"
#include "mode_motif_cis/motif_cis_data.h"
#include "mode_motif_trans/motif_trans_data.h"
#include "mode_span/span_data.h"
#include "mode_trans_filter/trans_filter_data.h"
#include "mode_stat/stat_data.h"
#include "mode_mbed/mbed_data.h"
#include "mode_merge/merge_data.h"
#include "mode_reptime/reptime_data.h"
#include "mode_overlap/overlap_data.h"
#include "mode_gpeak/gpeak_data.h"
#include "mode_vcfcollapse/vcfcollapse_data.h"
#include "mode_binding/binding_data.h"

void printModes(){
    vrb.ctitle("Usage:");
    vrb.print("  clomics [mode] [options]");
    vrb.print("  eg: clomics cis --help");
    vrb.ctitle("Available modes:");
    vrb.print("  corr     compute all correlations");
    vrb.print("  build    clusterize phenotype data");
    vrb.print("  score    score clusterization of phenotype data");
    vrb.print("  locate   localize tree nodes");
    vrb.print("  hic      compute all hiC interaction");
    vrb.print("  topo     characterize topology");
    vrb.print("  call     characterize topology");
    vrb.print("  motif    characterize topology");
    vrb.print("  stat     characterize topology");
    vrb.print("  mbed     modules 2 bed");
    vrb.print("  merge    merge multiple BEDs together");
    vrb.print("  reptime  extrapolate replication times for intervals from single points estimates");
    vrb.print("  overlap  count overlap by type between two BEDs");
    vrb.print("  gpeak    get module ids for peaks");
    vrb.print("  vcoll    count minor alleles falling in modules");
    vrb.print("  bind     generate sequences at binding positions");
}

int main(int argc, char ** argv) {

	//1. Start timing
	timer running_timer;

	//2. Open LOG file if necessary
	for (int a = 1 ; a < argc ; a ++) {
		if ((strcmp(argv[a], "--log") == 0) && !vrb.open_log(string(argv[a+1]))) vrb.error("Impossible to open log file!");
		if (strcmp(argv[a], "--silent") == 0) vrb.set_silent();
	}

	//3. Print header on screen
	vrb.ctitle("clomics");
	vrb.bullet("Authors : Olivier DELANEAU / Emmanouil DERMITZAKIS");
	vrb.bullet("Contact : olivier.delaneau@gmail.com");
	vrb.bullet("Webpage : XXX");
	vrb.bullet("Version : XXX");
	vrb.bullet("Date    : " + running_timer.date());

	//4. Switch mode
	vector < string > args;
    if (argc < 2){
        printModes();
        vrb.error("Not enough options, check command line!");
    }
	for (int a = 2 ; a < argc ; a ++) args.push_back(string(argv[a]));

	//5.1. BUILD mode
	if (strcmp(argv[1], "build") == 0) build_main(args);

	//5.2. SCORE mode
	else if (strcmp(argv[1], "score") == 0) score_main(args);

	//5.3. LOCATE mode
	else if (strcmp(argv[1], "locate") == 0) locate_main(args);

	//5.4. TOPO mode
	else if (strcmp(argv[1], "topo") == 0) topo_main(args);

	//5.5. TOPO mode
	else if (strcmp(argv[1], "corr") == 0) corr_main(args);

	//5.6. QUANTIFY mode
	else if (strcmp(argv[1], "quantify") == 0) quantify_main(args);

	//5.7. QTL mode
	else if (strcmp(argv[1], "hic") == 0) hic_main(args);

	//5.7. CALL mode
	else if (strcmp(argv[1], "call") == 0) call_main(args);

	//5.8. MOTIF mode
	else if (strcmp(argv[1], "motif-cis") == 0) motif_cis_main(args);

	//5.8. MOTIF mode
	else if (strcmp(argv[1], "motif-trans") == 0) motif_trans_main(args);

	//5.9. TRANS mode
	else if (strcmp(argv[1], "trans") == 0) trans_main(args);

	//5.9. TRANS mode
	else if (strcmp(argv[1], "span") == 0) span_main(args);

	//5.10. TRANS_FILTER mode
	else if (strcmp(argv[1], "trans_filter") == 0) trans_filter_main(args);

	//5.11. TRANS_FILTER mode
	else if (strcmp(argv[1], "stat") == 0) stat_main(args);

	//5.12. MBED mode
	else if (strcmp(argv[1], "mbed") == 0) mbed_main(args);

	//5.13. MERGE mode
	else if (strcmp(argv[1], "merge") == 0) merge_main(args);

	//5.14. REPTIME mode
	else if (strcmp(argv[1], "reptime") == 0) reptime_main(args);

	//5.15. REPTIME mode
	else if (strcmp(argv[1], "overlap") == 0) overlap_main(args);

	//5.16. REPTIME mode
	else if (strcmp(argv[1], "gpeak") == 0) gpeak_main(args);

	//5.17. VCOLL mode
	else if (strcmp(argv[1], "vcoll") == 0) vcfcollapse_main(args);

	//5.18. BIND mode
	else if (strcmp(argv[1], "bind") == 0) binding_main(args);

	//5.19. UNRECOGNIZED mode
    else{
        printModes();
        vrb.error("Unrecognized clomics mode!");
    }

	//5. Terminate
	vrb.title("Running time: " + stb.str(running_timer.abs_time()) + " seconds");
	vrb.close_log();
}
