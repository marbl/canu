
/**************************************************************************
 * Copyright (C) 2005, J Craig Venter Institute. All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received (LICENSE.txt) a copy of the GNU General Public 
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
/*************************************************
 * Module:  BuildUnitigs.c
 * Description:
 * 
 *    Programmer:  K. Li
 *       Started:  09/19/2005
 * 
 * Assumptions:
 * 
 * Notes:
 *
 *************************************************/

/* RCS info
 * $Id: BuildUnitigs.cc,v 1.40 2008-04-16 10:26:27 brianwalenz Exp $
 * $Revision: 1.40 $
 */

static const char BUILD_UNITIGS_MAIN_CM_ID[] = "$Id: BuildUnitigs.cc,v 1.40 2008-04-16 10:26:27 brianwalenz Exp $";

//  System include files

#include<map>
#include<set>
#include<vector>
#include<cmath>   // for ceil and log10
#include<cstdlib> // for abs(int)
#include<iostream>
#include<sstream>

#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_MateChecker.hh"

using std::cout;
using std::endl;
using std::map;
using std::set;
using std::vector;

using AS_BOG::BestOverlapGraph;
using AS_BOG::BogOptions;

//  Local include files
extern "C" {
#include "getopt.h"
#include "AS_OVS_overlapStore.h"
#include "AS_CGB_histo.h"
}

void outputHistograms(AS_BOG::UnitigGraph *, FILE *);

int
main (int argc, char * argv []) {

    OverlapStore  * ovlStore;
    OVSoverlap      olap;

    fprintf(stderr, "%s\n\n", BUILD_UNITIGS_MAIN_CM_ID);

    argc = AS_configure(argc, argv);

    // Get path/names of olap and frg stores from command line
    const char* OVL_Store_Path;
    const char* GKP_Store_Path;

    vector<float> erates;
    bool          eratesDefault = true;

    erates.push_back(1.5);

    long genome_size=0;
    int ch;
    bool argsDone=false;

    int   fragment_count_target = 0;
    char *output_prefix         = "bog";

    optarg = NULL;
    while(!argsDone && (ch = getopt(argc, argv,"B:o:O:G:e:m:s:bk"))) {
        switch(ch) {
            case -1:
                argsDone=true;
                break;

            case 'B':
                fragment_count_target = atoi(optarg);
                break;
            case 'o':
                output_prefix = optarg;
                break;

            case 'G':
                GKP_Store_Path = strdup(optarg);
                assert( GKP_Store_Path != NULL ); break;
            case 'O':
                OVL_Store_Path = strdup(optarg);
                assert( OVL_Store_Path != NULL ); break;
            case 'b':
                BogOptions::unitigIntersectBreaking = true; break;
            case 'e':
                if (eratesDefault) {
                  eratesDefault = false;
                  erates.clear();
                }
                if (NULL != strstr(optarg,"."))
                  //  Is floating point -- parts per hundred (aka percent error)
                  erates.push_back(atof(optarg) / 100.0);
                else
                  //  Is integer -- parts per thousand
                  erates.push_back(atof(optarg) / 1000.0);
                break;
            case 'k':
                BogOptions::ejectUnhappyContained = true; break;
            case 'm':
                BogOptions::badMateBreakThreshold = -atoi(optarg); break;
            case 's':
                genome_size = atol(optarg); break;
            default:
                fprintf(stderr,"Unrecognized option -%c optarg %s\n\n",optopt, optarg);
                fprintf(stderr, "usage: %s -O <OVL Store Path> -G <GKP Store Path>\n", argv[0]);
                fprintf(stderr, "\n");
                fprintf(stderr, "[-B b]       Target number of fragments per IUM batch.\n");
                fprintf(stderr, "\n");
                fprintf(stderr, "[-o prefix]  Output prefix name\n");
                fprintf(stderr, "\n");
                fprintf(stderr, "\n");
                fprintf(stderr, "[-s <genome size>]\n");
                fprintf(stderr, "  If the genome size is set to 0, this will cause the unitigger\n");
                fprintf(stderr, "  to try to estimate the genome size based on the constructed\n");
                fprintf(stderr, "  unitig lengths.\n");
                fprintf(stderr, "[-b] Break promisciuous unitigs at unitig intersection points\n");
                fprintf(stderr, "[-e] Erate to generate unitigs for; default is 1.5\n");
                fprintf(stderr, "     Multiple -e switches may be supplied; all erates will be computed.\n");
                fprintf(stderr, "[-k] Kick out unhappy contained mated reads into singleton unitigs\n");
                fprintf(stderr, "[-m] Number of bad mates in a region required to break a unitig\n");
                fprintf(stderr, "     default is 7\n");
                fprintf(stderr, " \n");
                exit(1);
        }
    }

    std::cerr << "Genome Size: " << genome_size << std::endl;

    ovlStore = AS_OVS_openOverlapStore(OVL_Store_Path);

    AS_BOG::MateChecker mateChecker;
    int numFrgsInGKP = mateChecker.readStore(GKP_Store_Path);    
    GateKeeperStore *gkpStore = openGateKeeperStore(GKP_Store_Path, FALSE);
    int numRandFrgInGKP = getNumGateKeeperRandomFragments(gkpStore);
    closeGateKeeperStore(gkpStore);
    
    
    // must be before creating the scoring objects, because it sets their size
    AS_BOG::BOG_Runner bogRunner(numFrgsInGKP);

    // Initialize our three different types of Best Overlap Graphs
    for(int i=0; i < erates.size(); i++)
        bogRunner.push_back( new AS_BOG::BestOverlapGraph( erates[i] ) );

    bogRunner.processOverlapStream(ovlStore, GKP_Store_Path);


    for(int i=0; i<bogRunner.size(); i++){
        //AS_BOG::ChunkGraph *cg = new AS_BOG::ChunkGraph();
	AS_BOG::PromiscuousChunkGraph *cg = new AS_BOG::PromiscuousChunkGraph();

	//cg.checkInDegree(bogRunner.metrics[i]);
	cg->build(bogRunner.metrics[i]);

	std::cout << "Num Fragments: " << cg->getNumFragments() << std::endl;
	std::cout << "Num Singletons: " << cg->countSingletons() << std::endl;
	std::cout << "Num Containees: " << bogRunner.metrics[i]->_best_containments.size() << std::endl;
	std::cout << std::endl;

	AS_BOG::UnitigGraph utg(bogRunner.metrics[i]);

	utg.build(cg);

        mateChecker.checkUnitigGraph(utg);        

        float globalARate = utg.getGlobalArrivalRate(numRandFrgInGKP, genome_size);
        AS_BOG::Unitig::setGlobalArrivalRate(globalARate);


        //  Ugh.  If we're doing multiple error rates, make a new
        //  prefix for each, which references the percent error
        //  allowed.
        //
        if (bogRunner.size() > 1) {
          char  prefix[FILENAME_MAX] = {0};
          sprintf(prefix, "%s_%05.2f", output_prefix,
                  AS_OVS_decodeQuality(bogRunner.metrics[i]->mismatchCutoff) * 100.0);
          utg.writeIUMtoFile(prefix, fragment_count_target);
        } else {
          utg.writeIUMtoFile(output_prefix, fragment_count_target);
        }


        {
          char  filename[FILENAME_MAX] = {0};
          if (bogRunner.size() > 1)
            sprintf(filename, "%s_%05.2f.cga.0", output_prefix,
                    AS_OVS_decodeQuality(bogRunner.metrics[i]->mismatchCutoff) * 100.0);
          else
            sprintf(filename, "%s.cga.0", output_prefix);
          FILE *stats = fopen(filename,"w");
          assert(NULL != stats);

          fprintf(stats, "Global Arrival Rate: %f\n", globalARate);
          fprintf(stats, "There were %d unitigs generated.\n", utg.unitigs->size());

          AS_BOG::BestEdgeCounts cnts = utg.countInternalBestEdges();

          fprintf(stats, "Overall best edge counts: dovetail %d oneWayBest %d neither %d\n",
                  cnts.dovetail,
                  cnts.oneWayBest,
                  cnts.neither);

          outputHistograms( &utg, stats );

          fclose(stats);
        }

        delete bogRunner.metrics[i];
        delete cg;
    }

    AS_OVS_closeOverlapStore(ovlStore);

    // Shouldn't these both be n a  destructor in BOG?
    delete[] BestOverlapGraph::fragLength;

    fprintf(stderr, "Bye.\n");

    return  0;
}


void outputHistograms(AS_BOG::UnitigGraph *utg, FILE *stats) {
    const int nsample=500;
    const int nbucket=500;
    MyHistoDataType zork;

    Histogram_t *length_of_unitigs_histogram = create_histogram(nsample,nbucket,0,TRUE);
    Histogram_t *covg_histogram = create_histogram(nsample,nbucket,0,TRUE);
    Histogram_t *arate_histogram = create_histogram(nsample,nbucket,0,TRUE);

    extend_histogram(length_of_unitigs_histogram, sizeof(MyHistoDataType),
                     myindexdata,mysetdata,myaggregate,myprintdata);
    extend_histogram(covg_histogram, sizeof(MyHistoDataType),
                     myindexdata,mysetdata,myaggregate,myprintdata);
    extend_histogram(arate_histogram, sizeof(MyHistoDataType),
                     myindexdata,mysetdata,myaggregate,myprintdata);

    AS_BOG::UnitigVector::const_iterator uiter = utg->unitigs->begin();
    for(;uiter != utg->unitigs->end(); uiter++) {

        AS_BOG::Unitig *u = *uiter;
        if (u == NULL)
            continue;
        zork.nsamples = 1;
        int numFrags = u->getNumFrags();
        zork.sum_frags = zork.min_frags = zork.max_frags = numFrags;

        int bases = u->getLength();
        zork.sum_bp = zork.min_bp = zork.max_bp = bases;

        int rho = static_cast<int>(u->getAvgRho());
        zork.sum_rho = zork.min_rho = zork.max_rho = rho;

        //long covg = static_cast<long>(rint(u->getCovStat() * 1000));
        int covg = static_cast<int>(rintf(u->getCovStat()));
        zork.sum_discr = zork.min_discr = zork.max_discr = covg;

        float arateF = u->getLocalArrivalRate() * 10000;
        int arate = static_cast<int>(rintf(arateF));
        if (arate < 0)
            fprintf(stats, "Negative Local ArrivalRate %f id %d arate %d\n", arateF, u->id(), arate);
        zork.sum_arrival = zork.min_arrival = zork.max_arrival = arate;

        zork.sum_rs_frags=zork.min_rs_frags=zork.max_rs_frags=0;
        zork.sum_nr_frags=zork.min_nr_frags=zork.max_nr_frags=0;

        add_to_histogram(length_of_unitigs_histogram, u->getLength(), &zork);
        add_to_histogram(covg_histogram, covg, &zork);
        add_to_histogram(arate_histogram, arate, &zork);
    }

    fprintf(stats, "Length of Unitigs histogram\n");
    print_histogram(stats,length_of_unitigs_histogram, 0, 1);

    fprintf(stats, "Unitig Coverage Stat histogram\n");
    print_histogram(stats,covg_histogram, 0, 1);

    fprintf(stats, "Unitig Arrival Rate histogram\n");
    print_histogram(stats,arate_histogram, 0, 1);

    free_histogram( length_of_unitigs_histogram );
    free_histogram( covg_histogram );
    free_histogram( arate_histogram );
}
