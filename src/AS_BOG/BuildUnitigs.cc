
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
 * $Id: BuildUnitigs.cc,v 1.29 2007-10-23 13:40:01 eliv Exp $
 * $Revision: 1.29 $
*/

static const char BUILD_UNITIGS_MAIN_CM_ID[] = "$Id: BuildUnitigs.cc,v 1.29 2007-10-23 13:40:01 eliv Exp $";

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

void outputHistograms(AS_BOG::UnitigGraph *);

int  main (int argc, char * argv [])

{
   OverlapStore  * my_store;
   OVSoverlap      olap;

   fprintf(stderr, "%s\n\n", BUILD_UNITIGS_MAIN_CM_ID);

   // Get path/names of olap and frg stores from command line
   const char* OVL_Store_Path;
   const char* GKP_Store_Path;

   vector<float> erates;
//   erates.push_back(2.5);
   erates.push_back(1.5);
//   erates.push_back(1.0);
   long genome_size=0;
   int ch;
   bool argsDone=false;
   optarg = NULL;
   while(!argsDone && (ch = getopt(argc, argv,"O:G:e:s:bk"))) {
       switch(ch) {
           case -1: argsDone=true;break;
           case 'G':
               GKP_Store_Path = strdup(optarg);
               assert( GKP_Store_Path != NULL ); break;
           case 'O':
               OVL_Store_Path = strdup(optarg);
               assert( OVL_Store_Path != NULL ); break;
           case 'b':
               BogOptions::unitigIntersectBreaking = true; break;
           case 'e': {
               erates.clear();
               std::string eOpt(optarg);
               std::istringstream iStr(eOpt);
               float t;
               while (iStr >> t) erates.push_back(t);
               } break;
           case 'k':
               BogOptions::ejectUnhappyContained = true; break;
           case 's':
               genome_size = atol(optarg); break;
           default:
        fprintf(stderr,"Unrecognized option -%c optarg %s\n\n",optopt, optarg);
        fprintf(stderr, "%s: -O <OVL Store Path> -G <GKP Store Path>\n", argv[0]);
        fprintf(stderr, "[-s <genome size>]\n");
        fprintf(stderr, "  If the genome size is set to 0, this will cause the unitigger\n");
        fprintf(stderr, "  to try to estimate the genome size based on the constructed\n");
        fprintf(stderr, "  unitig lengths.\n");
        fprintf(stderr, "[-b] Break promisciuous unitigs at unitig intersection points\n");
        fprintf(stderr, "[-e] Space separated list of erates to generate unitig for\n");
        fprintf(stderr, "     default is -e '2.5 1.5 1.0'\n");
        fprintf(stderr, "[-k] Kick out unhappy contained mated reads into singleton unitigs\n");
        fprintf(stderr, " \n");
           exit(1);
       }
   }

   std::cerr << "Genome Size: " << genome_size << std::endl;

   my_store = AS_OVS_openOverlapStore(OVL_Store_Path);

   AS_BOG::MateChecker mateChecker;
   int numFrgsInGKP = mateChecker.readStore(GKP_Store_Path);
   // must be before creating the scoring objects, because it sets their size
//   AS_BOG::BOG_Runner bogRunner(getLastElemFragStore() need to fix to get size of gkpStore;
   //AS_BOG::BOG_Runner bogRunner(AS_OVS_lastFragInStore(my_store));
   AS_BOG::BOG_Runner bogRunner(numFrgsInGKP);

   // Initialize our three different types of Best Overlap Graphs
   //AS_BOG::ErateScore erScore;
   //AS_BOG::LongestEdge lenScore;

   int i;
   for(i=0; i < erates.size(); i++) {
       bogRunner.push_back( new AS_BOG::LongestHighIdent( erates[i] ) );
   }
   bogRunner.processOverlapStream(my_store, GKP_Store_Path);

   ////////////////////////////////////////////////////////////////////////////

   for(i=0; i<bogRunner.size(); i++){
//	AS_BOG::ChunkGraph *cg = new AS_BOG::ChunkGraph();
	AS_BOG::PromiscuousChunkGraph *cg = new AS_BOG::PromiscuousChunkGraph();
	//cg.checkInDegree(bogRunner.metrics[i]);
	cg->build(bogRunner.metrics[i]);
	std::cout << "Num Fragments: " << cg->getNumFragments() << std::endl;
	std::cout << "Num Singletons: " << cg->countSingletons() << std::endl;
	std::cout << "Num Containees: " << bogRunner.metrics[i]->_best_containments.size() << std::endl;
	std::cout << std::endl;

	AS_BOG::UnitigGraph utg(bogRunner.metrics[i]);
	std::cerr << "Building Unitigs.\n" << std::endl;
	utg.build(cg);

	std::cerr << "Reporting.\n" << std::endl;
    mateChecker.checkUnitigGraph(utg);

    char fileStr[16];
    int mismatch = bogRunner.metrics[i]->getThreshold();
    sprintf( fileStr, "len%d.ium",mismatch);
    fprintf( stderr, "Made unitigs at %d mismatch rate\n",mismatch);
    std::cerr << "Setting Global Arrival Rate.\n";
    // should be number of Random frags when that's supported
    float globalARate = utg.getGlobalArrivalRate(cg->getNumFragments(), genome_size);
    AS_BOG::Unitig::setGlobalArrivalRate(globalARate);
    std::cerr << "Global Arrival Rate: " << globalARate << std::endl;
    std::cerr << std::endl<< "There were " << utg.unitigs->size() << " unitigs generated.\n";

    utg.writeIUMtoFile(fileStr);
    outputHistograms( &utg );
	std::cerr << "///////////////////////////////////////////////////////////\n" << std::endl;

    delete bogRunner.metrics[i];
    delete cg;
   }


   std::cerr<<"Done with unitiger specifics.\n"<<std::endl;

   ////////////////////////////////////////////////////////////////////////////

   std::cerr << "Cleaning up." << std::endl;
   // Free/clean up the frag/overlap store/stream handles
   AS_OVS_closeOverlapStore(my_store);
   // Shouldn't these both be n a  destructor in BOG?
   delete[] BestOverlapGraph::fragLength;

   return  0;
}


void outputHistograms(AS_BOG::UnitigGraph *utg) {
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
        if (arate < 0) {
            std::cerr << "Negative Local ArrivalRate " << arateF << " id " << u->id()
                << " arate " << arate << std::endl;
        }
        zork.sum_arrival = zork.min_arrival = zork.max_arrival = arate;

        zork.sum_rs_frags=zork.min_rs_frags=zork.max_rs_frags=0;
        zork.sum_nr_frags=zork.min_nr_frags=zork.max_nr_frags=0;

        add_to_histogram(length_of_unitigs_histogram, u->getLength(), &zork);
        add_to_histogram(covg_histogram, covg, &zork);
        add_to_histogram(arate_histogram, arate, &zork);
    }
    std::cerr << "Length of Unitigs histogram" << std::endl;
    print_histogram(stderr,length_of_unitigs_histogram, 0, 1);
    std::cerr << "Unitig Coverage Stat histogram" << std::endl;
    print_histogram(stderr,covg_histogram, 0, 1);
    std::cerr << "Unitig Arrival Rate histogram" << std::endl;
    print_histogram(stderr,arate_histogram, 0, 1);
    free_histogram( length_of_unitigs_histogram );
    free_histogram( covg_histogram );
    free_histogram( arate_histogram );
}
