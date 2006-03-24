
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
 * $Id: BuildUnitigs.cc,v 1.6 2006-03-24 21:38:26 eliv Exp $
 * $Revision: 1.6 $
*/

static const char BUILD_UNITIGS_MAIN_CM_ID[] = "$Id: BuildUnitigs.cc,v 1.6 2006-03-24 21:38:26 eliv Exp $";

//  System include files

#include<map>
#include<set>
#include<vector>
#include<cmath>   // for ceil and log10
#include<cstdlib> // for abs(int)
#include<iostream>

#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_UnitigGraph.hh"

using std::cout;
using std::endl;
using std::map;
using std::set;
using std::vector;
using AS_BOG::BestOverlapGraph;

//  Local include files
extern "C" {
#include "OlapStoreOVL.h"
#include "AS_PER_fragStore.h"
#include "AS_CGB_myhisto.h"
}

void outputHistograms(AS_BOG::UnitigGraph *);

int  main (int argc, char * argv [])

{
   OVL_Store_t  * my_store;
   OVL_Stream_t  * my_stream;
   Long_Olap_Data_t  olap;
   uint32  first,last;			// IUID of first/last fragments in olap store
   first = 1;

   fprintf(stderr, "%s\n\n", BUILD_UNITIGS_MAIN_CM_ID);

   if(argc !=4){
      fprintf(stderr, "%s: <OVL Store Path> <FRG Store Path> <genome size>\n\n", argv[0]);
      fprintf(stderr, "  If the genome size is set to 0, this will cause the unitigger\n");
      fprintf(stderr, "  to try to estimate the genome size based on the constructed\n");
      fprintf(stderr, "  unitig lengths.\n");
      fprintf(stderr, " \n");
      return(-1);
   }

   // Get path/names of olap and frg stores from command line
   const char* OVL_Store_Path = argv[1];
   const char* FRG_Store_Path = argv[2];

   long genome_size;
   genome_size=atol(argv[3]);
   std::cerr << "Genome Size: " << genome_size << std::endl;

   my_store = New_OVL_Store ();
   my_stream = New_OVL_Stream ();

   // Open and initialize Overlap store
   Open_OVL_Store (my_store, OVL_Store_Path);
   last = Last_Frag_In_OVL_Store(my_store);
   Init_OVL_Stream (my_stream, first, last, my_store);

   // must be before creating the scoring objects, because it sets their size
   AS_BOG::BOG_Runner bogRunner(last);

   // Initialize our three different types of Best Overlap Graphs
   //AS_BOG::ErateScore erScore;
   //AS_BOG::LongestEdge lenScore;
   bogRunner.push_back( new AS_BOG::LongestHighIdent(2.5));
   bogRunner.push_back( new AS_BOG::LongestHighIdent(1.5));
   bogRunner.push_back( new AS_BOG::LongestHighIdent(1.0));

   bogRunner.processOverlapStream(my_store, my_stream, FRG_Store_Path);

   ////////////////////////////////////////////////////////////////////////////

   int i;
   for(i=0; i<3; i++){
	AS_BOG::ChunkGraph cg;
	//cg.checkInDegree(bogRunner.metrics[i]);
	cg.build(bogRunner.metrics[i]);
	std::cout << "Num Fragments: " << cg.getNumFragments() << std::endl;
	std::cout << "Num Singletons: " << cg.countSingletons() << std::endl;
	std::cout << "Num Containees: " << bogRunner.metrics[i]->_best_containments.size() << std::endl;
	std::cout << std::endl;

	//std::cerr << er_cg << std::endl;
	AS_BOG::UnitigGraph utg;
	std::cerr << "Building Unitigs.\n" << std::endl;
	utg.build(&cg, bogRunner.metrics[i], cg.getNumFragments(), genome_size=0);

	std::cerr << "Reporting.\n" << std::endl;
	//std::cout<< utg << endl;
	switch(i){
		case 0:
		utg.writeIUMtoFile("len25.ium");
		break;
		case 1:
		utg.writeIUMtoFile("len15.ium");
        outputHistograms( &utg );
		break;
		case 2:
		utg.writeIUMtoFile("len10.ium");
		break;
	}
	std::cerr << "///////////////////////////////////////////////////////////\n" << std::endl;

   }


   std::cerr<<"Done with unitiger specifics.\n"<<std::endl;

   ////////////////////////////////////////////////////////////////////////////

   std::cerr << "Cleaning up." << std::endl;
   // Free/clean up the frag/overlap store/stream handles
   Free_OVL_Stream( my_stream );
   Free_OVL_Store( my_store );
   // Shouldn't these both be n a  destructor in BOG?
   delete[] BestOverlapGraph::fragLength;

   return  0;
}


void outputHistograms(AS_BOG::UnitigGraph *utg) {
    const int nsample=500;
    const int nbucket=500;
    MyHistoDataType zork;

    Histogram_t *length_of_unitigs_histogram = create_histogram(nsample,nbucket,0,TRUE);
    extend_histogram(length_of_unitigs_histogram, sizeof(MyHistoDataType),
            myindexdata,mysetdata,myaggregate,myprintdata);

    AS_BOG::UnitigVector::const_iterator uiter = utg->unitigs.begin();
    for(;uiter != utg->unitigs.end(); uiter++) {

        AS_BOG::Unitig *u = *uiter;
        zork.nsamples = 1;
        long numFrags = u->getNumFrags();
        zork.sum_frags = numFrags;
        zork.min_frags = numFrags;
        zork.max_frags = numFrags;
        long bases = u->getSumFragLength();
        zork.sum_bp = bases;
        zork.min_bp = bases;
        zork.max_bp = bases;
        long rho = static_cast<long>(u->getAvgRho());
        zork.sum_rho = rho;
        zork.min_rho = rho;
        zork.max_rho = rho;
        //long covg = static_cast<long>(rint(u->getCovStat() * 1000));
        long covg = static_cast<long>(u->getCovStat());
        zork.sum_discr = covg;
        zork.min_discr = covg;
        zork.max_discr = covg;
        zork.sum_rs_frags=zork.min_rs_frags=zork.max_rs_frags=0;
        zork.sum_nr_frags=zork.min_nr_frags=zork.max_nr_frags=0;
        zork.sum_arrival=zork.min_arrival=zork.max_arrival=0;
        add_to_histogram(length_of_unitigs_histogram,
                u->getLength(),&zork);
    }
    print_histogram(stderr,length_of_unitigs_histogram, 0, 1);
}
