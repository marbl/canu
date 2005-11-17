
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
 * $Id: BuildUnitigs.cc,v 1.2 2005-11-17 22:08:48 kli1000 Exp $
 * $Revision: 1.2 $
*/

static const char BUILD_UNITIGS_MAIN_CM_ID[] = "$Id: BuildUnitigs.cc,v 1.2 2005-11-17 22:08:48 kli1000 Exp $";

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
}

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

   // Open Frag store
   BestOverlapGraph::fragStoreHandle = openFragStore( FRG_Store_Path, "r");

   // Open and initialize Overlap store
   Open_OVL_Store (my_store, OVL_Store_Path);
   last = Last_Frag_In_OVL_Store (my_store);
   Init_OVL_Stream (my_stream, first, last, my_store);

   // Allocate and Initialize fragLength array
   // Seems like this should be done privately
   BestOverlapGraph::fragLength = new uint16[last+1];
   memset( BestOverlapGraph::fragLength, 0, sizeof(uint16)*(last+1));

   // Initialize our three different types of Best Overlap Graphs
   AS_BOG::ErateScore erScore(last);
   AS_BOG::LongestEdge lenScore(last);
   AS_BOG::LongestHighIdent lenIdent(last,2.0);

   // Put the three graphs into a vector, so we can step through them
   AS_BOG::BOG_Runner bogRunner;
   bogRunner.push_back(&erScore);
   bogRunner.push_back(&lenScore);
   bogRunner.push_back(&lenIdent);

   bogRunner.processOverlapStream(my_store, my_stream);


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
		utg.writeIUMtoFile("er.ium");
		break;
		case 1:
		utg.writeIUMtoFile("len.ium");
		break;
		case 2:
		utg.writeIUMtoFile("lenid.ium");
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
   closeFragStore( BestOverlapGraph::fragStoreHandle ); 
   // Shouldn't these both be n a  destructor in BOG?
   delete[] BestOverlapGraph::fragLength;
   delete_ReadStruct(BestOverlapGraph::fsread);

   return  0;
}
