
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
* Module:  FindContainedfragments.c
* Description:
*   Dump (to stdout) the contents of an overlap store
* 
*    Programmer:  E. Venter
*       Started:  11 July 2005
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: FindContainedFrags.cc,v 1.15 2005-10-06 21:39:20 eliv Exp $
 * $Revision: 1.15 $
*/

static const char CM_ID[] = "$Id: FindContainedFrags.cc,v 1.15 2005-10-06 21:39:20 eliv Exp $";

//  System include files

#include<map>
#include<set>
#include<vector>
#include<cmath>   // for ceil and log10
#include<cstdlib> // for abs(int)
#include<iostream>

#include"AS_BOG_BestOverlapGraph.hh"

using std::cout;
using std::endl;
using std::map;
using std::set;
using AS_BOG::BestOverlapGraph;

//  Local include files
extern "C" {
#include "OlapStoreOVL.h"
#include "AS_PER_fragStore.h"
}

int  main
    (int argc, char * argv [])

{
   OVL_Store_t  * my_store;
   OVL_Stream_t  * my_stream;
   uint32  first,last;			// IUID of first/last fragments in olap store
   first = 1;

   // Get path/names of olap and frg stores from command line
   const char* OVL_Store_Path = argv[1];
   const char* FRG_Store_Path = argv[2];

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

   // Go through the overlap stream, and populate the 3 overlap graphs
   bogRunner.processOverlapStream( my_stream );

   // Free/clean up the frag/overlap store/stream handles
   Free_OVL_Stream( my_stream );
   Free_OVL_Store( my_store );
   closeFragStore( BestOverlapGraph::fragStoreHandle ); 

   // Compute the width (number of digits) in the largest IUID
   int pad = static_cast<int>(ceil( log10( last )));

   // For each IUID
   for(int i = 1; i <= last; i++) {
       AS_BOG::BestEdgeOverlap* five[bogRunner.metrics.size()];
       AS_BOG::BestEdgeOverlap* three[bogRunner.metrics.size()];

       // For each graph type
       for(int j = 0; j < bogRunner.metrics.size(); j++)  { // output olap graph

           // Retrieve the best overlaps from the BOG for the current IUID
           five[j] = bogRunner.metrics[j]->getBestEdge( i, AS_BOG::FIVE_PRIME );
           three[j] = bogRunner.metrics[j]->getBestEdge( i, AS_BOG::THREE_PRIME );

           // Why isn't this just outside of the for j loop?
           if ( j == bogRunner.metrics.size()-1 ) {

               // Get the 5' overlap for the current IUID
               CDS_IID_t b0 = five[0]->frag_b_id;
               AS_BOG::BestEdgeOverlap* f1 = five[1];
               AS_BOG::BestEdgeOverlap* f2 = five[2];

               // Set up output formatting
               cout.flags(std::ios_base::left);
               cout.width(pad);

               // Output the fragments 5' ends in_degree, for each graph
               cout << i <<" 5' "<< five[0]->in_degree << " "
                    << f1->in_degree <<" "<< f2->in_degree;

               // Report whether best overlaps agree between graphs
               if (b0 == five[1]->frag_b_id && b0 == five[2]->frag_b_id) {
                   cout <<" best "; cout.width(pad * 3 + 2);
                   cout << b0 ;
               } else { 
                   cout<< " diff "          ; cout.width(pad);
                   cout<< b0 <<" "          ; cout.width(pad);
                   cout<< f1->frag_b_id<<" "; cout.width(pad);
                   cout<< f2->frag_b_id ;
               }

               // Get the 3' overlap for the current IUID
               b0 = three[0]->frag_b_id;
               f1 = three[1];
               f2 = three[2];

               // Output the fragments 3' ends in_degree, for each graph
               cout << " 3' "<< three[0]->in_degree << " "
                    << f1->in_degree <<" "<< f2->in_degree;

               // Report whether best overlaps agree between graphs
               if (b0 == three[1]->frag_b_id && b0 == three[2]->frag_b_id) {
                   cout <<" best "; cout.width(pad);
                   cout << b0 << endl;
               } else { 
                   cout<<" diff ";            cout.width(pad);
                   cout<< b0 <<" ";           cout.width(pad);
                   cout<< f1->frag_b_id<<" "; cout.width(pad);
                   cout<< f2->frag_b_id << endl;
               }

           }
       } // end for each metric
   } // end for each fragment

   // Remove transitive containments
   for(int j = 0; j < bogRunner.metrics.size(); j++)  {
       bogRunner.metrics[j]->transitiveContainment();
   }

   // We should typedef this map for best_containments
   map<CDS_IID_t,AS_BOG::BestContainment> c1 = bogRunner.metrics[0]->_best_containments;

   // Iterate through all the containees, this may miss containees that exists
   //   by other bogRunner.metrics/graph types
   for(map<CDS_IID_t,AS_BOG::BestContainment>::const_iterator it = c1.begin();
           it != c1.end(); it++)
   {

       // First is containee IUID, Second is container IUID
       CDS_IID_t id = it->first;
       AS_BOG::BestContainment bst = it->second;


	//  1: erScore
       cout << id << " c1 by "<< bst.container << " " << bst.score << " " << bst.a_hang <<
            " " << bst.b_hang <<" sameOrient " << bst.sameOrientation << endl;

	//  2: lenScore
       AS_BOG::BestContainment* b2 = bogRunner.metrics[1]->getBestContainer( id );
       if ( b2 != NULL )
           cout << id << " c2 by "<< b2->container << " " << b2->score << " " << b2->a_hang <<
            " " << b2->b_hang <<" sameOrient " << b2->sameOrientation << endl;
       
        //  3: lenIdent
       AS_BOG::BestContainment* b3 = bogRunner.metrics[2]->getBestContainer( id );
       if ( b3 != NULL )
           cout << id << " c3 by "<< b3->container << " " << b3->score << " " << b3->a_hang <<
           " " << b3->b_hang << " sameOrient " << b3->sameOrientation << endl;

   }

   // Shouldn't these both be n a  destructor in BOG?
   delete[] BestOverlapGraph::fragLength;
   delete_ReadStruct(BestOverlapGraph::fsread);

   return  0;
}
