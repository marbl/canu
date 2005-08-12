
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
 * $Id: FindContainedFrags.cc,v 1.7 2005-08-12 20:53:38 eliv Exp $
 * $Revision: 1.7 $
*/

static const char CM_ID[] = "$Id: FindContainedFrags.cc,v 1.7 2005-08-12 20:53:38 eliv Exp $";

//  System include files

#include<map>
#include<set>
#include<vector>
#include<cstdlib> // for abs(int)
#include<iostream>

#include"AS_BOG_BestOverlapGraph.hh"

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

struct MultiContain {
   int iid;
   std::set<int> in;
   std::map<int,float> equal;
};

int  main
    (int argc, char * argv [])

{
   OVL_Store_t  * my_store;
   OVL_Stream_t  * my_stream;
   Long_Olap_Data_t  olap;
   uint32  first,last;
   first = 1;

   const char* OVL_Store_Path = argv[1];
   const char* FRG_Store_Path = argv[2];

   my_store = New_OVL_Store ();
   my_stream = New_OVL_Stream ();

   BestOverlapGraph::fragStoreHandle = openFragStore( FRG_Store_Path, "r");

   Open_OVL_Store (my_store, OVL_Store_Path);
   last = Last_Frag_In_OVL_Store (my_store);
   Init_OVL_Stream (my_stream, first, last, my_store);

   BestOverlapGraph::fragLength = new uint16[last+1];
   memset( BestOverlapGraph::fragLength, 0, sizeof(uint16)*(last+1));

   MultiContain *multiContain;
   multiContain = new MultiContain[last+1];
   //vector<MultiContain> multiContain(last+1);
   int i;
   for(i=0; i <= last; i++) {
       MultiContain m;
       m.iid = i;
       multiContain[i] = m;
   }

   AS_BOG::ErateScore erScore(last);
   AS_BOG::LongestEdge lenScore(last);
   AS_BOG::LongestHighIdent lenIdent(last,2.0);
   vector<BestOverlapGraph *> metrics;

   metrics.push_back(&erScore);
   metrics.push_back(&lenScore);
   metrics.push_back(&lenIdent);

   int j;
   while  (Next_From_OVL_Stream (&olap, my_stream))
     {
        float erate = Expand_Quality(olap.corr_erate) * 100;
         if ( olap.a_hang == 0 && olap.b_hang == 0 )
         {
             multiContain[ olap.a_iid ].equal[ olap.b_iid ] = erate;
         }
         else if ( olap.a_hang >= 0 && olap.b_hang <= 0 )
         {
             multiContain[ olap.b_iid ].in.insert( olap.a_iid );
         }
         else if ( olap.a_hang <= 0 && olap.b_hang >= 0 )
         {
             multiContain[ olap.a_iid ].in.insert( olap.b_iid );
         } else {
             // no containment, so score
            for( j = 0; j < metrics.size(); j++) 
                metrics[j]->score( olap );

         }
     }
   Free_OVL_Stream( my_stream );
   Free_OVL_Store( my_store );
   closeFragStore( BestOverlapGraph::fragStoreHandle ); 

   for( i = 1; i <= last; i++) {
       AS_BOG::BestEdgeOverlap* five[metrics.size()];
       AS_BOG::BestEdgeOverlap* three[metrics.size()];

       for( j = 0; j < metrics.size(); j++)  {
           five[j] = metrics[j]->getBestEdge( i, AS_BOG::FIVE_PRIME );
           three[j] = metrics[j]->getBestEdge( i, AS_BOG::THREE_PRIME );

           if ( j == metrics.size()-1 ) {
               CDS_IID_t b0 = five[0]->frag_b_id;
               AS_BOG::BestEdgeOverlap* f1 = five[1];
               AS_BOG::BestEdgeOverlap* f2 = five[2];
               cout.width(6);
               cout << i <<" 5' in "<< five[0]->in_degree << " "
                    << f1->in_degree <<" "<< f2->in_degree;

               if (b0 == five[1]->frag_b_id && b0 == five[2]->frag_b_id) {
                   cout <<" best "; cout.width(17);
                   cout << b0 ;
               } else { 
                   cout<< " diff "          ; cout.width(5);
                   cout<< b0 <<" "          ; cout.width(5);
                   cout<< f1->frag_b_id<<" "; cout.width(5);
                   cout<< f2->frag_b_id ;
               }
               b0 = three[0]->frag_b_id;
               f1 = three[1];
               f2 = three[2];
               cout << " ;3' in "<< three[0]->in_degree << " "
                    << f1->in_degree <<" "<< f2->in_degree;

               if (b0 == three[1]->frag_b_id && b0 == three[2]->frag_b_id) {
                   cout <<" best "; cout.width(5);
                   cout << b0 << endl;
               } else { 
                   cout<<" diff ";            cout.width(5);
                   cout<< b0 <<" ";           cout.width(5);
                   cout<< f1->frag_b_id<<" "; cout.width(5);
                   cout<< f2->frag_b_id << endl;
               }
           }
       }
       if ( -1 != multiContain[i].iid && multiContain[i].in.size() > 0 ) {
           set<int> iids = multiContain[i].in;
//           cout << i << " in ";
           for(set<int>::const_iterator it = iids.begin(); it != iids.end(); it++) {
//               cout << *it << " ";
               if (multiContain[*it].in.size() > 0) {} // multi contain
           }
 //          cout << endl;
       }
       if ( -1 != multiContain[i].iid && multiContain[i].equal.size() > 0 ) {
//           cout << i << " ident ";
           map<int,float> iids = multiContain[i].equal;
           for(map<int,float>::const_iterator it = iids.begin(); it != iids.end(); it++){
//               cout << it->first << "," << it->second << " ";
           }
//           cout << endl;
       }
   }
   delete[] multiContain;
   delete[] BestOverlapGraph::fragLength;
   delete_ReadStruct(BestOverlapGraph::fsread);
   return  0;
}
