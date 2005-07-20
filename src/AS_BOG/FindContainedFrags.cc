
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
 * $Id: FindContainedFrags.cc,v 1.1 2005-07-20 15:04:03 eliv Exp $
 * $Revision: 1.1 $
*/

static char CM_ID[] = "$Id: FindContainedFrags.cc,v 1.1 2005-07-20 15:04:03 eliv Exp $";

//  System include files

#include<map>
#include<set>
#include<vector>
#include<cstdlib> // for abs(int)
#include<iostream>

using std::cout;
using std::endl;
using std::map;
using std::set;
using std::vector;
using std::ostream; 

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

static FragStoreHandle fragStoreHandle;
static uint16 *fragLength;
uint16 fragLen(int iid)
{
    static ReadStructp fsread = new_ReadStruct();
    if (fragLength[ iid ] == 0) {
        uint32 clrBgn, clrEnd;
        getFragStore( fragStoreHandle, iid, FRAG_S_SEQUENCE, fsread);
        getClearRegion_ReadStruct( fsread, &clrBgn, &clrEnd, READSTRUCT_LATEST);
        fragLength[ iid ] = clrEnd - clrBgn;
    }
    return fragLength[ iid ];
}

struct ScoreOverlap {
    int curFrag;
    float  bestScore;
    int *bestOverlap;
    ScoreOverlap(int num =0) : curFrag(0) {
        bestOverlap = new int[num+1];
    }
    ~ScoreOverlap() { delete bestOverlap; }
    inline short olapLength(const Long_Olap_Data_t& olap) {
        uint16 alen = fragLen(olap.a_iid);
        if (olap.a_hang < 0) 
            return alen - abs(olap.b_hang); 
        else
            return alen - olap.a_hang; 
    }
    inline virtual bool checkForNext(const Long_Olap_Data_t& olap, double scoreReset) {
        if (curFrag != olap.a_iid) {
            curFrag  = olap.a_iid;
            bestScore = scoreReset;
            bestOverlap[ curFrag ] = 0;
            return true;
        }
        return false;
    }
    virtual float score(const Long_Olap_Data_t& olap) =0;
};

struct ErateScore : public ScoreOverlap {
    short bestLength;
    ErateScore(int num) : ScoreOverlap(num), bestLength(0) {}

    inline bool checkForNext(const Long_Olap_Data_t& olap, double scoreReset)  {
        if (ScoreOverlap::checkForNext(olap, scoreReset))
            bestLength = 0;
    }
    float score(const Long_Olap_Data_t& olap) {

        float erate = Expand_Quality(olap.corr_erate) * 100;
        checkForNext(olap,100);
        short olapLen = olapLength(olap);
        //cout << olap.a_iid<<" "<<olap.b_iid<<" "<<erate<<" "<<alen<<" "<<fragLen(olap.b_iid);
        //const char *f = olap.flipped ? "I" : "N";
        //cout <<" "<< olap.a_hang<<" "<< olap.b_hang<<" "<< f<<" "<<olapLen << endl;
        if (erate < bestScore || erate == bestScore && olapLen > bestLength ) {
            bestOverlap[ curFrag ] = olap.b_iid;
            bestScore = erate;
            bestLength = olapLen;
        }
    }
};

struct LongestEdge : public ScoreOverlap {
    LongestEdge(int num) : ScoreOverlap(num) {}
    float score(const Long_Olap_Data_t& olap) {

        uint16 alen = fragLen(olap.a_iid);
        short olapLen = olapLength(olap);
        checkForNext(olap,0);
        if (bestScore < olapLen) {
            bestOverlap[ olap.a_iid ] = olap.b_iid;
            bestScore = olapLen;
        }
    }
};
struct LongestHighIdent : public ScoreOverlap {
    float mismatchCutoff;
    LongestHighIdent(int num, float maxMismatch) : ScoreOverlap(num),
                                                   mismatchCutoff(maxMismatch) {}
    float score(const Long_Olap_Data_t& olap) {

        uint16 alen = fragLen(olap.a_iid);
        short olapLen = olapLength(olap);
        float erate = Expand_Quality(olap.corr_erate) * 100;
        checkForNext(olap,0);

        if (erate > mismatchCutoff)
            return 0;
        if (bestScore < olapLen) {
            bestOverlap[ olap.a_iid ] = olap.b_iid;
            bestScore = olapLen;
        }
    }
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

   fragStoreHandle = openFragStore( FRG_Store_Path, "r");

   Open_OVL_Store (my_store, OVL_Store_Path);
   last = Last_Frag_In_OVL_Store (my_store);
   Init_OVL_Stream (my_stream, first, last, my_store);

   fragLength = new uint16[last+1];
   memset( fragLength, 0, sizeof(uint16)*(last+1));

   MultiContain *multiContain;
   multiContain = new MultiContain[last+1];
   //vector<MultiContain> multiContain(last+1);
   int i;
   for(i=0; i <= last; i++) {
       MultiContain m;
       m.iid = i;
       multiContain[i] = m;
   }

   ErateScore erScore(last);
   LongestEdge lenScore(last);
   LongestHighIdent lenIdent(last,2.0);
   vector<ScoreOverlap *> metrics;

   metrics.push_back(&erScore);
   metrics.push_back(&lenScore);
   metrics.push_back(&lenIdent);

   int j;
   while  (Next_From_OVL_Stream (& olap, my_stream))
     {
        float erate = Expand_Quality(olap.corr_erate) * 100;
         if ( olap.a_hang >= 0 && olap.b_hang < 0 )
         {
             multiContain[ olap.b_iid ].in.insert( olap.a_iid );
         }
         else if ( olap.a_hang < 0 && olap.b_hang >= 0 )
         {
             multiContain[ olap.a_iid ].in.insert( olap.b_iid );
         }
         else if ( olap.a_hang == 0 && olap.b_hang == 0 )
         {
             multiContain[ olap.a_iid ].equal[ olap.b_iid ] = erate;
         } else {
             // no containment, so score
            for( j = 0; j < metrics.size(); j++) 
                metrics[j]->score( olap );

         }
     }
   Free_OVL_Stream (my_stream);
   Free_OVL_Store (my_store);

   for( i = 1; i <= last; i++) {
       int pick[metrics.size()];
       for( j = 0; j < metrics.size(); j++)  {
           //cout << i << " "<< typeid(*metrics[j]).name()+2 << " " << metrics[j]->bestOverlap[i] << endl;
           pick[j] = metrics[j]->bestOverlap[i];
           if ( j == metrics.size()-1 ) 
               if (pick[0] == pick[1] && pick[0] == pick[2]) {
                   //cout << i << pick[0];
               } else { 
                   cout << i << " disagree " << pick[0]<<" "<<pick[1]<<" "<<pick[2]<<endl;
               }
       }
       if ( -1 != multiContain[i].iid && multiContain[i].in.size() > 0 ) {
           set<int> iids = multiContain[i].in;
           cout << i << " in ";
           for(set<int>::const_iterator it = iids.begin(); it != iids.end(); it++) {
               cout << *it << " ";
               if (multiContain[*it].in.size() > 0) {} // multi contain
           }
           cout << endl;
       }
       if ( -1 != multiContain[i].iid && multiContain[i].equal.size() > 0 ) {
           cout << i << " ident ";
           map<int,float> iids = multiContain[i].equal;
           for(map<int,float>::const_iterator it = iids.begin(); it != iids.end(); it++){
               cout << it->first << "," << it->second << " ";
           }
           cout << endl;
       }
   }
   return  0;
}
