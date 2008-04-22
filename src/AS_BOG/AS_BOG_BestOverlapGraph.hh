
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

#ifndef INCLUDE_AS_BOG_BESTOVERLAPGRAPH
#define INCLUDE_AS_BOG_BESTOVERLAPGRAPH

#include <set>

#include "AS_BOG_Datatypes.hh"

extern "C" {
#include "AS_OVS_overlapStore.h"
}

namespace AS_BOG{

    struct BestEdgeOverlap{
        iuid frag_b_id;
        float score;
        int in_degree;
        fragment_end_type bend;                
        short ahang;
        short bhang;
    };

    struct BestFragmentOverlap{

        // Contains information on what a known fragment overlaps.
        // It is assumed that an index into an array of BestOverlap
        // will tell us what fragment has this best overlap

        BestEdgeOverlap five_prime; 
        BestEdgeOverlap three_prime;
    };


    struct BestContainment{

        // Contains what kind of containment relationship exists between
        // fragment a and fragment b

        iuid   container;
        float  contain_score;

        short  a_hang;
        short  b_hang;
        bool   sameOrientation;
        bool   isPlaced;

        bool overlapsAreSorted;
        std::vector<iuid> overlaps;
    };

    //#warning why a map and not a vector or just an array
    typedef std::map<iuid, BestContainment> BestContainmentMap;


    struct BestOverlapGraph {

        friend class BOG_Runner;

        // Constructor, parametrizing maximum number of overlaps
        BestOverlapGraph(double AS_UTG_ERROR_RATE);

        // Destructor
        ~BestOverlapGraph();

        // Interface to graph visitor
        // accept(BestOverlapGraphVisitor bog_vis){
        // bog_vis.visit(this);
        // }

        // Accessor Functions
        BestEdgeOverlap *getBestEdgeOverlap(iuid frag_id, fragment_end_type which_end);
        BestEdgeOverlap *getBestEdgeOverlap(FragmentEnd*);
        // given a FragmentEnd sets it to the next FragmentEnd after following the
        // best edge
        void followOverlap(FragmentEnd*);
        void setBestEdgeOverlap(const OVSoverlap& olap, float newScore);
        void setBestContainer(const OVSoverlap& olap, float newScore);
        iuid getNumFragments() { return lastFrg; }

        bool isContained(const iuid fragid) {
            return(_best_containments.find(fragid) != _best_containments.end());
        };

       // Given a containee, returns pointer to BestContainment record
        BestContainment *getBestContainer(const iuid fragid) {
            return((isContained(fragid)) ? &_best_containments[fragid] : NULL);
        };

        void addContainEdge( iuid, iuid);
        bool containHaveEdgeTo( iuid, iuid);

        // Graph building methods
        fragment_end_type AEnd(const OVSoverlap& olap);
        fragment_end_type BEnd(const OVSoverlap& olap);
        void processOverlap(const OVSoverlap& olap);

        float scoreOverlap(const OVSoverlap& olap) {

#if 0
          // Computes the score for a Error Rate BOG based on overlap
          // corrected error rate.
          //
          // Error rate is normalized so that the higher the error
          // rate, the lower the score.
          //
          return(100.0 - AS_OVS_decodeQuality(olap.dat.ovl.corr_erate) * 100.0);
#endif

#if 0
          // Computes the score for a Longest Edge BOG based on
          // overlap length only.
          //
          return(olapLength(olap));
#endif

#if 0
          // The log for this is:
          //   Add alternate scoring schema that wasn't tested much for posterity's sake.
          //
          if (olap.dat.ovl.orig_erate > consensusCutoff)
              return 0;
          if (olap.dat.ovl.corr_erate > mismatchCutoff)
              return 0;
          return(olapLength(olap) / sqrt(1 + olap.dat.ovl.corr_erate));
#endif

#if 0
          // Computes the score for a Longest Edge BOG based on
          // overlap length but after applying an an error rate
          // cutoff.
          //
          if (olap.dat.ovl.orig_erate > consensusCutoff)
            return 0;
          if (olap.dat.ovl.corr_erate > mismatchCutoff)
            return 0;
          return olapLength(olap);
#endif

#if 1
          //  BPW's new score
          if (olap.dat.ovl.orig_erate > consensusCutoff)
              return 0;
          if (olap.dat.ovl.corr_erate > mismatchCutoff)
              return 0;

          int a_hang = olap.dat.ovl.a_hang;
          int b_hang = olap.dat.ovl.b_hang;

          int alen = fragLen(olap.a_iid);
          int blen = fragLen(olap.b_iid);

          //  anti-normal dovetail; score is length of overlap.
          if ((a_hang < 0) && (b_hang < 0))
              return(MIN((alen + b_hang), (blen + a_hang)));
          
          //  normal dovetail; score is length of overlap.
          if ((a_hang >= 0) && (b_hang >= 0))
              return(MIN((alen - a_hang), (blen - b_hang)));

          //  A contains B
          if ((a_hang < 0) && (b_hang >= 0))
              return(blen * (1.0 - AS_OVS_decodeQuality(olap.dat.ovl.corr_erate)));

          //  B contains A
          if ((a_hang >= 0) && (b_hang < 0))
              return(alen * (1.0 - AS_OVS_decodeQuality(olap.dat.ovl.corr_erate)));
#endif
        };

        // FragStore related variables
        //These should be moved to protected
        static uint16 *fragLength;

        static uint16 fragLen( iuid iid ) {
            return fragLength[ iid ];
        };

        static uint16 olapLength(iuid a_iid, iuid b_iid, short a_hang, short b_hang) {
            int alen = fragLen(a_iid);
            int blen = fragLen(b_iid);

            if (a_hang < 0) {
                if (b_hang < 0 )
                    return alen + b_hang;
                else
                    return blen + a_hang - b_hang; // spur or containment
            } else {
                if (b_hang < 0 )
                    return alen + b_hang - a_hang; // spur or containment
                else
                    return alen - a_hang;
            }
        };
        static uint16 olapLength(const OVSoverlap& olap) {
            return(olapLength(olap.a_iid, olap.b_iid, olap.dat.ovl.a_hang, olap.dat.ovl.b_hang));
        };

        BestContainmentMap _best_containments;

        bool checkForNextFrag(const OVSoverlap& olap);
        void scoreContainment(const OVSoverlap& olap);
        void scoreEdge(const OVSoverlap& olap);
        void updateInDegree(void);
        void removeTransitiveContainment();

    protected:
        static iuid lastFrg;
        iuid curFrag;
        int bestLength;
        BestFragmentOverlap* _best_overlaps;

    public:
        uint64 mismatchCutoff;
        uint64 consensusCutoff;

    }; //BestOverlapGraph


    struct BOG_Runner {
        BOG_Runner(int lastFrag) { BestOverlapGraph::lastFrg = lastFrag; }
        void push_back(BestOverlapGraph *bog) { metrics.push_back(bog); }
        int  size() { return metrics.size(); }
        void processOverlapStream(OverlapStore *, const char*);

        std::vector<BestOverlapGraph *> metrics;
    };
}; //AS_BOG namespace


#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH
