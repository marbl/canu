
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
/*************************************************
 * Module:  AS_BOG_BestOverlapGraph.hh
 * Description:
 *        Data structure to contain the best overlaps and containments
 *        based on a defined metric.
 * 
 *    Programmer:  K. Li
 *       Started:  20 July 2005
 * 
 * Assumptions:
 * 
 * Notes:
 *
 *************************************************/

/* RCS info
 * $Id: AS_BOG_BestOverlapGraph.hh,v 1.35 2007-12-05 23:54:42 brianwalenz Exp $
 * $Revision: 1.35 $
 */

//  System include files

#ifndef INCLUDE_AS_BOG_BESTOVERLAPGRAPH
#define INCLUDE_AS_BOG_BESTOVERLAPGRAPH

#include <set>

#include "AS_BOG_Datatypes.hh"

extern "C" {
#include "AS_OVS_overlapStore.h"
}

namespace AS_BOG{

    ///////////////////////////////////////////////////////////////////////

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

    ///////////////////////////////////////////////////////////////////////

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

    ///////////////////////////////////////////////////////////////////////

    struct BestOverlapGraph {

        friend class BOG_Runner;

        // Constructor, parametrizing maximum number of overlaps
        BestOverlapGraph();

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
        bool isContained(const iuid);
        BestContainment *getBestContainer(iuid frag_id);
        void addContainEdge( iuid, iuid);
        bool containHaveEdgeTo( iuid, iuid);
        void printFrom(iuid begin, iuid end=0);

        // Graph building methods
        fragment_end_type AEnd(const OVSoverlap& olap);
        fragment_end_type BEnd(const OVSoverlap& olap);
        void processOverlap(const OVSoverlap& olap);
        static overlap_type getType(const OVSoverlap & olap);
        virtual float scoreOverlap(const OVSoverlap& olap)=0;
        virtual int getThreshold()=0;

        // FragStore related variables
        //These should be moved to protected
        static uint16 *fragLength;
        static uint16 fragLen( iuid );
        static uint16 olapLength( iuid, iuid, short, short);
        static uint16 olapLength(const OVSoverlap& olap);

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

    }; //BestOverlapGraph

    ///////////////////////////////////////////////////////////////////////////

    struct ErateScore : public BestOverlapGraph {
        ErateScore() : BestOverlapGraph() {}
        float scoreOverlap( const OVSoverlap& olap);
    };

    struct LongestEdge : public BestOverlapGraph {
        LongestEdge() : BestOverlapGraph() {}
        float scoreOverlap( const OVSoverlap& olap);
    };

    struct LongestHighIdent : public BestOverlapGraph {
        int mismatchCutoff;
        int consensusCutoff;
        LongestHighIdent( float maxMismatch) : BestOverlapGraph() {
            mismatchCutoff  = AS_OVS_encodeQuality( maxMismatch / 100.0 );
            consensusCutoff = AS_OVS_encodeQuality( AS_CNS_ERROR_RATE );
            assert( consensusCutoff >= 0 ); // Set in AS_configure
        }
        int getThreshold() { return mismatchCutoff; }
        float scoreOverlap( const OVSoverlap& olap);
    };

    ///////////////////////////////////////////////////////////////////////////
    struct BOG_Runner {
        BOG_Runner(int lastFrag) { BestOverlapGraph::lastFrg = lastFrag; }
        void push_back(BestOverlapGraph *bog) { metrics.push_back(bog); }
        int  size() { return metrics.size(); }
        void processOverlapStream(OverlapStore *, const char*);

        std::vector<BestOverlapGraph *> metrics;
    };
}; //AS_BOG namespace


#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH
