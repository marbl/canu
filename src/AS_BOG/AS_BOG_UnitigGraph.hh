
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
* Module:  AS_BOG_UnitigGraph.hh
* Description:
*	Data structure to contain the unitig paths and how they connect to each
*	other
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
 * $Id: AS_BOG_UnitigGraph.hh,v 1.32 2007-10-17 14:35:06 eliv Exp $
 * $Revision: 1.32 $
*/


#ifndef INCLUDE_AS_BOG_UNITIGGRAPH
#define INCLUDE_AS_BOG_UNITIGGRAPH

static char AS_BOG_UNITIG_GRAPH_HH_CM_ID[] = "$Id: AS_BOG_UnitigGraph.hh,v 1.32 2007-10-17 14:35:06 eliv Exp $";

#include <set>
#include <iostream>
#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_ChunkGraph.hh"

namespace AS_BOG{

	///////////////////////////////////////////////////////////////////////
	// For the sake of self-documenting code

	typedef iuid container_id;
	typedef iuid containee_id;
	typedef iuid fragment_id;
	typedef std::vector<fragment_id> FragmentList;
	typedef std::map<iuid, FragmentList> FragmentEdgeList;

	///////////////////////////////////////////////////////////////////////

    typedef IntMultiPos DoveTailNode;
	typedef std::vector<DoveTailNode> DoveTailPath;
	typedef DoveTailPath::iterator    DoveTailIter;
	typedef DoveTailPath::const_iterator      DoveTailConstIter;
	typedef DoveTailPath::reverse_iterator    DoveTailRIter;

	typedef std::vector<iuid> ContaineeList;	
	typedef std::map<container_id, ContaineeList> ContainerMap;

	///////////////////////////////////////////////////////////////////////

	typedef std::map<fragment_id, SeqInterval> FragmentPositionMap;
	std::ostream& operator<< (std::ostream& os, FragmentPositionMap *fpm_ptr);

    inline bool isReverse( SeqInterval pos ) {
        return(pos.bgn > pos.end);
    }

    struct UnitigBreakPoint {
        int fragNumber;       // the number of the fragment in the unitig
        int inSize;           // the size of the incoming unitig
        int inFrags;          // the number of fragments in incoming unitig
        SeqInterval position; // coordinates in unitig
        FragmentEnd fragEnd;  // frag id and which end to break on 

        UnitigBreakPoint(iuid id=0, fragment_end_type end=FIVE_PRIME) :
            fragEnd(id, end), fragNumber(0), inSize(0), inFrags(0)
        {
            position.bgn = position.end = 0;
        }
    };
	typedef std::list<UnitigBreakPoint> FragmentEnds;

	///////////////////////////////////////////////////////////////////////

	struct Unitig{

		Unitig(iuid =nextId++);		
		~Unitig(void);		

        // Sort frags by position on the unitig
        void sort();

		// Compute unitig based on given dovetails and containments
		void recomputeFragmentPositions(ContainerMap *,BestContainmentMap*, BestOverlapGraph*);
		void computeFragmentPositions(BestOverlapGraph*);

        void shiftCoordinates(int);
        void reverseComplement( );
        void reverseComplement( int offset, BestOverlapGraph *);
		// Accessor methods
		float getAvgRho(void);
		static void setGlobalArrivalRate(float global_arrival_rate);
		void setLocalArrivalRate(float local_arrival_rate);
		float getLocalArrivalRate(void);
		float getCovStat(void);
		long getLength(void);
		long getNumFrags(void);
		long getSumFragLength(void);
		long getNumRandomFrags(void); // For now, same as numFrags, but should be randomly sampled frag count
        DoveTailNode getLastBackboneNode(iuid&);
		// For proto I/O messages
		IntUnitigMesg *getIUM_Mesg();
		void freeIUM_Mesg(IntUnitigMesg *ium_ptr);

		friend std::ostream& operator<< (std::ostream& os, Unitig& utg);

        iuid id();
        static iuid getNextId();
        static void setNextId(iuid);

        void addFrag( DoveTailNode );
        static iuid fragIn(iuid);
        static void resetFragUnitigMap(iuid numFrags);

		// Public Member Variables
        DoveTailPath *dovetail_path_ptr;

		private:
            void placeContains( const ContainerMap *, BestContainmentMap*,
                            const iuid , const SeqInterval, const int level );

			// Do not access these private variables directly, they may not be
			//  computed yet, use accessors!
			float _avgRho;
			float _covStat;
			long _length;
			long _numFrags;
			long _numRandomFrags;
			float _localArrivalRate;
            static iuid nextId;
            iuid _id;
			static float _globalArrivalRate;
            // records if that frg has be incorporated into a unitg or not
            static iuid *_inUnitig;

	};

	///////////////////////////////////////////////////////////////////////

	struct UnitigOverlap{
		// Contains the description of how unitigs are connected to
		// each other
			
		iuid unitig_a;
		orientation_type ori_a;

		iuid unitig_b;
		orientation_type ori_b;

		overlap_type ovl_type;
		int ovl_score;
	};

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	typedef std::vector<Unitig*> UnitigVector;
	typedef UnitigVector::iterator UnitigsIter;
	typedef UnitigVector::const_iterator UnitigsConstIter;
	typedef UnitigVector::reverse_iterator UnitigsRIter;

	struct UnitigGraph{
		// This will store the entire set of unitigs that are generated
		// It's just a unitig container.
        UnitigGraph( BestOverlapGraph*);
        ~UnitigGraph();

		// Call this on a chunk graph pointer to build a unitig graph
		void build(ChunkGraph *cg_ptr);

		// Debugging output operator
		friend std::ostream& operator<< (std::ostream& os, UnitigGraph& utgrph);

		// For compatibility with the rest of the assembler
		void writeIUMtoFile(char *filename);

		// For compatibility with the rest of the assembler
		void readIUMsFromFile(const char *filename);

		float getGlobalArrivalRate(long total_random_frags_in_genome=0, long genome_size=0);

        void breakUnitigs();
        void printUnitigBreaks();

        void filterBreakPoints( Unitig *, FragmentEnds &);

        UnitigBreakPoint selectSmall(const Unitig *,
                                     const FragmentEnds &,
                                     const UnitigBreakPoint &);

        UnitigVector* breakUnitigAt( Unitig *, FragmentEnds &);
        // Single call point for breakUnitigAt(), saves new unitigs in last arg
        void accumulateSplitUnitigs( UnitigsIter, FragmentEnds *, UnitigVector *);

		///////////////////////////////////////////////////////////////
		// Member Variables

		// Unitigs are the dove tails and their contained fragments
		UnitigVector *unitigs;

		// Overlaps are unitig overlaps
		std::vector<UnitigOverlap*> overlaps;

        BestOverlapGraph *bog_ptr;

		private:
			// Given a fragment, it will follow it's overlaps until 
			//   the end, and add them to the unitig
			void populateUnitig(
				Unitig* unitig, 
				iuid src_frag_id, 
				fragment_end_type whichEnd,
				ChunkGraph *cg_ptr,
                int offset);

			// Inverts the containment map to key by container, instead of containee
			void _build_container_map(BestContainmentMap*);

			// Build containee list
			ContainerMap *_extract_containees(DoveTailPath *dtp_ptr, 
				ContainerMap *cntnrmap_ptr);
            // Merge Unitigs by follow the bid's
            void mergeAllUnitigs( std::map<iuid,iuid> *);
            void mergeUnitigs(Unitig*, std::set<iuid> *, std::map<iuid,iuid> *);

            // handles next iteration setup for mergUnitigs
            BestEdgeOverlap* nextJoiner( Unitig*,
                    iuid &aPrev, iuid &fragA, int &tigEnd, bool &begRev,
                    BestEdgeOverlap *&fivePrime, BestEdgeOverlap *&threePrime );

			// Compute the global arrival rate based on the unitig rho's.
			float _compute_global_arrival_rate(void);

            FragmentEdgeList *unitigIntersect;
			ContainerMap *cntnrmap_ptr;
	};
		
	///////////////////////////////////////////////////////////////////////

    struct BogOptions {
        static bool unitigIntersectBreaking;
        static bool ejectUnhappyContained;
    };

};

#endif

