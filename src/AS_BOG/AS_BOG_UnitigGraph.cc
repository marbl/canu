
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
* Module:  AS_BOG_UnitigGraph.cc
* Description:
*	Data structure to contain the unitig paths and how they connect to each
*	other
* 
*    Programmer:  K. Li
*       Started:  2 Aug 2005
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: AS_BOG_UnitigGraph.cc,v 1.63 2007-10-11 19:41:28 eliv Exp $
 * $Revision: 1.63 $
*/

//static char AS_BOG_UNITIG_GRAPH_CC_CM_ID[] = "$Id: AS_BOG_UnitigGraph.cc,v 1.63 2007-10-11 19:41:28 eliv Exp $";
static char AS_BOG_UNITIG_GRAPH_CC_CM_ID[] = "gen> @@ [0,0]";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include <set>
#include <limits>
#include <cstdio>
#include <cmath>
#include <cassert>

extern "C" {
	#include "AS_global.h"
}
#undef max
namespace AS_BOG{

	//////////////////////////////////////////////////////////////////////////////
    UnitigGraph::UnitigGraph(BestOverlapGraph *in_bog_ptr) : bog_ptr(in_bog_ptr)
    {
        unitigs = new UnitigVector;
        unitigIntersect = new FragmentEdgeList;
    }

    UnitigGraph::~UnitigGraph() {
		UnitigVector::iterator utg_itr;
		for(utg_itr=unitigs->begin(); utg_itr!=unitigs->end(); utg_itr++)
            delete *utg_itr;
        delete unitigs;
        delete unitigIntersect;
    }

    // various class static methods and variables
    static std::map<iuid,int>* containPartialOrder;
    static int maxContainDepth = 0;
    iuid Unitig::nextId        = 1;
    iuid Unitig::getNextId()           { return nextId; }
    void Unitig::setNextId(iuid newId) { nextId = newId; }
    iuid Unitig::id()                  { return _id; }
    iuid* Unitig::_inUnitig = NULL;
    void Unitig::resetFragUnitigMap(iuid num_frags) {
        if (_inUnitig == NULL)
            _inUnitig = new iuid[num_frags+1];
        memset( _inUnitig, 0, (num_frags+1)*sizeof(iuid));
    }
    iuid Unitig::fragIn(iuid fragId){
        if (_inUnitig == NULL)
            return 0;
        return _inUnitig[fragId];
    }

    bool BogOptions::unitigIntersectBreaking = false;
    bool BogOptions::ejectUnhappyContained   = false;

	void UnitigGraph::build(ChunkGraph *cg_ptr) {

		iuid num_frags=cg_ptr->getNumFragments();
        Unitig::setNextId( 1 );
		iuid unitig_id=0;

		iuid frag_idx;
		iuid fp_dst_frag_id, tp_dst_frag_id;

		// Initialize where we've been to nowhere; "Do not retraverse list"
        Unitig::resetFragUnitigMap( num_frags );

		BestContainmentMap *best_cntr = &(bog_ptr->_best_containments);
		_build_container_map(best_cntr);

		// Step through all the fragments 
		std::cerr << "Building Unitigs from " << num_frags << " fragments.\n"; 
        while( frag_idx = cg_ptr->nextFragByChunkLength() ) {
            if (BestOverlapGraph::fragLen(frag_idx) == std::numeric_limits<uint16>::max())
                continue; // Deleted frag

            // Check the map to so we don't visit a unitig twice (once from
            //   both ends)
            if( !Unitig::fragIn( frag_idx ) && 
                    best_cntr->find(frag_idx) == best_cntr->end() )
            { 
                cg_ptr->getChunking(
                        frag_idx, 
                        fp_dst_frag_id, 
                        tp_dst_frag_id);

                if(!(unitig_id%10000)){
                    std::cerr << ".";
                }

                // Allocated a new unitig node
                Unitig *utg=new Unitig;
                unitig_id = utg->id();

                // create it going off the 5' end
                populateUnitig( utg, frag_idx, FIVE_PRIME, cg_ptr, 0);
                // now check if we can also go off 3' end
                fragment_end_type whichEnd = THREE_PRIME;
                BestEdgeOverlap *tpBest = bog_ptr->getBestEdgeOverlap( frag_idx, whichEnd);
                iuid tpId = tpBest->frag_b_id;
                // if the other end is already in another unitig, don't try to extend
                if (Unitig::fragIn(tpId)) {
                    // question, should we include the circle in the intersection list?
                    if (tpId != 0 && unitig_id != Unitig::fragIn(tpId))
                    {
                      // save this outgoing intersection point for future use
                      (*unitigIntersect)[tpId].push_back( frag_idx );
                      fprintf(stderr,"Unitig %5d 1st frag %7d -> Unitig %5d frag %7d\n",
                            unitig_id, frag_idx, Unitig::fragIn(tpId), tpId );
                    }
                } else {
                    // if other end is also 3' we need to walk of it's 5'
                    if (tpBest->bend == whichEnd)
                        whichEnd = FIVE_PRIME;

                    int len = utg->getLength();
                    int flen = BestOverlapGraph::fragLen(frag_idx);
                    int ahng = tpBest->ahang;
                    int offset = len - flen + ahng;
                    fprintf(stderr,
                        "Join unitig %d len %d at %d len %d and %d ahang %d offset %d\n",
                         unitig_id, len, frag_idx, flen, tpBest->frag_b_id, ahng, offset );

                    utg->reverseComplement();
                    populateUnitig( utg,
                            tpBest->frag_b_id, whichEnd, cg_ptr, offset);
                }

                // Store unitig in unitig graph
                unitigs->push_back(utg);

            }else{
                // We are either in the middle of a dovetail sequence,
                // it is a singleton, or we are going from 3' to 5'
                // across a dovetail.

            }
		}
        std::cerr << std::endl;

        // Pick up frags missed above, possibly circular unitigs
        for(frag_idx=1; frag_idx<=num_frags; frag_idx++){
            if (Unitig::fragIn( frag_idx ) == 0) {
				cg_ptr->getChunking( frag_idx, fp_dst_frag_id, tp_dst_frag_id);

                if ( fp_dst_frag_id != NULL_FRAG_ID &&
                     tp_dst_frag_id != NULL_FRAG_ID )
                {
					Unitig *utg=new Unitig;
                    // need to make populateUnitig break circle
                    populateUnitig( utg, frag_idx, FIVE_PRIME, cg_ptr, 0 );

                    fprintf(stderr,"Circular unitig %d 1st frag %d\n",
                            unitig_id,frag_idx);

					unitigs->push_back(utg);

                } else {
                    // Should both be null or neither null
                    // otherwise main loop failed to find it
                    assert( fp_dst_frag_id == NULL_FRAG_ID );
                    assert( tp_dst_frag_id == NULL_FRAG_ID );
                }
            }
        }
        if (BogOptions::unitigIntersectBreaking) { 
            printUnitigBreaks();
            breakUnitigs();
        }

        UnitigsConstIter uIter = unitigs->begin();
        for (; uIter != unitigs->end(); uIter++ ) {
            Unitig* tig = *uIter;
            if (tig == NULL)
                continue;
            tig->recomputeFragmentPositions(cntnrmap_ptr, best_cntr, bog_ptr);
        }
        fprintf(stderr,"Max contain depth is %d\n",maxContainDepth);
        delete cntnrmap_ptr;
        cntnrmap_ptr = NULL;
	}
    //////////////////////////////////////////////////////////////////////////////
    BestEdgeOverlap* UnitigGraph::nextJoiner( Unitig* tig, 
            iuid &aPrev, iuid &fragA, int &tigEnd, bool &begRev,
                       BestEdgeOverlap *&fivePrime, BestEdgeOverlap *&threePrime )
    {
        aPrev = 0;
        DoveTailNode joinerNode = tig->getLastBackboneNode(aPrev);
        fragA = joinerNode.ident;
        assert( aPrev != fragA );
        // Never join on a contain
        assert( joinerNode.contained == 0 );
        fivePrime = bog_ptr->getBestEdgeOverlap( fragA, FIVE_PRIME );
        threePrime = bog_ptr->getBestEdgeOverlap( fragA, THREE_PRIME );

        BestEdgeOverlap* bestEdge;
        begRev = joinerNode.position.bgn > joinerNode.position.end;
        tigEnd = begRev ? joinerNode.position.bgn : joinerNode.position.end;
        if ( aPrev == fivePrime->frag_b_id ) {
            bestEdge = threePrime;
        } else if ( aPrev == threePrime->frag_b_id ) {
            bestEdge = fivePrime;
        } else {
            fprintf(stderr,"A Disagree: A %d PREV %d 5' %d,%d 3' %d,%d\n",
                    fragA, aPrev, 
                    fivePrime->frag_b_id, fivePrime->in_degree,
                    threePrime->frag_b_id, threePrime->in_degree
            );
            if ( begRev )
                bestEdge = fivePrime;
            else
                bestEdge = threePrime;
        }
        return bestEdge;
    }

    //////////////////////////////////////////////////////////////////////////////
    void UnitigGraph::mergeUnitigs(Unitig* tig, std::set<iuid> *joined,
            std::map<iuid,iuid> *visited_map)
    {
        iuid beforeLast, beginId, joiner;
        int tigEnd;
        bool begRev;
        BestEdgeOverlap *bestEdge,*fivePrime,*threePrime;
        bestEdge = nextJoiner( tig, beforeLast,
                beginId, tigEnd, begRev, fivePrime, threePrime
        );
        joiner = bestEdge->frag_b_id;
        iuid tigIdToAdd = (*visited_map)[ joiner ];
        while (joiner != 0 && joined->find(tigIdToAdd) == joined->end() &&
                 visited_map->find(joiner) != visited_map->end())
        {
            Unitig* tigToAdd = unitigs->at( tigIdToAdd - 1);
            // Code assumes joining only from end of unitig
            if ( beforeLast != 0 ) // ok to join singleton though
                assert( beginId != tig->dovetail_path_ptr->front().ident);
            BestEdgeOverlap* addFivePrime = bog_ptr->getBestEdgeOverlap(
                    joiner, FIVE_PRIME );
            BestEdgeOverlap* addThreePrime = bog_ptr->getBestEdgeOverlap(
                    joiner, THREE_PRIME );

            bool bDisagree = false;
            if ( bestEdge->bend == FIVE_PRIME  && beginId != addFivePrime->frag_b_id
              || bestEdge->bend == THREE_PRIME && beginId != addThreePrime->frag_b_id)
            { 
                fprintf(stderr,"B Disagree: A %d,%d B %d 5' %d,%d 3' %d,%d\n",
                    beginId, bestEdge->in_degree, joiner,
                    addFivePrime->frag_b_id, addFivePrime->in_degree,
                    addThreePrime->frag_b_id, addThreePrime->in_degree
                );
                 bDisagree = true;
            }

            bool reverse = false;
            DoveTailNode first = tigToAdd->dovetail_path_ptr->front();
            iuid notUsed;
            DoveTailNode last  = tigToAdd->getLastBackboneNode(notUsed);
            if (joiner == first.ident)
                last = first;
            bool lastReverse = last.position.bgn > last.position.end;
            // asserts one of the in_degree's should always be none zero at a unitig break
            // then make sure the both frags agree on which ends are being used
            if ( bestEdge == threePrime ) {
                if ( bestEdge->bend == FIVE_PRIME ) {
                    if (begRev ^ lastReverse)
                        reverse = true;
                    fprintf(stderr,"3'->5' %d %d bR %d lR %d R %d\n", beginId, joiner,
                                            begRev, lastReverse, reverse);
                    assert( bDisagree || bestEdge->in_degree != 1 ||
                            addFivePrime->in_degree != 1);
                } else if ( bestEdge->bend == THREE_PRIME ) {
                    if (!(begRev ^ lastReverse))
                        reverse = true;
                    fprintf(stderr,"3'->3' %d %d bR %d lR %d R %d\n", beginId, joiner,
                                            begRev, lastReverse, reverse);
                    assert( bDisagree || bestEdge->in_degree != 1 ||
                                addThreePrime->in_degree != 1);
                } else {
                    assert(0);
                }
            } else {
                if ( bestEdge->bend == FIVE_PRIME ) {
                    if (!(begRev ^ lastReverse))
                        reverse = true;
                    fprintf(stderr,"5'->5' %d %d bR %d lR %d R %d\n", beginId, joiner,
                                            begRev, lastReverse, reverse);
                    assert( bDisagree || bestEdge->in_degree != 1 ||
                                addFivePrime->in_degree != 1 );
                } else if ( bestEdge->bend == THREE_PRIME ) {
                    if (begRev ^ lastReverse)
                        reverse = true;
                    fprintf(stderr,"5'->3' %d %d bR %d lR %d R %d\n", beginId, joiner,
                                            begRev, lastReverse, reverse);
                    assert( bDisagree || bestEdge->in_degree != 1 ||
                                addThreePrime->in_degree != 1);
                } else {
                    assert(0);
                }
                if (!begRev) {
                    if ( beforeLast == 0 ) // skip singleton, doesn't add much
                        return;
                    fprintf(stderr, "Begin needs reverse: %d to %d.\n",beginId,joiner);
                    // should only ever try to join this end on a singleton
                    assert( first.ident == last.ident);
                    return; // need to reverse the first tig before joinging
                }
            }
            // if the start frag is reversed, we join on it's 5' edge, else 3'
            if ( begRev )
                assert( bestEdge == fivePrime );
            else
                assert( bestEdge == threePrime );

            int offset  = 1 + tigEnd - BestOverlapGraph::olapLength( beginId,
                    joiner, bestEdge->ahang, bestEdge->bhang
            );
            if (!reverse && joiner == first.ident) {
                DoveTailPath::iterator addIter = tigToAdd->dovetail_path_ptr->begin();
                for(;addIter != tigToAdd->dovetail_path_ptr->end(); addIter++) {
                    addIter->position.bgn += offset;
                    addIter->position.end += offset;
                    tig->addFrag( *addIter );
                }
            } else { // reverse complement and append
                tigToAdd->reverseComplement( offset, bog_ptr);
                tig->dovetail_path_ptr->insert( tig->dovetail_path_ptr->end(),
                                           tigToAdd->dovetail_path_ptr->begin(),
                                           tigToAdd->dovetail_path_ptr->end() );
            }
            delete tigToAdd;
            unitigs->at( tigIdToAdd - 1) = NULL;
            joined->insert(tigIdToAdd);
            bestEdge = nextJoiner( tig, beforeLast,
                    beginId, tigEnd, begRev, fivePrime, threePrime
            );
            joiner = bestEdge->frag_b_id;
            if ( joiner == 0)
                break;
            tigIdToAdd = (*visited_map)[ joiner ];
        }
    }
    //////////////////////////////////////////////////////////////////////////////
    void UnitigGraph::mergeAllUnitigs( std::map<iuid,iuid> *visited_map) {

        std::set<iuid>* joined = new std::set<iuid>;
        UnitigVector* newTigs = new UnitigVector;
        UnitigVector::iterator tigIter = unitigs->begin();
        for(;tigIter != unitigs->end(); tigIter++) {
            Unitig* tig = *tigIter;
            if (tig == NULL)
                continue;
            if (joined->find( tig->id() ) == joined->end()) {
                joined->insert( tig->id() );
                if (tig->getNumFrags() == 1) {
                    newTigs->push_back( tig );
                    continue;
                }
                mergeUnitigs( tig, joined, visited_map );
                newTigs->push_back( tig );
            }
        }
        delete unitigs;
        unitigs = newTigs;
        delete joined;
    }
	//////////////////////////////////////////////////////////////////////////////
    DoveTailNode Unitig::getLastBackboneNode(iuid &prevId) {
        DoveTailPath::reverse_iterator rIter = dovetail_path_ptr->rbegin();
        DoveTailNode lastNonContain;
        memset(&lastNonContain, 0, sizeof(lastNonContain));
        for(;rIter != dovetail_path_ptr->rend(); rIter++) {
            DoveTailNode node = *rIter;
            if (node.contained == 0) {
                if (lastNonContain.ident == 0)
                    lastNonContain = node;
                else
                {
                    prevId = node.ident;
                    return lastNonContain;
                }
            }
        }
        return lastNonContain;
    }
	//////////////////////////////////////////////////////////////////////////////

    void Unitig::computeFragmentPositions(BestOverlapGraph* bog_ptr) {
        if (dovetail_path_ptr == NULL || dovetail_path_ptr->empty())
            return;
        // we need to determine which orientation the first frag is in
        // so we peek ahead to the 2nd frag to find out
        DoveTailIter iter = dovetail_path_ptr->begin();
        int frag_begin,fragPrevEnd,fragNextEnd;
        frag_begin = fragPrevEnd = fragNextEnd = 0;
        iuid first = iter->ident;
        int frag_end = BestOverlapGraph::fragLen(first);
        fragPrevEnd = frag_end;
        if(  dovetail_path_ptr->size() == 1) {
            // still set singleton's coords
            dovetail_path_ptr->front().position.bgn = frag_begin;
            dovetail_path_ptr->front().position.end = frag_end;
            return;
        }
        BestEdgeOverlap* bestEdge = bog_ptr->getBestEdgeOverlap(
                first, FIVE_PRIME );
        BestEdgeOverlap* threeP = bog_ptr->getBestEdgeOverlap(
                first, THREE_PRIME );
        fprintf(stderr,"1stFrag %d beg %d end %d\n",
                    first, frag_begin, frag_end );
        iter++;
        if (iter->ident == bestEdge->frag_b_id) {
            dovetail_path_ptr->front().position.bgn = frag_end;
            dovetail_path_ptr->front().position.end = frag_begin;
            fprintf(stderr,"Go off 5' edge and get frag %d\n",iter->ident);
        } else {
            // doesn't handle contianed, only backbone
            fprintf(stderr,"Go off 3' edge and get frag %d %d\n",iter->ident,threeP->frag_b_id);
            assert( iter->ident == threeP->frag_b_id );
            bestEdge = threeP;
            dovetail_path_ptr->front().position.bgn = frag_begin;
            dovetail_path_ptr->front().position.end = frag_end;
        }
        fprintf(stderr,"ahang %d bhang %d\n", bestEdge->ahang, bestEdge->bhang);
        if ( bestEdge->ahang > 0 )
            frag_begin = bestEdge->ahang;
        else 
            frag_begin = -bestEdge->bhang;

        fragment_end_type whichEnd = FIVE_PRIME; // the end to walk off for next frag
        if (bestEdge->bend == FIVE_PRIME)
            whichEnd = THREE_PRIME;

        for (; iter != dovetail_path_ptr->end(); iter++ ) {
            DoveTailNode* node = &*iter;
            //DoveTailNode* nextNode = &*(iter+1);
			bestEdge = bog_ptr->getBestEdgeOverlap( node->ident, whichEnd);
            int currLen = BestOverlapGraph::fragLen(node->ident);
            // The end of the fragment can be calulated as start + length
            // or as end of previous frag + b_hang. They should be the same
            // if there are no gaps, but with gaps the 2nd method should be better
            frag_end = frag_begin + currLen;

            fprintf(stderr,"Frag %d len %d beg %d end %d ahang %d bhang %d end %s bfrg %d bend %s\n",
                    node->ident, currLen, frag_begin, frag_end,
                     bestEdge->ahang, bestEdge->bhang, whichEnd == FIVE_PRIME ? "5'":"3'",
                     bestEdge->frag_b_id, bestEdge->bend == FIVE_PRIME ? "5'":"3'");

            int end;
            // pick the smallest end that's greater then the previous end
            // this is critical for preserving the correct frag order after
            // a reverse complement in the unitig merging code
            if (frag_end < fragPrevEnd) {
                end = fragNextEnd;
            } else if (fragNextEnd < fragPrevEnd) {
                end = frag_end;
            } else {
                end = fragNextEnd > frag_end ? frag_end : fragNextEnd;
            }
/* only if I can't find a better way to fix up reverse complemented unitigs
           if (nextNode != dovetail_path_ptr->end() &&
                bestEdge->frag_b_id  != nextNode->ident) {
                // disagreement, use next node's end and hangs
                BestEdgeOverlap tP = bog_ptr->getBestEdgeOverlap( nextNode->ident,THREE_PRIME);
                BestEdgeOverlap fP = bog_ptr->getBestEdgeOverlap( nextNode->ident,FIVE_PRIME);
                if (node->ident == fP->frag_b_id) {
                    if (whichEnd == FIVE_PRIME) {
                       if ( fp->bend == FIVE_PRIME) {
                       }
                    }
                }
            }
*/

			if(whichEnd == FIVE_PRIME){
				node->position.bgn = end;
				node->position.end = frag_begin;
			}else {
				node->position.bgn = frag_begin;
				node->position.end = end;
			}
            fragNextEnd = frag_end;
            fragPrevEnd = end;

            // Prep the start position of the next fragment
            if (bestEdge->ahang < 0 && bestEdge->bhang < 0 ) {
                fragNextEnd -= bestEdge->ahang ;
			    frag_begin  -= bestEdge->bhang ;
            } else {
                fragNextEnd += bestEdge->bhang ;
			    frag_begin  += bestEdge->ahang ;
            }
            if ( bestEdge->bend == FIVE_PRIME ) {
                whichEnd = THREE_PRIME;
            } else {
                whichEnd = FIVE_PRIME;
            }
        }
    }
	void UnitigGraph::populateUnitig( Unitig* unitig,
		iuid src_frag_id, fragment_end_type firstEnd, ChunkGraph *cg_ptr, int offset){

	// Note:  I only need BestOverlapGraph for it's frag_len and olap_length

		iuid current_frag_id=src_frag_id;
		iuid next_frag_id = 0;
		fragment_end_type whichEnd = firstEnd;
        int frag_end,frag_begin,fragNextEnd,fragPrevEnd;
        frag_begin = fragNextEnd = fragPrevEnd = offset;

		while(current_frag_id != NULL_FRAG_ID && !Unitig::fragIn( current_frag_id )) {
			// Store the current fragment into dovetail path
			BestEdgeOverlap* bestEdge = bog_ptr->getBestEdgeOverlap(
				current_frag_id, whichEnd
                );
            DoveTailNode dt_node;
            dt_node.type         = AS_READ; /* Isolated fragtype VR */
            dt_node.ident        = current_frag_id;
            dt_node.contained    = 0;
            dt_node.delta_length = 0;
            dt_node.delta        = NULL;
            next_frag_id         = bestEdge->frag_b_id;
#ifdef NEW_UNITIGGER_INTERFACE
            dt_node.ident2       = next_frag_id;
            // consensus wants positive hangs, so swap
            if (bestEdge->ahang < 0 && bestEdge->bhang < 0 ) {
                dt_node.ahang = -bestEdge->bhang;
                dt_node.bhang = -bestEdge->ahang;
            } else {
                dt_node.ahang = bestEdge->ahang;
                dt_node.bhang = bestEdge->bhang;
            }
#endif
            int currLen = BestOverlapGraph::fragLen(current_frag_id);
            // The end of the fragment can be calulated as start + length
            // or as end of previous frag + b_hang. They should be the same
            // if there are no gaps, but with gaps the 2nd method should be better
            frag_end = frag_begin + currLen;
            if ( fragNextEnd == offset )
                fragNextEnd = currLen + offset;

//            fprintf(stderr,"Frag %d len %d beg %d end %d ahang %d bhang %d nEnd %d\n",
//                    current_frag_id, currLen, frag_begin, frag_end,
//                     bestEdge->ahang, bestEdge->bhang, fragNextEnd);

            int end;
            // pick the smallest end that's greater then the previous end
            // this is critical for preserving the correct frag order after
            // a reverse complement in the unitig merging code
            if (frag_end < fragPrevEnd) {
                end = fragNextEnd;
            } else if (fragNextEnd < fragPrevEnd) {
                end = frag_end;
            } else {
                end = fragNextEnd > frag_end ? frag_end : fragNextEnd;
            }
			if(whichEnd == FIVE_PRIME){
				dt_node.position.bgn = end;
				dt_node.position.end = frag_begin;
			}else {
				dt_node.position.bgn = frag_begin;
				dt_node.position.end = end;
			}
            fragNextEnd = frag_end;
            fragPrevEnd = end;
			unitig->addFrag(dt_node);

            // Prep the start position of the next fragment
            if (bestEdge->ahang < 0 && bestEdge->bhang < 0 ) {
                fragNextEnd -= bestEdge->ahang ;
			    frag_begin  -= bestEdge->bhang ;
            } else {
                fragNextEnd += bestEdge->bhang ;
			    frag_begin  += bestEdge->ahang ;
            }
            int chunkNextId = cg_ptr->getChunking(current_frag_id, whichEnd); 
            if ( chunkNextId != NULL_FRAG_ID )
                assert( chunkNextId == next_frag_id );

			// Set current to next
			next_frag_id = current_frag_id;
			current_frag_id = chunkNextId;
            if ( current_frag_id == src_frag_id )
                break; // Break a circle
            if ( bestEdge->bend == FIVE_PRIME ) {
                whichEnd = THREE_PRIME;
            } else {
                whichEnd = FIVE_PRIME;
            }
        }
		if ( next_frag_id != 0 && Unitig::fragIn( current_frag_id ) ) {
            // save this outgoing intersection point for future use
            (*unitigIntersect)[current_frag_id].push_back( next_frag_id );
            fprintf(stderr,"Unitig %5d frag %7d -> Unitig %5d frag %7d\n",
                    unitig->id(), next_frag_id, Unitig::fragIn(current_frag_id),
                    current_frag_id );
        }
	}
	//////////////////////////////////////////////////////////////////////////////
    void UnitigGraph::accumulateSplitUnitigs( UnitigsIter tig, FragmentEnds* breaks,
            UnitigVector* splits)
    {
        UnitigVector* newUs = breakUnitigAt( *tig, *breaks );
        if (newUs == NULL)
            return;
        if (!newUs->empty()) {
            // we have some new split up unitigs
            int frgCnt = 0;
            UnitigsConstIter newIter = newUs->begin();
            // Number of frags should stay the same
            for(;newIter != newUs->end(); newIter++) {
                Unitig* tigTmp = *newIter;
                frgCnt += tigTmp->dovetail_path_ptr->size();
            }
            fprintf(stderr,"Num frgs after splits is %d\n",frgCnt);
            // store new unitigs
            splits->insert( splits->end(), newUs->begin(), newUs->end());
            // remove old unsplit unitig
            delete *tig;
            *tig = NULL;
        }
        delete newUs;
    }
	//////////////////////////////////////////////////////////////////////////////
    void UnitigGraph::breakUnitigs() {

        const char * fiveP  = "5'";
        const char * threeP = "3'";
        UnitigVector* splits = new UnitigVector(); // save split locations
        UnitigsIter  tigIter = unitigs->begin();
        for(;tigIter != unitigs->end(); tigIter++)
        {
            Unitig* tig = *tigIter;
            tig->sort();
            FragmentEnds breaks;
            int fragCount = 1;
            DoveTailNode lastBackbone;
            DoveTailPath::const_iterator dt_itr;
            for( dt_itr  = tig->dovetail_path_ptr->begin();
                 dt_itr != tig->dovetail_path_ptr->end(); dt_itr++,fragCount++)
            {
                DoveTailNode node = *dt_itr;
                iuid dtFrag = node.ident;
                // Include contained frags in count for coverage calculation
                if ( cntnrmap_ptr->find( dtFrag ) != cntnrmap_ptr->end()) {
                    fragCount += (*cntnrmap_ptr)[ dtFrag ].size();
                }
                BestEdgeOverlap *bestEdge;
                // Detect first or last frg in unitig
                if (dt_itr   == tig->dovetail_path_ptr->begin() ||
                    dt_itr+1 == tig->dovetail_path_ptr->end() )
                {
                    // default to 3' end of untig for the last frag
                    fragment_end_type dtEnd = THREE_PRIME;
                    const char *endStr = threeP;
                    // want the backbone, not a contained
                    if (node.contained)
                        node = lastBackbone;

                    if ( isReverse( node.position ) )
                        dtEnd = FIVE_PRIME;

                    // At the begining we want the 5' end of unitig, so switch
                    if (dt_itr == tig->dovetail_path_ptr->begin()) {
                        dtEnd = dtEnd == FIVE_PRIME ? THREE_PRIME : FIVE_PRIME;
                        endStr = fiveP;
                    }

                    bestEdge = bog_ptr->getBestEdgeOverlap(node.ident,dtEnd);
                    iuid bFrg = bestEdge->frag_b_id;
                    fprintf(stderr,"Unitig %5d %s frag %7d points to Unitig %5d frag %7d\n",
                            tig->id(), endStr, node.ident, Unitig::fragIn(bFrg), bFrg);
                }
                // store previous backbone frag for next loop iteration
                if (!node.contained)
                    lastBackbone = node;

                FragmentEdgeList::const_iterator edge_itr =
                                            unitigIntersect->find( dtFrag );
                if (edge_itr != unitigIntersect->end() ) {
                // We have a set of best edges incoming from other unitigs
                    FragmentList::const_iterator fragItr;
                    for( fragItr  = edge_itr->second.begin();
                         fragItr != edge_itr->second.end();  fragItr++)
                    {
                        iuid inFrag = *fragItr;
                        BestEdgeOverlap *bestEdge;
                        // check if it's incoming frag's 5' best edge
                        fragment_end_type bestEnd = FIVE_PRIME;
                        bestEdge = bog_ptr->getBestEdgeOverlap( inFrag, bestEnd );
                        if (bestEdge->frag_b_id != dtFrag)
                        { // not 5' best edge, so 3'
                            bestEnd = THREE_PRIME;
                            bestEdge = bog_ptr->getBestEdgeOverlap( inFrag, bestEnd );
                        }
                        // dtFrag must be either 3' or 5' best edge of inFrag
                        assert( bestEdge->frag_b_id == dtFrag );
                        const char *inEnd = bestEnd == FIVE_PRIME ? fiveP : threeP;
                        // Guess dtFrag's end is 3' or end pos
                        const char *dtEnd = threeP;
                        int pos = node.position.end;
                        if ( bestEdge->bend == FIVE_PRIME) {
                            // nope it's the 5'/begining of dtFrag
                            dtEnd = fiveP;
                            pos = node.position.bgn;
                        }
                        // get incoming unitig, unitig id's start at 1, storage at 0
                        Unitig *inTig = (*unitigs)[Unitig::fragIn(inFrag)-1];
                        if (inTig == NULL) {
                            fprintf(stderr, "  NullBreak tig %5d at frag %7d\n",
                                    tig->id(), dtFrag );
                        }
                        else {
                            UnitigBreakPoint breakPoint(dtFrag, bestEdge->bend);
                            breakPoint.position = node.position;
                            breakPoint.fragNumber = fragCount ;
                            breakPoint.inSize = inTig->getLength();
                            breakPoint.inFrags = inTig->getNumFrags();
                            breaks.push_back( breakPoint );
                        }
                        fprintf(stderr,
"  Utig %4d frgs %5d len %4d frag %6d end %s into Utig %4d frag %6d end %s pos %6d\n",
                                Unitig::fragIn(inFrag), inTig->getNumFrags(),
                                inTig->getLength(), inFrag, inEnd,
                                tig->id(), dtFrag, dtEnd, pos );
                    }
                }
            }
            if ( !breaks.empty() ) {
                // create a final fake bp for the last frag so
                // we have a reference point to the end of the tig for
                // filtering the breakpoints
                fragment_end_type lastEnd = THREE_PRIME;
                if ( isReverse( lastBackbone.position ) )
                        lastEnd = FIVE_PRIME;
                UnitigBreakPoint breakPoint(lastBackbone.ident, lastEnd);
                breakPoint.position   = lastBackbone.position;
                breakPoint.fragNumber = fragCount;
                breakPoint.inSize     = std::numeric_limits<int>::max();
                breaks.push_back( breakPoint );

                accumulateSplitUnitigs( tigIter, &breaks, splits );
            }
        }
        fprintf(stderr,"Num splits %d\n",splits->size());
        unitigs->insert( unitigs->end(), splits->begin(), splits->end());
        delete splits;
    }
	//////////////////////////////////////////////////////////////////////////////
    void UnitigGraph::printUnitigBreaks() {

        UnitigVector::const_iterator tigIter = unitigs->begin();
        for(;tigIter != unitigs->end(); tigIter++) {
            Unitig* tig = *tigIter;
            DoveTailNode first = tig->dovetail_path_ptr->front();
            iuid prev = 0;
            DoveTailNode last = tig->getLastBackboneNode(prev);
            if ( prev == 0 )
                continue; // skip singletons
            iuid id1 = first.ident;
            iuid id2 = last.ident;
            assert( prev != id2 );
            assert( last.contained == 0 );
            BestEdgeOverlap *firstBest, *lastBest;
            if ( first.position.bgn < first.position.end) {
                firstBest = bog_ptr->getBestEdgeOverlap( id1, FIVE_PRIME );
            } else {
                firstBest = bog_ptr->getBestEdgeOverlap( id1, THREE_PRIME );
            }
            if ( last.position.bgn < last.position.end) {
                lastBest = bog_ptr->getBestEdgeOverlap( id2, THREE_PRIME );
            } else {
                lastBest = bog_ptr->getBestEdgeOverlap( id2, FIVE_PRIME );
            }
            if (lastBest->frag_b_id == 0 || firstBest->frag_b_id == 0)
                continue; // Skip non bubbles, anchored ends
            fprintf(stderr, "Bubble %d to %d len %d\n",firstBest->frag_b_id,
                    lastBest->frag_b_id, tig->getLength() );  
            if (firstBest->frag_b_id == lastBest->frag_b_id) {  
                fprintf(stderr, "Self Bubble %d to %d\n",id1,id2);  
            }
        }
    }
	//////////////////////////////////////////////////////////////////////////////

	void UnitigGraph::_build_container_map(BestContainmentMap *best_cntr)
    {
		cntnrmap_ptr = new ContainerMap;
		int containees=0;

		BestContainmentMap::const_iterator bstcnmap_itr;
		for( bstcnmap_itr  = best_cntr->begin(); 
		     bstcnmap_itr != best_cntr->end(); bstcnmap_itr++)
        {
			iuid ctnee_id     = bstcnmap_itr->first;
			iuid container_id = bstcnmap_itr->second.container;

			(*cntnrmap_ptr)[container_id].push_back(ctnee_id);
			containees++;
		}
		
		std::cerr << "BestContainments Inverted into ContainerMap." << std::endl;
		std::cerr << "Num containees: " << containees << std::endl;
		std::cerr << "Num containers: " << cntnrmap_ptr->size() << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////

	std::ostream& operator << (std::ostream& os, Unitig& utg){
		
		os << "Dovetails:" << std::endl;

		DoveTailPath::const_iterator dt_itr;
		for(dt_itr=utg.dovetail_path_ptr->begin();
			dt_itr!=utg.dovetail_path_ptr->end(); 
			dt_itr++){
			
			os << "  " << dt_itr->ident << std::endl;
		}

		os << "Containments:" << std::endl;
	
		os << "avgRho: " << utg.getAvgRho() << std::endl;
		os << "covStat: " << utg.getCovStat() << std::endl;
		os << "numFrags: " << utg.getNumFrags() << std::endl;
		os << "numRandomFrags: " << utg.getNumRandomFrags() << std::endl;
		os << "length: " << utg.getLength() << std::endl;

		return(os);

	}

	//////////////////////////////////////////////////////////////////////////////

	std::ostream& operator << (std::ostream& os, UnitigGraph& utgrph){
		
		UnitigVector::const_iterator utg_itr;
		iuid num_utgs=0;
		for(utg_itr=utgrph.unitigs->begin();
			utg_itr!=utgrph.unitigs->end();
			utg_itr++){
		
			std::cerr << num_utgs << std::endl;
			std::cerr << (**utg_itr) << std::endl;
			num_utgs++;
		
		}

		return(os);
	}

	//////////////////////////////////////////////////////////////////////////////

	float UnitigGraph::getGlobalArrivalRate(long total_random_frags_in_genome, long genome_size){
		
		float _globalArrivalRate;

		//  If the genome size has not been specified, estimate the GAR.
		if(genome_size == 0){

			float total_rho=0, avg_rho;
			float total_arrival_frags=0;
            size_t rho_gt_10000 = 0;

			// Go through all the unitigs to sum rho and unitig arrival frags
			UnitigVector::const_iterator iter;
			for(
			    iter=unitigs->begin();
			    iter!=unitigs->end();
			    iter++){
				
                if (*iter == NULL)
                    continue;

				avg_rho = (*iter)->getAvgRho();
				total_rho += avg_rho;
                if (avg_rho > 10000.0)
                    rho_gt_10000 += (size_t)avg_rho / 10000;
            
				float unitig_random_frags = (*iter)->getNumRandomFrags();
                if (--unitig_random_frags < 0)
                    unitig_random_frags = 0;

				total_arrival_frags += unitig_random_frags;
                (*iter)->setLocalArrivalRate( unitig_random_frags / avg_rho );
			}
			// Estimate GAR
            _globalArrivalRate = (total_rho > 0) ? (total_arrival_frags / total_rho): 0;

            std::cerr << "Calculated Global Arrival rate " << _globalArrivalRate <<
                std::endl;
            // Now recalculate based on big unitigs, copied from AS_CGB/AS_CGB_cgb.c
            if (rho_gt_10000 * 20000 > total_rho) {
                float min_10_local_arrival_rate          = _globalArrivalRate;
                float median_local_arrival_rate          = _globalArrivalRate;
                float max_local_arrival_rate             = _globalArrivalRate;
                float recalibrated_fragment_arrival_rate = _globalArrivalRate;
                size_t num_arrival_rates=0;
                int median_index;
                std::vector<float> arrival_rate_array(rho_gt_10000);
                
                for( iter=unitigs->begin(); iter!=unitigs->end(); iter++)
                {
                    if (*iter == NULL)
                        continue;
                    avg_rho = (*iter)->getAvgRho();
                    if (avg_rho > 10000.0) {
                        const int num_10000 = (size_t)avg_rho / 10000;
                        const float local_arrival_rate =
                            (*iter)->getNumRandomFrags() / avg_rho;
                        assert(num_10000 > 0);
                        int i;
                        for(i=0;i<num_10000;i++){
                            assert(i < rho_gt_10000);
                            arrival_rate_array[i] = local_arrival_rate;
                            num_arrival_rates++;
                        }
                        if(num_arrival_rates > 0){
                            float tmp_fragment_arrival_rate, max_diff_arrival_rate;
                            float prev_arrival_rate, cur_arrival_rate, diff_arrival_rate;
                            int max_diff_index;
                            std::sort(arrival_rate_array.begin(),arrival_rate_array.end());
                            min_10_local_arrival_rate = arrival_rate_array[num_arrival_rates / 10];
                            median_index = (num_arrival_rates * 5) / 10;
                            median_local_arrival_rate = arrival_rate_array[median_index];
                            max_local_arrival_rate = arrival_rate_array[num_arrival_rates-1];
                            recalibrated_fragment_arrival_rate =
                                arrival_rate_array[(num_arrival_rates * 19) / 20];
                            prev_arrival_rate = min_10_local_arrival_rate;
                            max_diff_arrival_rate = 0.0;
                            for(i=num_arrival_rates / 10;i<median_index;i++){
                                cur_arrival_rate = arrival_rate_array[i];
                                diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
                                prev_arrival_rate = cur_arrival_rate;
                                if(diff_arrival_rate > max_diff_arrival_rate){
                                    max_diff_arrival_rate = diff_arrival_rate;
                                }
                            }
                            max_diff_arrival_rate *= 2.0;
                            max_diff_index = num_arrival_rates - 1;
                            for(i=median_index;i<num_arrival_rates;i++){
                                cur_arrival_rate = arrival_rate_array[i];
                                diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
                                prev_arrival_rate = cur_arrival_rate;
                                if(diff_arrival_rate > max_diff_arrival_rate){
                                    max_diff_arrival_rate = diff_arrival_rate;
                                    max_diff_index = i - 1;
                                    break;
                                }
                            }
                            max_diff_arrival_rate = arrival_rate_array[max_diff_index];
                            tmp_fragment_arrival_rate =MIN(min_10_local_arrival_rate * 2.0,
                                                           median_local_arrival_rate * 1.25);
                            if(tmp_fragment_arrival_rate < recalibrated_fragment_arrival_rate){
                                recalibrated_fragment_arrival_rate = tmp_fragment_arrival_rate;
                            }
                            if(max_diff_arrival_rate < recalibrated_fragment_arrival_rate){
                                recalibrated_fragment_arrival_rate = max_diff_arrival_rate;
                            }
                        }
                        if(recalibrated_fragment_arrival_rate > _globalArrivalRate){
                            _globalArrivalRate = recalibrated_fragment_arrival_rate;
                            std::cerr <<
                                  "Used recalibrated global_fragment_arrival_rate="
                                << _globalArrivalRate << std::endl
                                << "Used recalibrated global_fragment_arrival_distance="
                                <<((_globalArrivalRate > 0.) ? 1./(_globalArrivalRate) : 0.)
                                << std::endl
                                << "Chunk arrival rates sorted at 1/100s"
                                << std::endl ;
                            for(i=0;i<100;i++) {
                               std::cerr<<arrival_rate_array[((num_arrival_rates * i) / 100)]
                                        << std::endl;
                            }
                            std::cerr << max_local_arrival_rate << std::endl;
                        }
                    }
                }
            }
		}else{
			// Compute actual GAR
			_globalArrivalRate = (float)total_random_frags_in_genome / (float)genome_size;
		}

		return(_globalArrivalRate);		

	}

	//////////////////////////////////////////////////////////////////////////////

	Unitig::Unitig(void){
		// Initialize values to unlikely values
		_localArrivalRate=-1;
		_covStat=FLT_MAX;
		_length=-1;
		_numFrags=-1;
		_numRandomFrags=-1;
        _avgRho = -1;
		dovetail_path_ptr = new DoveTailPath;
        _id=nextId++;
	}

	//////////////////////////////////////////////////////////////////////////////

	Unitig::~Unitig(void){
		if(dovetail_path_ptr!=NULL) delete dovetail_path_ptr;
	}

	//////////////////////////////////////////////////////////////////////////////

    void Unitig::addFrag( DoveTailNode node) {
        dovetail_path_ptr->push_back( node );
        _inUnitig[ node.ident ] = id();
    }

	//////////////////////////////////////////////////////////////////////////////

	float Unitig::getAvgRho(void){

        if(dovetail_path_ptr->size() == 1) 
            _avgRho = 1;
		if(_avgRho!=-1)
			return(_avgRho);
		
		
		// We will compute the average rho.
		//
		// Since rho is the length(unitig) - length(last fragment),
		//   and the length(last fragment) is ambiguous depending on which
		//   direction we are walking the unitig from.  We will take the average 
		//   of the rhos through both directions.

		DoveTailPath::const_iterator dtp_iter;
        if (dovetail_path_ptr == NULL || dovetail_path_ptr->empty()) {
            fprintf(stderr,"NULL dovetailpath\n");
            return 1;
        }

		// Get first fragment's length
		dtp_iter=dovetail_path_ptr->begin();
        int ident1 = dtp_iter->ident;
		long first_frag_len = BestOverlapGraph::fragLen(dtp_iter->ident);

		// Get last fragment's length
		dtp_iter=dovetail_path_ptr->end();
		dtp_iter--;
        int ident2 = dtp_iter->ident;
		long last_frag_len = BestOverlapGraph::fragLen(dtp_iter->ident);

		// Get average of first and last fragment lengths
		double avg_frag_len = (last_frag_len + first_frag_len)/2.0;
		
		// Compute average rho
		long unitig_length=getLength();
		_avgRho = unitig_length - avg_frag_len;
		
        if (_avgRho <= 0 ) {
            std::cerr << "Negative Rho fid1 " << ident1 << " fid2 " << ident2 <<
                   " ulen " << unitig_length << " len1 " << first_frag_len <<
                   " len2 " << last_frag_len << " avg " << avg_frag_len <<std::endl;
            _avgRho = 1;
        }
		return(_avgRho);
	}

	//////////////////////////////////////////////////////////////////////////////

	float Unitig::_globalArrivalRate = -1;

	void Unitig::setGlobalArrivalRate(float global_arrival_rate){
		_globalArrivalRate=global_arrival_rate;
	}
	void Unitig::setLocalArrivalRate(float local_arrival_rate){
        
        if ( local_arrival_rate < std::numeric_limits<float>::epsilon())
            _localArrivalRate = 0;
        else
		    _localArrivalRate = local_arrival_rate;
	}
	float Unitig::getLocalArrivalRate(){
		if (_localArrivalRate != -1 )
            return _localArrivalRate;
        setLocalArrivalRate((getNumFrags() - 1) / getAvgRho());
        return _localArrivalRate;
    }

	//////////////////////////////////////////////////////////////////////////////

	float Unitig::getCovStat(void){

		const float ln2=0.69314718055994530941723212145818;

		// Note that we are using numFrags in this calculation.
		//   If the fragments in the unitig are not randomly sampled
		//   from the genome, this calculation will be wrong.
		//   Ie. if fragments being assembled are from a locally
		//   sequenced batch, the region may look artificially repetitive.
		//   The value should really be "number of randomly sampled
		//   fragments in the unitig".

		if(_globalArrivalRate==-1){
			std::cerr << "You have not set the _globalArrivalRate variable." << std::endl;
		}

		if(_covStat == FLT_MAX){
			if(_globalArrivalRate > 0.0){
				_covStat = (getAvgRho() * _globalArrivalRate) - 
					(ln2 * (getNumFrags() -1));
			}else{
				_covStat = 0.0;
			}
		}

		return(_covStat);

	}

	//////////////////////////////////////////////////////////////////////////////

	long Unitig::getLength(void){
	// The length of the unitig is just the position of the fragment 
	//   with the large begin or end value 

		if(dovetail_path_ptr->size()==0){
			std::cerr << "This Unitig has an empty fragPositions." << std::endl;	
		}
        // used to use getLastBackboneNode(), but that's not required to be the 
        // highest base, we can have a read stick out further, but not contain the last
        int bgn = 0;
        int end = 0;
        DoveTailPath::reverse_iterator rIter = dovetail_path_ptr->rbegin();
        for(;rIter != dovetail_path_ptr->rend(); rIter++) {
            DoveTailNode node = *rIter;
            if (node.contained)
                continue;
            if (isReverse(node.position))
            {
                if (node.position.bgn > end) {
                    end = node.position.bgn;
                    bgn = node.position.end;
                }
                else if ( node.position.bgn < bgn ) {
                    break; // past all potenially longer reads, stop
                }
            } else { // forward read
                if (node.position.end > end) {
                    end = node.position.end;
                    bgn = node.position.bgn;
                }
                else if ( node.position.end < bgn ) {
                    break;
                }
            }
        }
		return end;
	}

	//////////////////////////////////////////////////////////////////////////////
	long Unitig::getSumFragLength(void){
	// The length of the unitig is just the position of the fragment 
	//   with the large begin or end value 

		long sum=0;

		if(dovetail_path_ptr->size()==0){
			std::cerr << "This Unitig has an empty dovetail." << std::endl;	
		}

        DoveTailPath::const_iterator fpm_itr;
		for(
		    fpm_itr=dovetail_path_ptr->begin();
		    fpm_itr!=dovetail_path_ptr->end();
		    fpm_itr++){

			SeqInterval intrvl=fpm_itr->position;
            sum += abs( intrvl.end - intrvl.bgn );
		}
		return(sum);
	}
	//////////////////////////////////////////////////////////////////////////////

	long Unitig::getNumFrags(void){

//		if(_numFrags!=-1){
//			return(_numFrags);
//		}
        if (!dovetail_path_ptr->empty()) {
            _numFrags = dovetail_path_ptr->size();
            return _numFrags;
        }

		long num_dovetailing_frags = dovetail_path_ptr->size();
		long num_contained_frags=0;

		// Total frags are the sum
		_numFrags = num_dovetailing_frags + num_contained_frags;
		return(_numFrags);

	}

	//////////////////////////////////////////////////////////////////////////////

	long Unitig::getNumRandomFrags(void){
	// This is a placeholder, random frags should not contain guides, or other
	//   fragments that are not randomly sampled across the whole genome.
		
		if(_numRandomFrags!=-1){
			return(_numRandomFrags);
		}	

		_numRandomFrags=getNumFrags();
		return(_numRandomFrags);
	}
	//////////////////////////////////////////////////////////////////////////////
    void Unitig::shiftCoordinates(int offset) {
        //simple version
        if (dovetail_path_ptr->empty()) {
        } else {
//        fprintf(stderr,"shift unitig by %d 1st frg %d\n",
//                offset,dovetail_path_ptr->front().ident);
        }
        DoveTailIter iter = dovetail_path_ptr->begin();
        for(; iter != dovetail_path_ptr->end(); iter++) {
            iter->position.bgn += offset;
            assert( iter->position.bgn >= 0);
            iter->position.end += offset;
            assert( iter->position.end >= 0);
        }
    }
	//////////////////////////////////////////////////////////////////////////////
    void Unitig::reverseComplement() {
        //simple version
        int length = this->getLength();
        DoveTailNode first = dovetail_path_ptr->front();
        DoveTailNode last = dovetail_path_ptr->back();
        fprintf(stderr,"Before reverse 1st %d %d %d last %d %d %d\n",
                first.ident,first.position.bgn,first.position.end,
                last.ident,last.position.bgn,last.position.end);
        DoveTailIter iter = dovetail_path_ptr->begin();
        for(; iter != dovetail_path_ptr->end(); iter++) {
            iter->position.bgn = length - iter->position.bgn;
            assert( iter->position.bgn >= 0);
            iter->position.end = length - iter->position.end;
            assert( iter->position.end >= 0);
        }
        reverse(dovetail_path_ptr->begin(),dovetail_path_ptr->end());
        first = dovetail_path_ptr->front();
        last = dovetail_path_ptr->back();
        fprintf(stderr,"After reverse 1st %d %d %d last %d %d %d\n",
                first.ident,first.position.bgn,first.position.end,
                last.ident,last.position.bgn,last.position.end);
    }
	//////////////////////////////////////////////////////////////////////////////
    void Unitig::reverseComplement(int offset, BestOverlapGraph *bog_ptr)
    {
        iuid notUsed;
        DoveTailNode last  = getLastBackboneNode(notUsed);
        int lastEnd = last.position.end > last.position.bgn ? 
                      last.position.end : last.position.bgn;
        int prevBeg = 0;
        DoveTailPath contains;
        DoveTailPath *revP = new DoveTailPath;
        DoveTailPath::reverse_iterator rend = dovetail_path_ptr->rend();
        DoveTailPath::reverse_iterator addIter = dovetail_path_ptr->rbegin();
        for(; addIter != rend; addIter++) {
            addIter->position.bgn = lastEnd - addIter->position.bgn + offset;
            assert( addIter->position.bgn >= 0);
            addIter->position.end = lastEnd - addIter->position.end + offset;
            assert( addIter->position.end >= 0);
            if (addIter->contained != 0) {
#ifdef NEW_UNITIGGER_INTERFACE
                int tmp = addIter->ahang;
                addIter->ahang = -addIter->bhang;
                addIter->bhang = -tmp;
#endif
                contains.push_back( *addIter );
            } else {
#ifdef NEW_UNITIGGER_INTERFACE
                int left = addIter->position.end < addIter->position.bgn ? 
                           addIter->position.end : addIter->position.bgn;
                assert( prevBeg <= left ); 
                DoveTailPath::reverse_iterator prev = addIter+1;
                while (prev != rend &&
                        prev->contained != 0) 
                {
                    prev = prev+1;
                }
                if (prev == rend) {

                    // After reverse, need to switch to other edge
                    BestEdgeOverlap* other;
                    iuid id = addIter->ident;
                    if (addIter->position.bgn > addIter->position.end) {
                        other = bog_ptr->getBestEdgeOverlap( id, FIVE_PRIME );
                    } else {
                        other = bog_ptr->getBestEdgeOverlap( id, THREE_PRIME );
                    }
                    addIter->ident2 = other->frag_b_id;
                    if (other->ahang < 0 && other->bhang < 0 ) {
                        addIter->ahang = -other->bhang;
                        addIter->bhang = -other->ahang;
                    } else {
                        addIter->ahang = other->ahang;
                        addIter->bhang = other->bhang;
                    }
                } else {
                    assert( prev->contained == 0 );
                    addIter->ident2 = prev->ident;
                    addIter->bhang = prev->ahang;
                    addIter->ahang = prev->bhang;
                }
#endif
                revP->push_back( *addIter );
                if (contains.size() > 0) {
                    revP->insert(revP->end(), contains.begin(), contains.end());
                    contains.clear();
                }
            }
        }
        delete dovetail_path_ptr;
        dovetail_path_ptr = revP;
    }

	//////////////////////////////////////////////////////////////////////////////
    //Recursively place all contains under this one into the FragmentPositionMap

    void Unitig::placeContains(const ContainerMap* cntnrp, BestContainmentMap *bestCtn,
                             const iuid container, const SeqInterval intvl, const int level)
    {
        if (cntnrp->size() == 0)
            return;

        ContainerMap::const_iterator ctmp_itr = cntnrp->find( container );
        if (ctmp_itr != cntnrp->end() ) {

                ContaineeList::const_iterator cntee_itr;
                for(cntee_itr  = ctmp_itr->second.begin();
                    cntee_itr != ctmp_itr->second.end(); cntee_itr++)
                {
                    iuid cntee = *cntee_itr;
                    BestContainment &best = (*bestCtn)[ cntee ];

                    if (best.isPlaced)
                        continue;

                    assert( best.container == container );
                    (*containPartialOrder)[ cntee ] = level;
                    if (level > maxContainDepth)
                        maxContainDepth = level;

                    int offset = best.a_hang;
                    int bhang  = best.b_hang;
                    DoveTailNode pos;
                    pos.type         = AS_READ; /* Isolated FT VR */
                    pos.ident        = cntee;
                    pos.contained    = container;
                    pos.delta_length = 0;
                    pos.delta        = NULL;

                    if(intvl.bgn < intvl.end) {
                        pos.position.bgn = intvl.bgn + offset;
                        pos.position.end = intvl.end + bhang;
#ifdef NEW_UNITIGGER_INTERFACE
                        pos.ident2       = container;
                        pos.ahang        = offset;
                        pos.bhang        = bhang;
#endif

                    } else if (intvl.bgn > intvl.end) {
                        pos.position.bgn = intvl.bgn - offset;
                        pos.position.end = intvl.end - bhang;
#ifdef NEW_UNITIGGER_INTERFACE
                        pos.ident2       = container;
                        // consensus seeems to need these reversed to get the
                        // correct alignment
                        pos.bhang        = -offset;
                        pos.ahang        = -bhang;
#endif
                    }else{
                        std::cerr << "Error, container size is zero." << std::endl;
                        assert(0);
                    }
                    // Swap ends if containee is not same strand as container
                    if(!best.sameOrientation){
                        int tmp          = pos.position.bgn;
                        pos.position.bgn = pos.position.end;
                        pos.position.end = tmp;
                    }

                    addFrag( pos );
                    best.isPlaced = true;
                    placeContains( cntnrp, bestCtn, cntee, pos.position, level+1);
                }
        }
    }

    void Unitig::recomputeFragmentPositions(ContainerMap *allcntnr_ptr,
                 BestContainmentMap *bestContain, BestOverlapGraph *bog_ptr)
    {
		long frag_ins_begin = 0;
		long frag_ins_end;
        iuid lastFrag = 0;
#ifdef NEW_UNITIGGER_INTERFACE
        iuid nextFrag = 0;
#endif
		//std::cerr << "Positioning dovetails." << std::endl;
		// place dovetails in a row
        if (dovetail_path_ptr == NULL)
            return;
//        computeFragmentPositions( bog_ptr );
        containPartialOrder = new std::map<iuid,int>;
        long numDoveTail = dovetail_path_ptr->size();
		for(long i=0; i < numDoveTail; i++) {
            DoveTailNode *dt_itr = &(*dovetail_path_ptr)[i];

            iuid fragId  = dt_itr->ident;
#ifdef NEW_UNITIGGER_INTERFACE
            if ( nextFrag != 0 )
                assert( nextFrag == fragId);
            nextFrag = dt_itr->ident2;
#endif
            lastFrag = fragId;

            placeContains( allcntnr_ptr, bestContain, fragId, dt_itr->position,1);
		}
        this->sort();
        delete containPartialOrder;
        // Compute assuming that containee is the same orientation as container
        //	if(cntnr_intvl.begin < cntnr_intvl.end)

        // Container is in forward direction
        //
        // |=========F0============|
        // 0          |=============CER=============>|
        //                    |=====CEE=====>|
        //
        // |---Cro----|
        //            |--Ceo--|
        // |-------Cep--------|
        //
        // Cro = container offset from beginning of unitig = cntnr_intvl.begin
        // Ceo = containee offset from 5' end of container = cntee->olap_offset
        // Cep = containee offset from beginning of unitig = cntee_intvl.begin
        // CEE fragment can be from either orientation since
        //   definition of olap_offset is based on 3' origin.

        // else if(cntnr_intvl.begin > cntnr_intvl.end)

        // Container is in reverse direction
        //
        // |=========F0============|
        // 0          |<============CER==============|
        //                    |<====CEE======|
        //
        // |---Cro----|
        //                                   |--Ceo--|
        // |-------Cep-----------------------|
        //
        // Cro = container offset from beginning of unitig = cntnr_intvl.end
        // Ceo = containee offset from 5' end of container = cntee->olap_offset
        // Cep = containee offset from beginning of unitig = cntee_intvl.end
        // CEE fragment can be from either orientation since
        //   definition of olap_offset is based on 3' origin.

	}
	//////////////////////////////////////////////////////////////////////////////
    int lastBPCoord   = 0;
    int lastBPFragNum = 0;
    UnitigBreakPoint UnitigGraph::selectSmall(const Unitig *tig,
           const FragmentEnds &smalls, const UnitigBreakPoint& big)
    {
        double difference = 0.0;
        UnitigBreakPoint selection;
        int rFrgs =  big.fragNumber; 
        bool bRev = isReverse( big.position);
        iuid bid =  big.fragEnd.id;
        int right,bContain;
        bContain = 0;
        if ( big.fragEnd.end == FIVE_PRIME) right = big.position.bgn;
        else                                right = big.position.end;

        if ( cntnrmap_ptr->find( bid ) != cntnrmap_ptr->end() )
            bContain = (*cntnrmap_ptr)[ bid ].size();

        if ( (bRev && big.fragEnd == THREE_PRIME) ||
            ( !bRev && big.fragEnd == FIVE_PRIME) ) 
        {
            rFrgs -= 1 + bContain;
        }
        FragmentEnds::const_iterator sIter = smalls.begin();
        for(; sIter != smalls.end(); sIter++) {
            UnitigBreakPoint small = *sIter;
            if (small.fragEnd.id == big.fragEnd.id)
                continue; 

            bool rev  = isReverse(small.position);
            int sContain = 0;
            int lFrgs = small.fragNumber;
            iuid sid = small.fragEnd.id;

            // use middle of frag instead of end, to get some overlap
            double bp = (small.position.end + small.position.bgn) / 2.0;

//            if (small.fragEnd.end == FIVE_PRIME) bp = small.position.bgn;
//            else                                 bp = small.position.end;

            if ( cntnrmap_ptr->find( sid ) != cntnrmap_ptr->end() )
                sContain = (*cntnrmap_ptr)[ sid ].size();

            if ( (rev && small.fragEnd == THREE_PRIME) ||
                (!rev && small.fragEnd == FIVE_PRIME ) ) 
            // left side of the frag in the unitig, don't count it
            {
                lFrgs -= 1 + sContain;
            }
            if (rFrgs - lFrgs == 1)
                continue;
            double lRate = (lFrgs - lastBPFragNum) / (bp - lastBPCoord); 
            double rRate = (rFrgs - lFrgs) / (right - bp);
            fprintf(stderr,
              "Break frg %7d b %4d l %4d pos b %5d e %5.0f lRate %.4f\n",
              sid, lastBPFragNum, lFrgs, lastBPCoord, bp, lRate );
            double ratio = lRate > rRate ? lRate / rRate : rRate / lRate;
            fprintf(stderr,
              "     diff %4d r %4d pos %5d rRate %.4f ratio %.2f to frag %7d\n",
               rFrgs - lFrgs, rFrgs, right, rRate, ratio, bid ); 

            double diff = fabs( lRate - rRate);
            if (ratio > 1.8 && diff > difference) {
                fprintf(stderr,"Select frg %d for break on arrival rate diff %.4f\n",
                        sid, diff);
                difference = diff;
                selection = small;
            }
        }
        lastBPCoord   = right;
        lastBPFragNum = rFrgs;
        return selection;
    }

    void UnitigGraph::filterBreakPoints(Unitig *tig, FragmentEnds &breaks)
    {
        UnitigBreakPoint fakeEnd = breaks.back();
        // Check for sentinal at end, not needed by MateChecker
        if (fakeEnd.inSize == std::numeric_limits<int>::max())
            breaks.pop_back();
        FragmentEnds smallBPs, newBPs;
        bool hadBig = false;
        FragmentEnds::iterator iter = breaks.begin();
        for(; iter != breaks.end(); iter++)
        {
            UnitigBreakPoint nextBP = *iter;
            //if (nextBP.inSize > 2500) 
            if (nextBP.inFrags > 1 && nextBP.inSize > 500) {
                hadBig = true;
                // big one, compare against smalls
                if (newBPs.empty() || !(nextBP.fragEnd == newBPs.back().fragEnd)) {
                    if (!smallBPs.empty()) {
                        UnitigBreakPoint small = selectSmall( tig, smallBPs, nextBP);
                        if ( small.fragNumber > 0 )
                            newBPs.push_back( small );
                        smallBPs.clear();
                    } else {
                        if ( nextBP.fragEnd.end == FIVE_PRIME)
                            lastBPCoord = nextBP.position.bgn;
                        else
                            lastBPCoord = nextBP.position.end;

                        iuid bid = nextBP.fragEnd.id;
                        int bContain = 0;
                        if ( cntnrmap_ptr != NULL &&
                             cntnrmap_ptr->find( bid ) != cntnrmap_ptr->end() )
                            bContain = (*cntnrmap_ptr)[ bid ].size();

                        lastBPFragNum = nextBP.fragNumber;
                        bool bRev = isReverse( nextBP.position );
                        if ( (bRev && nextBP.fragEnd == THREE_PRIME) ||
                           ( !bRev && nextBP.fragEnd == FIVE_PRIME) ) 
                        {
                            lastBPFragNum -= 1 + bContain;
                        }
                    }
                    newBPs.push_back( nextBP );
                }
            } else 
                if ( newBPs.empty() ||
                    !(nextBP.fragEnd == newBPs.back().fragEnd) &&
                    !(nextBP.fragEnd == smallBPs.back().fragEnd) ) 
                    smallBPs.push_back( nextBP );
        }
        if (!hadBig)
            assert( newBPs.size() == 0);
        if (hadBig && !smallBPs.empty()) {
            UnitigBreakPoint small = selectSmall( tig, smallBPs, fakeEnd);
            if ( small.fragNumber > 0 )
                newBPs.push_back( small );
            smallBPs.clear();
        }
        breaks = newBPs;
    }
    // Doesn't handle contained frags yet
    UnitigVector* UnitigGraph::breakUnitigAt(Unitig *tig, FragmentEnds &breaks)
    {
        if (breaks.empty())
            return NULL;
        fprintf(stderr,"Before break size is %d, unitig %d\n",breaks.size(),tig->dovetail_path_ptr->size());
        lastBPCoord = 0;
        lastBPFragNum = 0;
        // remove small break points
        filterBreakPoints( tig, breaks );
        fprintf(stderr,"After filter break size is %d\n",breaks.size());
        UnitigVector* splits = new UnitigVector(); // holds the new split unitigs
        // if we filtered all the breaks out return an empty list
        if (breaks.empty()) 
            return splits;
        Unitig* newTig = new Unitig(); // the first of the split unitigs
        splits->push_back( newTig );
        DoveTailIter dtIter = tig->dovetail_path_ptr->begin();
        UnitigBreakPoint breakPoint = breaks.front();
        breaks.pop_front();
        bool bothEnds = false;
        UnitigBreakPoint nextBP;
        if (!breaks.empty())
            nextBP = breaks.front();
        int frgCnt = 0;
        int offset = 0;
        for( ; dtIter != tig->dovetail_path_ptr->end(); dtIter++) {
            frgCnt++;
            DoveTailNode frg = *dtIter;
            // reduce multiple breaks at the same fragment end down to one
            while ( !breaks.empty() && nextBP.fragEnd.id == breakPoint.fragEnd.id ) {
                if (nextBP.fragEnd.end != breakPoint.fragEnd.end)
                    bothEnds = true;
                breaks.pop_front();
                nextBP = breaks.front();
            }
            bool reverse = isReverse(frg.position);
            if (breakPoint.fragEnd.id == frg.ident) {
                if (frg.contained) {
                    fprintf(stderr,"Breaking untig on contained frg %d\n",frg.ident);
                    frg.contained = NULL_FRAG_ID;
                }
                // At a fragment to break on
                if (bothEnds) {
                    // create singleton when breaking at both ends of a fragment
                    fprintf(stderr,"  Breaking tig %d at both ends of %d num %d\n",
                            tig->id(), breakPoint.fragEnd.id, breakPoint.fragNumber);
                    newTig->shiftCoordinates( offset ); // set coords for last piece
                    newTig = new Unitig();              // create the singlenton
                    splits->push_back( newTig );        // add to list of new untigs
                    frg.position.bgn = 0;               // set begin and end
                    frg.position.end = BestOverlapGraph::fragLen( frg.ident );
                    newTig->addFrag( frg );             // add frag to unitig
                    newTig = new Unitig();              // create next new unitig
                    bothEnds = false;
                }
                else if (breakPoint.fragEnd.end ==  FIVE_PRIME && !reverse ||
                         breakPoint.fragEnd.end == THREE_PRIME && reverse)
                {
                    // break at left end of frg, frg starts new tig
                    fprintf(stderr,"  Break tig %d before %d num %d\n",
                            tig->id(), breakPoint.fragEnd.id, breakPoint.fragNumber);
                    newTig->shiftCoordinates( offset ); // set coords for last piece
                    newTig = new Unitig();              // create next new unitig
                    // record offset of start of new unitig
                    offset = reverse ? -frg.position.end : -frg.position.bgn;
                    fprintf(stderr,"Offset is %d for frg %d %d,%d ",
                            offset,frg.ident,frg.position.bgn,frg.position.end);
                    newTig->addFrag( frg );             // add frag as 1st in unitig
                }
                else if (breakPoint.fragEnd.end ==  FIVE_PRIME && reverse ||
                         breakPoint.fragEnd.end == THREE_PRIME && !reverse )
                {
                    // break at right end of frg, frg goes in existing tig
                    fprintf(stderr,"  Break tig %d after %d num %d\n",
                            tig->id(), breakPoint.fragEnd.id, breakPoint.fragNumber);
                    newTig->addFrag( frg );             // add frag as last
                    newTig->shiftCoordinates( offset );
                    newTig = new Unitig();              // create new empty unitig
                } else {
                    // logically impossible!
                    assert(0);
                }
                splits->push_back( newTig );           // add new tig to list of splits
                if (breaks.empty()) {
                    // Done breaking, continue adding remaining frags to newTig
                    breakPoint.fragEnd.id = 0;
                    continue;
                }
                else {
                    breakPoint = breaks.front();
                    breaks.pop_front();
                    if (!breaks.empty())
                        nextBP = breaks.front();
                }
            } else {
                // not a unitig break point fragment, add to current new tig
                if (newTig->dovetail_path_ptr->empty()) {
                    // if it's empty, this is it's 1st frag, store the offset
                    offset = reverse ? -frg.position.end : -frg.position.bgn;
                    fprintf(stderr,"Offset is %d for frg %d %d,%d ",
                            offset,frg.ident,frg.position.bgn,frg.position.end);
                }
                if (frg.contained) {
                    if (Unitig::fragIn(frg.contained) != newTig->id())
                        frg.contained = NULL_FRAG_ID;
                }
                newTig->addFrag( frg );
            }
        }
        newTig->shiftCoordinates( offset );
        // should be zero
        fprintf(stderr,"After break size is %d\n",breaks.size());
        return splits;
    }

    //////////////////////////////////////////////////////////////////////////////

    void Unitig::freeIUM_Mesg(IntUnitigMesg *ium_ptr){
        delete ium_ptr;
    }

    //////////////////////////////////////////////////////////////////////////////

    int IntMultiPosCmp(const void *a, const void *b){
        IntMultiPos *impa=(IntMultiPos*)a;
        IntMultiPos *impb=(IntMultiPos*)b;
        long aleft,aright,bleft,bright;
        if (impa->position.bgn < impa->position.end) {
            aleft  = impa->position.bgn;
            aright = impa->position.end;
        } else {
            aright = impa->position.bgn;
            aleft  = impa->position.end;
        }
        if (impb->position.bgn < impb->position.end) {
            bleft  = impb->position.bgn;
            bright = impb->position.end;
        } else {
            bright = impb->position.bgn;
            bleft  = impb->position.end;
        }
        if(aleft!=bleft)
        {
            return(aleft - bleft);
        }
        else if (aright != bright)
        {
            if (impa->contained==0 && impb->contained!=0)
                return(-1);
            else if (impa->contained!=0 && impb->contained==0)
                return(1);
            if (impa->contained!=0 || impb->contained!=0)
                return(bright - aright);
            else 
                return(aright - bright);
        }
        else {
            if(impa->contained == impb->ident)
                return(1);
            if(impb->contained == impa->ident)
                return(-1);
            if(impa->contained!=0 && impb->contained!=0)
                return((*containPartialOrder)[impa->ident] - (*containPartialOrder)[impb->ident]);
            if(impa->contained!=0)
                return(1);
            if(impb->contained!=0)
                return(-1);
            return(0);
        }
    }
    //////////////////////////////////////////////////////////////////////////////
    void Unitig::sort() {
        qsort( &(dovetail_path_ptr->front()), getNumFrags(),
                sizeof(IntMultiPos), &IntMultiPosCmp );
    }
    //////////////////////////////////////////////////////////////////////////////
    IntUnitigMesg *Unitig::getIUM_Mesg(){

        IntUnitigMesg *ium_mesg_ptr=new IntUnitigMesg;
        IntMultiPos *imp_msg_arr   = &(dovetail_path_ptr->front());

        //void qsort(void *base, size_t nmemb, size_t size,
        //  int(*compar)(const void *, const void *));
        // Sort by 5' end of fragment position
        //qsort(imp_msg_arr, getNumFrags(), sizeof(IntMultiPos), &IntMultiPosCmp);

        // Populate the IUM message with unitig info
#ifdef AS_ENABLE_SOURCE
        /*char*/		ium_mesg_ptr->source=AS_BOG_UNITIG_GRAPH_CC_CM_ID;
#endif
        /*float  */		ium_mesg_ptr->coverage_stat=getCovStat();
        /*UnitigStatus*/	ium_mesg_ptr->status=AS_UNASSIGNED;
        /*UnitigFUR*/   ium_mesg_ptr->unique_rept=AS_FORCED_NONE;
        /*CDS_COORD_t*/		ium_mesg_ptr->length=getLength();
        /*char* */		ium_mesg_ptr->consensus="";
        /*char* */		ium_mesg_ptr->quality="";
        /*int32*/		ium_mesg_ptr->forced=0;
        /*int32*/		ium_mesg_ptr->num_frags=getNumFrags();
        /*IntMultiPos* */	ium_mesg_ptr->f_list=imp_msg_arr;

        return(ium_mesg_ptr);
    }

    //////////////////////////////////////////////////////////////////////////////

	void UnitigGraph::writeIUMtoFile(char *filename){
		
		// Open up the output file
		FILE *fptr=fopen(filename, "w");
		if(fptr==NULL){
			std::cerr << "Could not open " << filename << 
				"for writing IUM (IntUnitigMesg) output." << std::endl;
			
			exit(-1);
		}

		// Step through all the unitigs
        int iumId = 0;
		UnitigVector::iterator utg_itr;
		for(utg_itr=unitigs->begin(); utg_itr!=unitigs->end(); utg_itr++){
            if( *utg_itr == NULL)
                continue;

            if( (*utg_itr)->getNumFrags() == 0 )
                continue;

			IntUnitigMesg *ium_mesg_ptr;
			ium_mesg_ptr=(*utg_itr)->getIUM_Mesg();
            fprintf(stderr,"Internal tig %ld is %ld in IUM msg.\n",
                    (*utg_itr)->id(),iumId );
			ium_mesg_ptr->iaccession=iumId++; //IntChunk_ID

			GenericMesg mesg;
			mesg.m=ium_mesg_ptr;
			mesg.t=MESG_IUM;

			WriteProtoMesg_AS(fptr, &mesg);
			(*utg_itr)->freeIUM_Mesg(ium_mesg_ptr);
		}

		fclose(fptr);
	}

	//////////////////////////////////////////////////////////////////////////////


} //AS_BOG namespace
