
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

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include <set>
#include <limits>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cerrno>

extern "C" {
#include "AS_global.h"
}
#undef max
namespace AS_BOG{

    //////////////////////////////////////////////////////////////////////////////
    
    int BogOptions::badMateBreakThreshold    = -7;
    bool BogOptions::unitigIntersectBreaking = false;
    bool BogOptions::ejectUnhappyContained   = false;
    bool BogOptions::useGkpStoreLibStats     = false;

    //////////////////////////////////////////////////////////////////////////////
    UnitigGraph::UnitigGraph(BestOverlapGraph *in_bog_ptr) : bog_ptr(in_bog_ptr) {
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

    void UnitigGraph::build(ChunkGraph *cg_ptr) {

        iuid num_frags=cg_ptr->getNumFragments();
        Unitig::setNextId( 1 );

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
                best_cntr->find(frag_idx) == best_cntr->end() ) { 
                cg_ptr->getChunking(frag_idx, 
                                    fp_dst_frag_id, 
                                    tp_dst_frag_id);
                
                // Allocated a new unitig node
                Unitig *utg=new Unitig;

                // create it going off the 5' end
                populateUnitig( utg, frag_idx, FIVE_PRIME, cg_ptr, 0);

                // now check if we can also go off 3' end
                fragment_end_type whichEnd = THREE_PRIME;
                BestEdgeOverlap *tpBest = bog_ptr->getBestEdgeOverlap( frag_idx, whichEnd);
                iuid tpId = tpBest->frag_b_id;

                // if the other end is already in another unitig, don't try to extend
                if (Unitig::fragIn(tpId)) {
                    // question, should we include the circle in the intersection list?
                    if (tpId != 0 && utg->id() != Unitig::fragIn(tpId)) {
                        // save this outgoing intersection point for future use
                        (*unitigIntersect)[tpId].push_back( frag_idx );
                        fprintf(stderr,"Unitig %5d 1st frag %7d -> Unitig %5d frag %7d\n",
                                utg->id(), frag_idx, Unitig::fragIn(tpId), tpId );
                    }
                } else {
                    // if other end is also 3' we need to walk of it's 5'
                    if (tpBest->bend == whichEnd)
                        whichEnd = FIVE_PRIME;

                    int len = utg->getLength();
                    int flen = BestOverlapGraph::fragLen(frag_idx);
                    int ahng = tpBest->ahang;
                    int offset = len - flen + ahng;
                    fprintf(stderr, "Join unitig %d len %d at %d len %d and %d ahang %d offset %d\n",
                            utg->id(), len, frag_idx, flen, tpBest->frag_b_id, ahng, offset );

                    utg->reverseComplement();
                    populateUnitig( utg, tpBest->frag_b_id, whichEnd, cg_ptr, offset);
                }

                // Store unitig in unitig graph
                unitigs->push_back(utg);

            }else{
                // We are either in the middle of a dovetail sequence,
                // it is a singleton, or we are going from 3' to 5'
                // across a dovetail.

            }
        }

        // Pick up frags missed above, possibly circular unitigs
        for(frag_idx=1; frag_idx<=num_frags; frag_idx++){
            if (Unitig::fragIn( frag_idx ) == 0) {
                cg_ptr->getChunking( frag_idx, fp_dst_frag_id, tp_dst_frag_id);

                if ( fp_dst_frag_id != NULL_FRAG_ID &&
                     tp_dst_frag_id != NULL_FRAG_ID ) {
                    Unitig *utg=new Unitig;
                    // need to make populateUnitig break circle
                    populateUnitig( utg, frag_idx, FIVE_PRIME, cg_ptr, 0 );

                    fprintf(stderr,"Circular unitig %d 1st frag %d\n", utg->id(), frag_idx);

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
                                              BestEdgeOverlap *&fivePrime, BestEdgeOverlap *&threePrime ) {
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
                                   std::map<iuid,iuid> *visited_map) {
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
               visited_map->find(joiner) != visited_map->end()) {
#warning DANGEROUS assume unitig is at id-1 in vector
            Unitig* tigToAdd = unitigs->at( tigIdToAdd - 1);
            assert(tigToAdd->id() == tigIdToAdd);
            // Code assumes joining only from end of unitig
            if ( beforeLast != 0 ) // ok to join singleton though
                assert( beginId != tig->dovetail_path_ptr->front().ident);
            BestEdgeOverlap* addFivePrime = bog_ptr->getBestEdgeOverlap(
                                                                        joiner, FIVE_PRIME );
            BestEdgeOverlap* addThreePrime = bog_ptr->getBestEdgeOverlap(
                                                                         joiner, THREE_PRIME );

            bool bDisagree = false;
            if ( bestEdge->bend == FIVE_PRIME  && beginId != addFivePrime->frag_b_id
                 || bestEdge->bend == THREE_PRIME && beginId != addThreePrime->frag_b_id) { 
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

            int offset  = 1 + tigEnd - BestOverlapGraph::olapLength(beginId, joiner, bestEdge->ahang, bestEdge->bhang);

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
#warning DANGEROUS assume unitig is at id-1 in vector
            assert(unitigs->at(tigIdToAdd-1)->id() == tigIdToAdd);
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
        DoveTailNode lastNonContain;

        memset(&lastNonContain, 0, sizeof(DoveTailNode));

        for(DoveTailPath::reverse_iterator rIter = dovetail_path_ptr->rbegin(); rIter != dovetail_path_ptr->rend(); rIter++) {
            if (rIter->contained == 0) {
                if (lastNonContain.ident == 0) {
                    lastNonContain = *rIter;
                } else {
                    prevId = rIter->ident;
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
    void UnitigGraph::breakUnitigs() {

        //  Stop when we've seen all current unitigs.  Replace tiMax
        //  in the for loop below with unitigs->size() to recursively
        //  split unitigs.
        int  tiMax = unitigs->size();

        for (int  ti=0; ti<tiMax; ti++) {
            Unitig        *tig = (*unitigs)[ti];
            FragmentEnds   breaks;
            int            fragCount = 1;
            DoveTailNode   lastBackbone;

            tig->sort();

            for (DoveTailPath::const_iterator dt_itr = tig->dovetail_path_ptr->begin();
                 dt_itr != tig->dovetail_path_ptr->end();
                 dt_itr++,fragCount++) {

                DoveTailNode node = *dt_itr;
                iuid dtFrag = node.ident;

                // Include contained frags in count for coverage calculation
                if ( cntnrmap_ptr->find( dtFrag ) != cntnrmap_ptr->end())
                    fragCount += (*cntnrmap_ptr)[ dtFrag ].size();

                // Detect first or last frg in unitig
                if (dt_itr   == tig->dovetail_path_ptr->begin() ||
                    dt_itr+1 == tig->dovetail_path_ptr->end() ) {
                    // default to 3' end of untig for the last frag
                    fragment_end_type dtEnd = THREE_PRIME;
                    char              dtStr = '3';

                    // want the backbone, not a contained
                    if (node.contained)
                        node = lastBackbone;

                    if ( isReverse( node.position ) )
                        dtEnd = FIVE_PRIME;

                    // At the begining we want the 5' end of unitig, so switch
                    if (dt_itr == tig->dovetail_path_ptr->begin()) {
                        dtEnd = (dtEnd == FIVE_PRIME) ? THREE_PRIME : FIVE_PRIME;
                        dtStr = '5';
                    }

                    BestEdgeOverlap *bestEdge = bog_ptr->getBestEdgeOverlap(node.ident,dtEnd);

                    fprintf(stderr,"Unitig %5d %c' frag %7d points to Unitig %5d frag %7d\n",
                            tig->id(), dtStr, node.ident, Unitig::fragIn(bestEdge->frag_b_id), bestEdge->frag_b_id);
                }

                // store previous backbone frag for next loop iteration
                if (!node.contained)
                    lastBackbone = node;

                FragmentEdgeList::const_iterator edge_itr = unitigIntersect->find( dtFrag );

                if (edge_itr != unitigIntersect->end() ) {

                    // We have a set of best edges incoming from other unitigs
                    for (FragmentList::const_iterator fragItr = edge_itr->second.begin();
                         fragItr != edge_itr->second.end();
                         fragItr++) {
                        iuid inFrag = *fragItr;

                        // check if it's incoming frag's 5' best edge
                        fragment_end_type bestEnd = FIVE_PRIME;
                        BestEdgeOverlap  *bestEdge = bog_ptr->getBestEdgeOverlap( inFrag, bestEnd );

                        if (bestEdge->frag_b_id != dtFrag) {
                            // not 5' best edge, so 3'
                            bestEnd = THREE_PRIME;
                            bestEdge = bog_ptr->getBestEdgeOverlap( inFrag, bestEnd );
                        }

                        // dtFrag must be either 3' or 5' best edge of inFrag
                        assert( bestEdge->frag_b_id == dtFrag );

                        int pos = (bestEdge->bend == FIVE_PRIME) ? node.position.bgn : node.position.end;

                        // get incoming unitig, unitig id's start at 1, storage at 0
                        Unitig *inTig = (*unitigs)[Unitig::fragIn(inFrag)-1];

                        if (inTig == NULL) {
                            fprintf(stderr, "  NullBreak tig %5d at frag %7d\n", tig->id(), dtFrag );
                        } else {
                            UnitigBreakPoint breakPoint(dtFrag, bestEdge->bend);
                            breakPoint.position    = node.position;
                            breakPoint.fragNumber  = fragCount ;
                            breakPoint.inSize      = inTig->getLength();
                            breakPoint.inFrags     = inTig->getNumFrags();
                            breaks.push_back( breakPoint );
                        }

                        fprintf(stderr, "  Utig %4d frgs %5d len %4d frag %6d end %c' into Utig %4d frag %6d end %c' pos %6d\n",
                                Unitig::fragIn(inFrag), inTig->getNumFrags(),
                                inTig->getLength(), inFrag, (bestEnd == FIVE_PRIME) ? '5' : '3',
                                tig->id(), dtFrag, (bestEdge->bend == FIVE_PRIME) ? '5' : '3', pos );
                    }
                }
            }

            if (breaks.empty() == false) {

                //  create a final fake bp for the last frag so we
                //  have a reference point to the end of the tig for
                //  filtering the breakpoints.  fakeEnd, for searching.
                //
                UnitigBreakPoint breakPoint(lastBackbone.ident, (isReverse(lastBackbone.position)) ? FIVE_PRIME : THREE_PRIME);
                breakPoint.position   = lastBackbone.position;
                breakPoint.fragNumber = fragCount;
                breakPoint.inSize     = std::numeric_limits<int>::max();
                breaks.push_back( breakPoint );

                UnitigVector* newUs = breakUnitigAt( tig, breaks );

                if (newUs != NULL) {
                    delete tig;
                    (*unitigs)[ti] = NULL;
                    unitigs->insert(unitigs->end(), newUs->begin(), newUs->end());
                }

                delete newUs;
            }
        }
    }




    //////////////////////////////////////////////////////////////////////////////
    void UnitigGraph::printUnitigBreaks() {

        for(UnitigVector::const_iterator tigIter = unitigs->begin(); tigIter != unitigs->end(); tigIter++) {
            Unitig*      tig = *tigIter;
            DoveTailNode first = tig->dovetail_path_ptr->front();
            iuid         prev = 0;
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

            if (firstBest->frag_b_id == lastBest->frag_b_id)
                fprintf(stderr, "Self Bubble %d to %d\n",id1,id2);  
        }
    }
    //////////////////////////////////////////////////////////////////////////////

    void UnitigGraph::_build_container_map(BestContainmentMap *best_cntr) {
        cntnrmap_ptr = new ContainerMap;
        int containees=0;

        BestContainmentMap::const_iterator bstcnmap_itr;
        for( bstcnmap_itr  = best_cntr->begin(); 
             bstcnmap_itr != best_cntr->end(); bstcnmap_itr++) {
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
                
                for( iter=unitigs->begin(); iter!=unitigs->end(); iter++) {
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

    Unitig::Unitig(bool report){
        _localArrivalRate = -1;
        _covStat          = FLT_MAX;
        _length           = -1;
        _numFrags         = -1;
        _numRandomFrags   = -1;
        _avgRho           = -1;
        dovetail_path_ptr = new DoveTailPath;
        _id               = nextId++;

        if (report)
            fprintf(stderr, "Creating Unitig %d\n", _id);
    }
    //////////////////////////////////////////////////////////////////////////////

    Unitig::~Unitig(void){
        delete dovetail_path_ptr;
    }

    //////////////////////////////////////////////////////////////////////////////

    void Unitig::addFrag( DoveTailNode node, int offset, bool report) {

        node.position.bgn += offset;
        node.position.end += offset;

        assert(node.position.bgn >= 0);
        assert(node.position.end >= 0);

        // keep track of the unitig a frag is in
        _inUnitig[ node.ident ] = id();

        // keep track of max position in unitig
        int frgEnd = MAX( node.position.bgn, node.position.end);
        if ( frgEnd > _length)
            _length = frgEnd;

        //if ((node.contained) && (_inUnitig[node.contained] != id()))
        //    node.contained = NULL_FRAG_ID;

        dovetail_path_ptr->push_back(node);

        if (report)
            if (node.contained)
                fprintf(stderr, "Added frag %d to unitig %d at %d,%d (contained in %d)\n", node.ident, id(), node.position.bgn, node.position.end, node.contained);
            else
                fprintf(stderr, "Added frag %d to unitig %d at %d,%d\n", node.ident, id(), node.position.bgn, node.position.end);
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
            //fprintf(stderr,"NULL dovetailpath\n");
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
        return _length;
    }

    //////////////////////////////////////////////////////////////////////////////
    long Unitig::getSumFragLength(void){
        long sum=0;

        DoveTailPath::const_iterator fpm_itr;
        for(fpm_itr=dovetail_path_ptr->begin();
            fpm_itr!=dovetail_path_ptr->end();
            fpm_itr++){

            SeqInterval intrvl=fpm_itr->position;
            sum += abs( intrvl.end - intrvl.bgn );
        }
        return(sum);
    }
    //////////////////////////////////////////////////////////////////////////////

    long Unitig::getNumFrags(void){

        if (!dovetail_path_ptr->empty()) {
            _numFrags = dovetail_path_ptr->size();
        } else {
          long num_dovetailing_frags = dovetail_path_ptr->size();
          long num_contained_frags   = 0;

          // Total frags are the sum
          _numFrags = num_dovetailing_frags + num_contained_frags;
        }
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
#if 0
        fprintf(stderr,"Before reverse 1st %d %d %d last %d %d %d\n",
                first.ident,first.position.bgn,first.position.end,
                last.ident,last.position.bgn,last.position.end);
#endif
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
#if 0
        fprintf(stderr,"After reverse 1st %d %d %d last %d %d %d\n",
                first.ident,first.position.bgn,first.position.end,
                last.ident,last.position.bgn,last.position.end);
#endif
    }
    //////////////////////////////////////////////////////////////////////////////
    void Unitig::reverseComplement(int offset, BestOverlapGraph *bog_ptr) {
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
                       prev->contained != 0) {
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
    //
    // Compute assuming that containee is the same orientation as container
    //	if(cntnr_intvl.begin < cntnr_intvl.end)
    //
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
    //
    // else if(cntnr_intvl.begin > cntnr_intvl.end)
    //
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

    void Unitig::placeContains(const ContainerMap* cntnrp,
                               BestContainmentMap *bestCtn,
                               const iuid containerId,
                               const SeqInterval containerPos,
                               const int level) {
        if (cntnrp->size() == 0)
            return;

        ContainerMap::const_iterator ctmp_itr = cntnrp->find( containerId );

        if (ctmp_itr == cntnrp->end() )
            return;

        for(ContaineeList::const_iterator cntee_itr  = ctmp_itr->second.begin();
            cntee_itr != ctmp_itr->second.end();
            cntee_itr++) {

            iuid cntee = *cntee_itr;
            BestContainment &best = (*bestCtn)[ cntee ];

            if (best.isPlaced)
                continue;

            assert( best.container == containerId );

            (*containPartialOrder)[ cntee ] = level;
            if (level > maxContainDepth)
                maxContainDepth = level;

            if (cntee == 132732) {
                fprintf(stderr, "132732\n");
            }

            DoveTailNode pos;

            pos.type         = AS_READ; /* Isolated FT VR */
            pos.ident        = cntee;
            pos.contained    = containerId;
            pos.delta_length = 0;
            pos.delta        = NULL;

            if(containerPos.bgn < containerPos.end) {
                //  Container is forward
                pos.position.bgn = containerPos.bgn + best.a_hang;  //  BPW says "looks ok"
                pos.position.end = containerPos.end + best.b_hang;
#ifdef NEW_UNITIGGER_INTERFACE
                pos.ident2       = container;
                pos.ahang        = best.a_hang;
                pos.bhang        = best.b_hang;
#endif

            } else if (containerPos.bgn > containerPos.end) {
                //  Container is reverse
                pos.position.bgn = containerPos.bgn - best.a_hang;  //  BPW says "suspicious"
                pos.position.end = containerPos.end - best.b_hang;
#ifdef NEW_UNITIGGER_INTERFACE
                pos.ident2       = containerId;
                pos.ahang        = - best.b_hang;   //  consensus seems to want these reversed
                pos.bhang        = - best.a_hang;
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

            addFrag( pos, 0, false );
            best.isPlaced = true;
            placeContains( cntnrp, bestCtn, cntee, pos.position, level+1);
        }
    }

    void Unitig::recomputeFragmentPositions(ContainerMap *allcntnr_ptr,
                                            BestContainmentMap *bestContain,
                                            BestOverlapGraph *bog_ptr) {
        long frag_ins_begin = 0;
        long frag_ins_end;
        iuid lastFrag = 0;
#ifdef NEW_UNITIGGER_INTERFACE
        iuid nextFrag = 0;
#endif

        //fprintf(stderr, "==> PLACING CONTAINED FRAGMENTS\n");

        // place dovetails in a row
        if (dovetail_path_ptr == NULL)
            return;

        containPartialOrder = new std::map<iuid,int>;
        long numDoveTail = dovetail_path_ptr->size();
        for(long i=0; i < numDoveTail; i++) {
            DoveTailNode *dt_itr = &(*dovetail_path_ptr)[i];

#ifdef NEW_UNITIGGER_INTERFACE
            if ( nextFrag != 0 )
                assert( nextFrag == dt_itr->ident);
            nextFrag = dt_itr->ident2;
#endif
            lastFrag = dt_itr->ident;

            placeContains(allcntnr_ptr, bestContain, dt_itr->ident, dt_itr->position, 1);
        }
        this->sort();
        delete containPartialOrder;
        containPartialOrder = NULL;
    }

    //////////////////////////////////////////////////////////////////////////////




    UnitigBreakPoint UnitigGraph::selectSmall(const Unitig *tig,
                                              const FragmentEnds &smalls,
                                              const UnitigBreakPoint& big,
                                              int   &lastBPCoord,
                                              int   &lastBPFragNum) {
        UnitigBreakPoint selection;

        double difference = 0.0;
        int  rFrgs   = big.fragNumber; 
        bool bRev    = isReverse(big.position);
        int bContain = 0;

        int right    = (big.fragEnd.end == FIVE_PRIME) ? big.position.bgn : big.position.end;

        if (cntnrmap_ptr->find(big.fragEnd.id) != cntnrmap_ptr->end())
            bContain = (*cntnrmap_ptr)[big.fragEnd.id].size();

        if (((bRev == true)  && (big.fragEnd == THREE_PRIME)) ||
            ((bRev == false) && (big.fragEnd == FIVE_PRIME)))
            rFrgs -= 1 + bContain;

        for(FragmentEnds::const_iterator sIter = smalls.begin(); sIter != smalls.end(); sIter++) {
            UnitigBreakPoint small = *sIter;

            if (small.fragEnd.id == big.fragEnd.id)
                continue; 

            bool rev  = isReverse(small.position);
            int sContain = 0;
            int lFrgs = small.fragNumber;
            iuid sid  = small.fragEnd.id;

            // use middle of frag instead of end, to get some overlap
            double bp = (small.position.end + small.position.bgn) / 2.0;

            if (cntnrmap_ptr->find(sid) != cntnrmap_ptr->end())
                sContain = (*cntnrmap_ptr)[sid].size();
            
            // left side of the frag in the unitig, don't count it
            if (((rev == true)  && (small.fragEnd == THREE_PRIME)) ||
                ((rev == false) && (small.fragEnd == FIVE_PRIME)))
                lFrgs -= 1 + sContain;

            if (rFrgs - lFrgs == 1)
                continue;

            double lRate = (lFrgs - lastBPFragNum) / (bp - lastBPCoord); 
            double rRate = (rFrgs - lFrgs) / (right - bp);
            double ratio = (lRate > rRate) ? lRate / rRate : rRate / lRate;
            double diff  = fabs( lRate - rRate);

            if ((ratio > 1.8) && (diff > difference)) {
                //fprintf(stderr, "Break frg %7d b %4d l %4d pos b %5d e %5.0f lRate %.4f\n", sid, lastBPFragNum, lFrgs, lastBPCoord, bp, lRate );
                //fprintf(stderr, "     diff %4d r %4d pos %5d rRate %.4f ratio %.2f to frag %7d\n", rFrgs - lFrgs, rFrgs, right, rRate, ratio, big.fragEnd.id); 
                //fprintf(stderr,"     select frg %d for break on arrival rate diff %.4f\n", sid, diff);
                difference = diff;
                selection = small;
            }
        }

        lastBPCoord   = right;
        lastBPFragNum = rFrgs;

        return selection;
    }



    void UnitigGraph::filterBreakPoints(Unitig *tig,
                                        FragmentEnds &breaks) {

        int lastBPCoord = 0;
        int lastBPFragNum = 0;

        // Check for sentinal at end, not added by MateChecker
        UnitigBreakPoint fakeEnd = breaks.back();
        if (fakeEnd.inSize == std::numeric_limits<int>::max())
            breaks.pop_back();

        FragmentEnds smallBPs;
        FragmentEnds newBPs;

        bool hadBig = false;

        for(FragmentEnds::iterator iter = breaks.begin(); iter != breaks.end(); iter++) {
            UnitigBreakPoint nextBP = *iter;

            if ((nextBP.inFrags > 1) && (nextBP.inSize > 500)) {
                hadBig = true;

                // big one, compare against smalls -- if we haven't
                // seen this break point already.
                //
                if (newBPs.empty() || (nextBP.fragEnd != newBPs.back().fragEnd)) {

                    if (smallBPs.empty() == false) {
                        //  smalls exist, select one.
                        UnitigBreakPoint small = selectSmall(tig, smallBPs, nextBP, lastBPCoord, lastBPFragNum);
                        if (small.fragNumber > 0)
                            newBPs.push_back(small);
                        smallBPs.clear();

                    } else {
                        //  No smalls.  Update state to move past the current big breakpoint.
                        lastBPCoord   = (nextBP.fragEnd.end == FIVE_PRIME) ? nextBP.position.bgn : nextBP.position.end;
                        lastBPFragNum = nextBP.fragNumber;

                        int bContain = 0;

                        if ((cntnrmap_ptr != NULL) && (cntnrmap_ptr->find(nextBP.fragEnd.id) != cntnrmap_ptr->end()))
                            bContain = (*cntnrmap_ptr)[nextBP.fragEnd.id].size();

                        bool bRev = isReverse(nextBP.position);

                        if (((bRev == true)  && (nextBP.fragEnd == THREE_PRIME)) ||
                            ((bRev == false) && (nextBP.fragEnd == FIVE_PRIME)))
                            lastBPFragNum -= 1 + bContain;
                    }

                    newBPs.push_back( nextBP );
                }
            } else {
                //  Not a big breaker.  Save the small breaker if we've not seen it yet.
                if (newBPs.empty() || ((nextBP.fragEnd != newBPs.back().fragEnd) &&
                                       (nextBP.fragEnd != smallBPs.back().fragEnd)))
                    smallBPs.push_back(nextBP);
            }
        }

        //  If we've got small ones saved, select one.

        //  XXX  This used to only get triggered if hadBig was true!

        if ((hadBig) && (smallBPs.empty() == false)) {
            UnitigBreakPoint small = selectSmall(tig, smallBPs, fakeEnd, lastBPCoord, lastBPFragNum);
            if (small.fragNumber > 0)
                newBPs.push_back(small);
            smallBPs.clear();
        }

        breaks = newBPs;
    }

    // Doesn't handle contained frags yet
    UnitigVector* UnitigGraph::breakUnitigAt(Unitig *tig, FragmentEnds &breaks) {

        if (breaks.empty())
            return NULL;

        // remove small break points
        filterBreakPoints(tig, breaks);

        // if we filtered all the breaks out return an empty list
        if (breaks.empty()) 
            return NULL;

        //  This is to catch a bug in....something before us.
        //  Sometimes the breakpoint list contains nonsense that tells
        //  us to split a unitig before the first fragment, or after
        //  the last.
        //
        //  That code is pretty dense (and big) so instead we'll make
        //  one pass through the list of breaks and see if we'd do any
        //  splitting.
        //
        //  It's a huge code duplication of the main loop.
        //
        {
            FragmentEnds breakstmp = breaks;

            UnitigBreakPoint breakPoint = breakstmp.front();
            breakstmp.pop_front();

            UnitigBreakPoint nextBP;

            int   newUnitigsConstructed = 0;
            int   newUnitigExists = 0;

            if (!breakstmp.empty())
                nextBP = breakstmp.front();

            for(DoveTailIter dtIter = tig->dovetail_path_ptr->begin(); dtIter != tig->dovetail_path_ptr->end(); dtIter++) {
                DoveTailNode frg = *dtIter;

                bool bothEnds = false;
                while ( !breakstmp.empty() && nextBP.fragEnd.id == breakPoint.fragEnd.id ) {
                    if (nextBP.fragEnd.end != breakPoint.fragEnd.end)
                        bothEnds = true;
                    breakstmp.pop_front();
                    if (!breakstmp.empty())
                        nextBP = breakstmp.front();
                }

                bool reverse = isReverse(frg.position);

                if (breakPoint.fragEnd.id == frg.ident) {

                    if (bothEnds) {
                        newUnitigsConstructed++;
                        newUnitigExists = 0;
                    }

                    else if (breakPoint.fragEnd.end ==  FIVE_PRIME && !reverse ||
                             breakPoint.fragEnd.end == THREE_PRIME && reverse) {
                        newUnitigsConstructed++;
                        newUnitigExists = 1;
                    }

                    else if (breakPoint.fragEnd.end ==  FIVE_PRIME && reverse ||
                             breakPoint.fragEnd.end == THREE_PRIME && !reverse ) {
                        if (newUnitigExists == 0)
                            newUnitigsConstructed++;
                        newUnitigExists = 0;
                    } else {
                        assert(0);
                    }

                    if (breakstmp.empty()) {
                        breakPoint.fragEnd.id = 0;
                    } else {
                        breakPoint = breakstmp.front();
                        breakstmp.pop_front();
                        if (!breakstmp.empty())
                            nextBP = breakstmp.front();
                    }
                } else {
                    if (newUnitigExists == 0) {
                        newUnitigsConstructed++;
                        newUnitigExists = 1;
                    }
                }
            }

            if (newUnitigsConstructed < 2) {
                fprintf(stderr, "SPLITTING BUG DETECTED!  Adjusting for it.\n");
                return NULL;
            }

        }  //  end of explicit split test
        



        UnitigVector *splits = new UnitigVector();
        Unitig       *newTig = NULL;
        int           offset = 0;

        UnitigBreakPoint breakPoint = breaks.front();
        breaks.pop_front();

        UnitigBreakPoint nextBP;

        if (!breaks.empty())
            nextBP = breaks.front();

        for(DoveTailIter dtIter = tig->dovetail_path_ptr->begin(); dtIter != tig->dovetail_path_ptr->end(); dtIter++) {
            DoveTailNode frg = *dtIter;

            // reduce multiple breaks at the same fragment end down to one
            bool bothEnds = false;
            while ( !breaks.empty() && nextBP.fragEnd.id == breakPoint.fragEnd.id ) {
                if (nextBP.fragEnd.end != breakPoint.fragEnd.end)
                    bothEnds = true;
                breaks.pop_front();
                if (!breaks.empty())
                    nextBP = breaks.front();
            }

            bool reverse = isReverse(frg.position);

            if (breakPoint.fragEnd.id == frg.ident) {

                //  There is method to this madness.  The three if
                //  blocks below have the same two basic code blocks
                //  in various orders.
                //
                //  o One code block makes a new unitig.
                //
                //  o The other adds a fragment to it, possibly
                //    remembering the beginning offset.
                //
                //  Finally, we delay making a new unitig until we are
                //  absolutely sure we're going to put fragments into
                //  it.
                //

                if (bothEnds) {
                    //
                    //  Break at both ends, create a singleton for this fragment.
                    //
                    fprintf(stderr,"  Breaking tig %d at both ends of %d num %d\n", tig->id(), breakPoint.fragEnd.id, breakPoint.fragNumber);

                    newTig = new Unitig();  //  always make a new unitig, we put a frag in it now
                    splits->push_back( newTig );

                    if (newTig->dovetail_path_ptr->empty())
                        offset = reverse ? -frg.position.end : -frg.position.bgn;
                    newTig->addFrag( frg, offset );

                    newTig = NULL;  //  delay until we need to make it
                }

                else if (breakPoint.fragEnd.end ==  FIVE_PRIME && !reverse ||
                         breakPoint.fragEnd.end == THREE_PRIME && reverse) {
                    //
                    //  Break at left end of frg, frg starts new tig
                    //
                    fprintf(stderr,"  Break tig %d before %d num %d\n", tig->id(), breakPoint.fragEnd.id, breakPoint.fragNumber);

                    newTig = new Unitig();  //  always make a new unitig, we put a frag in it now
                    splits->push_back( newTig );

                    if (newTig->dovetail_path_ptr->empty())
                        offset = reverse ? -frg.position.end : -frg.position.bgn;
                    newTig->addFrag( frg, offset );
                }

                else if (breakPoint.fragEnd.end ==  FIVE_PRIME && reverse ||
                         breakPoint.fragEnd.end == THREE_PRIME && !reverse ) {
                    //
                    //  Break at right end of frg, frg goes in existing tig, then make new unitig
                    //

                    if (newTig == NULL) {  //  delayed creation?
                        newTig = new Unitig();
                        splits->push_back( newTig );
                    }

                    if (newTig->dovetail_path_ptr->empty())
                        offset = reverse ? -frg.position.end : -frg.position.bgn;
                    newTig->addFrag( frg, offset );

                    fprintf(stderr,"  Break tig %d after %d num %d\n", tig->id(), breakPoint.fragEnd.id, breakPoint.fragNumber);

                    //  Delay making a new unitig until we need it.
                    newTig = NULL;
                } else {
                    // logically impossible!
                    assert(0);
                }


                // Done breaking, continue adding remaining frags to newTig
                if (breaks.empty()) {
                    breakPoint.fragEnd.id = 0;
                } else {
                    breakPoint = breaks.front();
                    breaks.pop_front();
                    if (!breaks.empty())
                        nextBP = breaks.front();
                }
            } else {
                //
                //  frag is not a unitig break point fragment, add to current new tig
                //

                if (newTig == NULL) {  //  delayed creation?
                    newTig = new Unitig();
                    splits->push_back( newTig );
                }
                if (newTig->dovetail_path_ptr->empty())
                    offset = reverse ? -frg.position.end : -frg.position.bgn;
                newTig->addFrag( frg, offset );
            }
        }

        return splits;
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
        if(aleft!=bleft) {
            return(aleft - bleft);
        }
        else if (aright != bright) {
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
                if (containPartialOrder == NULL)
                    return(0);
                else
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
        qsort( &(dovetail_path_ptr->front()), getNumFrags(), sizeof(IntMultiPos), &IntMultiPosCmp );
    }
    //////////////////////////////////////////////////////////////////////////////

    void UnitigGraph::writeIUMtoFile(char *fileprefix, int fragment_count_target){
      int         fragment_count         = 0;
      int         file_count             = 1;
      char        filename[FILENAME_MAX] = {0};
      int         iumiid                 = 0;
      GenericMesg mesg;

      // Open up the initial output file

      sprintf(filename, "%s_%03d.cgb", fileprefix, file_count++);
      FILE *file = fopen(filename,"w");
      assert(NULL != file);

      sprintf(filename, "%s.iidmap", fileprefix, file_count++);
      FILE *iidm = fopen(filename,"w");
      assert(NULL != iidm);

      // Step through all the unitigs

      for (UnitigVector::iterator utg_itr=unitigs->begin(); utg_itr != unitigs->end(); utg_itr++) {
        Unitig        *utg = *utg_itr;

        if (utg == NULL)
          continue;

        if (utg->getNumFrags() == 0)
          continue;

        IntUnitigMesg *ium_mesg_ptr = new IntUnitigMesg;

        ium_mesg_ptr->iaccession    = iumiid++;
#ifdef AS_ENABLE_SOURCE
        ium_mesg_ptr->source        = "gen> @@ [0,0]";
#endif
        ium_mesg_ptr->coverage_stat = utg->getCovStat();
        ium_mesg_ptr->status        = AS_UNASSIGNED;
        ium_mesg_ptr->unique_rept   = AS_FORCED_NONE;
        ium_mesg_ptr->length        = utg->getLength();
        ium_mesg_ptr->consensus     = "";
        ium_mesg_ptr->quality       = "";
        ium_mesg_ptr->forced        = 0;
        ium_mesg_ptr->num_frags     = utg->getNumFrags();
        ium_mesg_ptr->f_list        = &(utg->dovetail_path_ptr->front());

        fprintf(iidm, "Unitig %d == IUM %d\n", (*utg_itr)->id(), ium_mesg_ptr->iaccession);

        fragment_count += ium_mesg_ptr->num_frags;

        if ((fragment_count_target >= 0) &&
            (fragment_count >= fragment_count_target)) {
          fclose(file);
          sprintf(filename, "%s_%03d.cgb", fileprefix, file_count++);
          file = fopen(filename,"w");
          assert(NULL != file);
          fragment_count = ium_mesg_ptr->num_frags;
        }

        mesg.m = ium_mesg_ptr;
        mesg.t = MESG_IUM;

        WriteProtoMesg_AS(file, &mesg);

        delete ium_mesg_ptr;
      }

      fclose(file);
      fclose(iidm);
    }

    //////////////////////////////////////////////////////////////////////////////

    void UnitigGraph::readIUMsFromFile(const char *filename, iuid maxIID){
        FILE *iumIn = fopen( filename, "r");
        if (errno) {
            fprintf(stderr, "Could not open '%s' for input: %s\n",
                    filename, strerror(errno));
            exit(1);
        }
        Unitig::resetFragUnitigMap(maxIID+1);
        BestOverlapGraph::fragLength = new uint16[maxIID+1];
        memset( BestOverlapGraph::fragLength, std::numeric_limits<uint16>::max(),
                sizeof(uint16)*(maxIID+1) );

        GenericMesg *pmesg = NULL;
        while ((ReadProtoMesg_AS( iumIn, &pmesg ) != EOF)) {

            if (pmesg->t == MESG_IUM) {
                IntUnitigMesg *intig = (IntUnitigMesg *)(pmesg->m);
                Unitig *tig = new Unitig(intig->iaccession);
                unitigs->push_back( tig );

                for(int i=0; i < intig->num_frags; i++) {
                    IntMultiPos frg = intig->f_list[ i ];
                    tig->addFrag( frg );
                    BestOverlapGraph::fragLength[ frg.ident ]= abs(frg.position.end - frg.position.bgn);

                    OVSoverlap olap;
                    memset( static_cast<void*>(&olap), 0, sizeof(olap));
                    olap.a_iid = frg.contained;
                    olap.b_iid = frg.ident;
                    bog_ptr->setBestContainer( olap, 0 );
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////

    BestEdgeCounts UnitigGraph::countInternalBestEdges( const Unitig *tig) {
        DoveTailConstIter iter = tig->dovetail_path_ptr->begin();
        DoveTailNode prevFrag;
        BestEdgeCounts cnts;
        for(; iter != tig->dovetail_path_ptr->end(); iter++) {
            DoveTailNode frag = *iter;
            if (iter == tig->dovetail_path_ptr->begin()){
                prevFrag = frag;
                continue;
            }
            if (frag.contained) {
                cnts.contained++;
                continue;
            }

            BestEdgeOverlap* prevBestEdge;
            if( isReverse( prevFrag.position ) )
                prevBestEdge = bog_ptr->getBestEdgeOverlap( prevFrag.ident, FIVE_PRIME);
            else
                prevBestEdge = bog_ptr->getBestEdgeOverlap( prevFrag.ident, THREE_PRIME);

            BestEdgeOverlap* bestEdge;
            if( isReverse( frag.position ) )
                bestEdge = bog_ptr->getBestEdgeOverlap( frag.ident, THREE_PRIME);
            else
                bestEdge = bog_ptr->getBestEdgeOverlap( frag.ident, FIVE_PRIME);

            if ( prevBestEdge->frag_b_id == frag.ident ) {
                // Either one way best or dovetail
                if (bestEdge->frag_b_id == prevFrag.ident) // dovetail
                    cnts.dovetail++;
                else
                    cnts.oneWayBest++;
            } else {
                // Either one way best or neither
                if (bestEdge->frag_b_id == prevFrag.ident)
                    cnts.oneWayBest++;
                else
                    cnts.neither++;
            }
            prevFrag = frag;
        }
        return cnts;
    }

    //////////////////////////////////////////////////////////////////////////////

    BestEdgeCounts UnitigGraph::countInternalBestEdges() {
        BestEdgeCounts allTigs;
        UnitigsConstIter iter = unitigs->begin();
        for(iter = unitigs->begin(); iter != unitigs->end(); iter++){
            if( *iter == NULL)
                continue;

            if( (*iter)->getNumFrags() == 0 )
                continue;

            Unitig* tig = *iter;
            BestEdgeCounts cnts = countInternalBestEdges(tig);
            //std::cerr << "Tig " << tig->id() << " overlap counts, dovetail " << cnts.dovetail
            //          << " oneWayBest " << cnts.oneWayBest << " neither " << cnts.neither
            //          << " contained " << cnts.contained << std::endl;
            allTigs += cnts;
        }
        return allTigs;
    }

    //////////////////////////////////////////////////////////////////////////////

} //AS_BOG namespace
