
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
* Module:  AS_BOG_BestOverlapGraph.cc
* Description:
*    Data structure to contain the best overlaps and containments
*    based on a defined metric.
* 
*    Programmer:  K. Li
*       Started:  1 August 2005
*
*   Implementation: E. Venter
*   Commented:    K. Li 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: AS_BOG_BestOverlapGraph.cc,v 1.21 2005-10-31 15:55:57 eliv Exp $
 * $Revision: 1.21 $
*/

static const char CM_ID[] = "$Id: AS_BOG_BestOverlapGraph.cc,v 1.21 2005-10-31 15:55:57 eliv Exp $";

//  System include files
#include<iostream>

#include "AS_BOG_BestOverlapGraph.hh"
//#include "AS_BOG_BestOverlapGraphVisitor.hh"

extern "C" {
#include "AS_PER_fragStore.h"
}

namespace AS_BOG{

    ///////////////////////////////////////////////////////////////////////////
    //
    //  Example:
    //
    //  The overlap, pi, exists between A and B:
    //
    //  A -------------->
    //         |||||||||
    //  B      ---------------->
    //
    //  AEnd(pi) is 3'
    //  BEnd(pi) is 5'
    //

    fragment_end_type BestOverlapGraph::AEnd(const Long_Olap_Data_t& olap) {
        if (olap.a_hang < 0 && olap.b_hang < 0)
            return FIVE_PRIME;
        if (olap.a_hang > 0 && olap.b_hang > 0)
            return THREE_PRIME;

        assert(0); // no contained
    }

    fragment_end_type BestOverlapGraph::BEnd(const Long_Olap_Data_t& olap) {
        if (olap.a_hang < 0 && olap.b_hang < 0)
            if ( olap.flipped )
                return FIVE_PRIME;
            else
                return THREE_PRIME;

        if (olap.a_hang > 0 && olap.b_hang > 0)
            if ( olap.flipped )
                return THREE_PRIME;
            else
                return FIVE_PRIME;

        assert(0); // no contained
    }

    ///////////////////////////////////////////////////////////////////////////
    // Constructor

    // Create BestOverlapGraph as an array of size max fragments.
    //     Assuming that our iuids start at index value of 1.

    BestOverlapGraph::BestOverlapGraph(int max_fragments)
        : _num_fragments(max_fragments), curFrag(0)
    {
        _best_overlaps = new BestFragmentOverlap[max_fragments+1];

        memset(_best_overlaps, 0, sizeof(BestFragmentOverlap)*(max_fragments+1));
    }

    // Destructor
    BestOverlapGraph::~BestOverlapGraph(){
        delete[] _best_overlaps;
    }

    // Interface to graph visitor
//    void BestOverlapGraph::accept(BestOverlapGraphVisitor bog_vis){
//        bog_vis.visit(this);
//    }

    ///////////////////////////////////////////////////////////////////////////
    // Accessor Get Functions

    //  Given a fragment IUID and which end, returns pointer to
    //  BestOverlap node.  
    BestEdgeOverlap *BestOverlapGraph::getBestEdgeOverlap(
        iuid frag_id, fragment_end_type which_end){

        if(which_end == FIVE_PRIME)
            return(&_best_overlaps[frag_id].five_prime);
        else if(which_end == THREE_PRIME){
            return(&_best_overlaps[frag_id].three_prime);
        }
    }

    //  Given an overlap, determines which record (iuid and end) and sets the newScore.
    void BestOverlapGraph::setBestEdgeOverlap(const Long_Olap_Data_t& olap, float newScore) {

        if (AEnd(olap) == THREE_PRIME) {
            _best_overlaps[ olap.a_iid ].three_prime.frag_b_id = olap.b_iid;
            _best_overlaps[ olap.a_iid ].three_prime.score     = newScore;
            _best_overlaps[ olap.a_iid ].three_prime.bend      = BEnd(olap);
            _best_overlaps[ olap.a_iid ].three_prime.olap_len  = olapLength(olap);

        }
        if (AEnd(olap) == FIVE_PRIME) {
            _best_overlaps[ olap.a_iid ].five_prime.frag_b_id = olap.b_iid;
            _best_overlaps[ olap.a_iid ].five_prime.score     = newScore;
            _best_overlaps[ olap.a_iid ].five_prime.bend      = BEnd(olap);
            _best_overlaps[ olap.a_iid ].five_prime.olap_len  = olapLength(olap);
        }
    }

    bool BestOverlapGraph::isContained(const iuid fragIID) {
        std::map<iuid,BestContainment>::iterator i =
                                _best_containments.find( fragIID ); 
        if ( i != _best_containments.end() ) 
            return true;
        else
            return false;
    }
    // Given a containee, returns pointer to BestContainment record
    BestContainment* BestOverlapGraph::getBestContainer(iuid containee)
    {
        if ( isContained( containee ) ) 
            return &_best_containments[containee];
        else
            return NULL;
    }

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    // Graph creation methods

    ///////////////////////////////////////////////////////////////////////////
    // Transitively removes redundant containments, so all containees in a container, refer
    // to the same container.  Algorithm will go through each element in the list of contained
    // fragments, and then for each element follow each container's container.
    void BestOverlapGraph::removeTransitiveContainment() {

        // Loop through each containee that has been stored in _best_containments
        for(std::map<iuid,BestContainment>::const_iterator it = _best_containments.begin();
                it != _best_containments.end(); it++)
        {
            iuid id = it->first;                     // Gets the iuid of the containee under analysis
            BestContainment bst = it->second;        // Gets the BestContainment record of the containee under analysis

            // Get container information based on containee id/containment record
            bool sameOrient = bst.sameOrientation;
            bool useAhang = true;

            // For this std::map:
            //   iuid is the containee iuid
            //   BestContainment.container is the best container iuid
            //       This finds the iterator if the containerOf(id) is contained, else returns
            //       _best_containments.end() if not found.
            //
            // When this find returns, the i2 will point to id's container's container
            std::map<iuid,BestContainment>::iterator i2 =
                                        _best_containments.find( bst.container);

            // Keep track of which container's we've looked at for each transitive path starting from
            //   the containee under analysis.
            std::map<iuid,BestContainment> found;
            found[bst.container] = bst;

            // Loop while the current container is a containee of another container 
            while ( i2 != _best_containments.end() ) {
                BestContainment nb = i2->second;
                std::cout << id <<" "<<bst.container<<" "<< nb.container<< std::endl;

                // Delete containee under analysis from _best_containments if its container is concontained
                //   by itself.  ie if id's container is contained by id.  This eliminates ciruclar containment.
                if ( nb.container == id ) {
                    _best_containments.erase( id );
                    std::cout << "Erase self" << std::endl;
                    break;
                }

                // Look through the list of containers we've already walked past to find
                //   circular containment greater than one degree away.
                std::map<iuid,BestContainment>::iterator seen= found.find( nb.container );

                // If i2's container has already been traversed
                if ( seen != found.end() ) { 

                    // Set id's container to the larger container.
                    _best_containments[ id ] = seen->second;
                    std::cout << "Circled " << seen->second.container<< std::endl;

                    // Remove the container of id's new larger container
                    _best_containments.erase( seen->second.container );
                    break;
                }

                // Set id's container to id's container's container, in this case
                //   one transitive step away.
                _best_containments[ id ] = nb;

                // Keep track of what containers have been used.
                found[ nb.container ] = nb;

                // Reset the orientation of id's containment.
                if (!sameOrient) {
                    sameOrient = _best_containments[id].sameOrientation = ! nb.sameOrientation;
                    useAhang = !useAhang;
                }
                found[nb.container].sameOrientation = sameOrient;
                // make hang relative to new container
                _best_containments[ id ].a_hang += useAhang ? bst.a_hang : abs(bst.b_hang); 

                // Iterator to the i2's container
                i2 = _best_containments.find( nb.container);
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    // Array of cached fragment lengths.  Initialized to 0.  Use BestOverlapGraph::fragLen(iuid)
    //   to access contents b/c method will read info from FragStore and populate
    //   record if not already read in.
    uint16         *BestOverlapGraph::fragLength;                        

    // Frag Store related data structures
    ReadStructp     BestOverlapGraph::fsread = new_ReadStruct();
    FragStoreHandle BestOverlapGraph::fragStoreHandle;

    uint16 BestOverlapGraph::fragLen( iuid iid ) {
        // If fragLength is not already cached, compute it after reading it
        //   in from the fragStore, store it and return it to the caller.
        if (BestOverlapGraph::fragLength[ iid ] == 0) {
            uint32 clrBgn, clrEnd;
            getFragStore( fragStoreHandle, iid, FRAG_S_SEQUENCE, fsread);
            getClearRegion_ReadStruct( fsread, &clrBgn, &clrEnd, READSTRUCT_LATEST);
            BestOverlapGraph::fragLength[ iid ] = clrEnd - clrBgn;
        }
        return BestOverlapGraph::fragLength[ iid ];
    }

    ///////////////////////////////////////////////////////////////////////////

    short BestOverlapGraph::olapLength(const Long_Olap_Data_t& olap) {
        // Computes overlap length given the amount of overhang in the
        //   overlap record and fragment length (which is not in the
        //   overlap record).
        uint16 alen = fragLen(olap.a_iid);
        if (olap.a_hang < 0)
            return alen - abs(olap.b_hang);
        else
            return alen - olap.a_hang;
    }

    ///////////////////////////////////////////////////////////////////////////

    void BestOverlapGraph::updateInDegree() {
        if ( ! isContained(curFrag) ) {
            // Update B's in degree on A's 3' End
            iuid bid = _best_overlaps[ curFrag ].three_prime.frag_b_id;
            switch(_best_overlaps[ curFrag ].three_prime.bend){
                case THREE_PRIME:
                    _best_overlaps[ bid ].three_prime.in_degree++; break;
                case FIVE_PRIME:
                    _best_overlaps[ bid ].five_prime.in_degree++; break;
                default: assert(0);
            }

            // Update B's in degree on A's 5' End
            bid = _best_overlaps[ curFrag ].five_prime.frag_b_id;
            switch(_best_overlaps[ curFrag ].five_prime.bend){
                case THREE_PRIME:
                    _best_overlaps[ bid ].three_prime.in_degree++; break;
                case FIVE_PRIME:
                    _best_overlaps[ bid ].five_prime.in_degree++; break;
                default: assert(0);
            }
        }
    }
    bool BestOverlapGraph::checkForNextFrag(const Long_Olap_Data_t& olap) {
        // Update the in_degrees whenever the incoming overlap's A's Fragment IUID changes.
        //   Returns true, if the olap's A fragment ID has changed since the previous call,
        //   else false.
        //
        // This method is only called by processOverlap, should make it private.
        //
        // Since the overlaps are coming in in order of A's iuid, 
        //   its safe for us to assume that once the incoming overlap A iuid changes
	//   that we are completely done processing overlaps from A's
        //   fragment.   This implies that we know that A's best overlap
        //   at this point will not change, so it's safe to update it's overlapping
        //   partner's (B's), in degree for both ends of A.
        
	// In this code, the olap.a_iid is considered the "next" frag, so 
	//   curFrag is the IID of the fragment prior to receiving this olap.
        if (curFrag != olap.a_iid) {

            // update in degree if not contained
            updateInDegree();

            // Set up the curFrag to refer to the incoming (next) fragment IUID	
            curFrag = olap.a_iid;

            // Initialize the overlap score for A
            _best_overlaps[ olap.a_iid ].three_prime.score = 0;
            _best_overlaps[ olap.a_iid ].five_prime.score = 0;
            bestLength = 0;

            // Means that A's fragment ID has changed since the previous
            return true;
        }

	// Means that A's fragment ID is the same
        return false;
    }
    void BestOverlapGraph::scoreContainment(const Long_Olap_Data_t& olap) {
        // This function builds the Containment graph by considering the specified
        //   overlap record.  It's important to remember that the overlaps
        //   must be passed in in sorted order of fragment A's iuid.  Also,
        //   it is also very important that overlaps are doubly listed, ie.
        //   if an overlap exists, an olap must be passed in where A overlaps with
        //   B and an olap must later be passed in where B overlaps A.
        //
        // It calls the virtual score function to score the overlap,
        //   determines whether the overlap is containment or dovetailing,
        //   then stores the overlap in the BestOverlapGraph member variables.

        // Compute the score for this overlap based on the virtual score function
        float newScr = scoreOverlap(olap);

        // If the score is 0, the overlap doesn't pass the scoring criteria at all
        //   so don't store the overlap whether or not it's dovetailing or containment.
        if ( newScr <= 0 )
            return;

        // A contains B
        if ( olap.a_hang >= 0 && olap.b_hang <= 0 ) {
             //handle a contains b

             //  Get the pointer to the container for the containee, B, if it exists
             BestContainment *best = getBestContainer( olap.b_iid );

             //  Set the containee, B, to the container, A, if any of the following conditions exist:
             //     1.) a containment relationship/record doesn't exist
             //     2.) the new relationship has a higher score
             //     3.) Container fragment ID is nominally less than the existing container ID
             if (NULL == best || newScr > best->score || newScr == best->score
                      && fragLen(best->container) < fragLen(olap.a_iid) )
             {
//                 std::cout << olap.a_iid << " contains " << olap.b_iid <<" "<< fragLen(olap.a_iid) << std::endl; 
                 BestContainment newBest;
                 newBest.container = olap.a_iid;
                 newBest.score     = newScr;
                 newBest.sameOrientation =  olap.flipped ? false : true;
                 newBest.a_hang =  olap.a_hang;
                 newBest.b_hang =  olap.b_hang;
                 _best_containments[ olap.b_iid ] = newBest;
             }
        
        // B contains A, this code is commented out/unnecessary because the overlap
        //   recorded is listed twice in the overlap store.
        } else if ( olap.a_hang <= 0 && olap.b_hang >= 0 ) {
             //handle b contains a
             
        //  Dove tailing overlap
        } else {
        }
    }
    void BestOverlapGraph::scoreEdge(const Long_Olap_Data_t& olap) {
        // This function builds the BestOverlapGraph by considering the specified
        //   overlap record.  It's important to remember that the overlaps
        //   must be passed in in sorted order of fragment A's iuid.  Also,
        //   it is also very important that overlaps are doubly listed, ie.
        //   if an overlap exists, an olap must be passed in where A overlaps with
        //   B and an olap must later be passed in where B overlaps A.
        //
        // It calls the virtual score function to score the overlap,
        //   determines whether the overlap is containment or dovetailing,
        //   then stores the overlap in the BestOverlapGraph member variables.

        // Compute the score for this overlap based on the virtual score function
        float newScr = scoreOverlap(olap);

        // If the score is 0, the overlap doesn't pass the scoring criteria at all
        //   so don't store the overlap whether or not it's dovetailing or containment.
        if ( newScr <= 0 )
            return;

        // A contains B
        if ( olap.a_hang >= 0 && olap.b_hang <= 0 ) {
             //handle a contains b

        // B contains A, this code is commented out/unnecessary because the overlap
        //   recorded is listed twice in the overlap store.
        } else if ( olap.a_hang <= 0 && olap.b_hang >= 0 ) {
             //handle b contains a
             
        //  Dove tailing overlap
        } else {
             // no containment, so score

             // Update the in degree for what A overlaps with if the current
             //   A fragment is different than the previous one.
             checkForNextFrag(olap);

             BestEdgeOverlap *best = getBestEdgeOverlap( olap.a_iid, AEnd(olap));
             short olapLen = olapLength(olap);

             // Store the overlap if:
             //   1.)  The score is better than what is already in the graph
             //   2.)  If the scores are identical, the one with the longer length
             //
             // Since the order of how the overlaps are read in from the overlap
             //   store are by A's increasing iuid, by default, if the score
             //   and length are the same, the iuid of the lower value will be
             //   kept.

             if (newScr > best->score || newScr == best->score &&
                olapLen > bestLength )
             {
                 setBestEdgeOverlap( olap, newScr );
                 bestLength = olapLen;
             }
         }
    }

    ///////////////////////////////////////////////////////////////////////////

    void BestOverlapGraph::printFrom(iuid begin, iuid end){
        iuid i;

	end=(end==0)?begin:end;

	for(i=begin; i<=end; i++){
	    std::cout << 
		_best_overlaps[i].five_prime.frag_b_id << 
		"<-" << i << "->" <<
	        _best_overlaps[i].three_prime.frag_b_id << 
		"\t" <<  _best_overlaps[i].five_prime.in_degree <<
		"/" <<  _best_overlaps[i].three_prime.in_degree <<
		std::endl;
	}
    }
    ///////////////////////////////////////////////////////////////////////////

    float ErateScore::scoreOverlap(const Long_Olap_Data_t& olap) {
    // Computes the score for a Error Rate BOG based on overlap corrected error rate.  
    //        Error rate is normalized so that the higher the error rate, the lower the score.

        return 100 - Expand_Quality(olap.corr_erate) * 100;
    }
    
    float LongestEdge::scoreOverlap(const Long_Olap_Data_t& olap) {
    // Computes the score for a Longest Edge BOG based on overlap length only.

        return olapLength(olap);
    }

    float LongestHighIdent::scoreOverlap(const Long_Olap_Data_t& olap) {
    // Computes the score for a Longest Edge BOG based on overlap length but
    //   after applying an an error rate cutoff.

        short olapLen = olapLength(olap);
        float erate = Expand_Quality(olap.corr_erate) * 100;

        // mismatchCuoff is specified in the constructor
        if (erate > mismatchCutoff)
            return 0;
        return olapLen;
    }

    ///////////////////////////////////////////////////////////////////////////

    void BOG_Runner::processOverlapStream(OVL_Store_t * my_store,OVL_Stream_t *my_stream) {

        // Go through the overlap stream in two passes:
        // first pass finds the containments
        // second pass builds the overlap graph, excluding contained frags
        
        Long_Olap_Data_t olap;
        while  (Next_From_OVL_Stream (&olap, my_stream))
        {
            for(int j = 0; j < metrics.size(); j++)
                metrics[j]->scoreContainment( olap );
        }
        for(int j = 0; j < metrics.size(); j++) {
            metrics[j]->removeTransitiveContainment();
        }
        Renew_OVL_Stream(my_stream);
        uint32 last = Last_Frag_In_OVL_Store( my_store );
        Init_OVL_Stream(my_stream, 1, last, my_store);

        while  (Next_From_OVL_Stream (&olap, my_stream))
        {
            for(int j = 0; j < metrics.size(); j++) {
                if (metrics[j]->isContained( olap.b_iid ) )
                    continue; // no contained frags in BestOverlapGraph

                metrics[j]->scoreEdge( olap );
            }
        }
        // Update degree on final frag
        for(int j = 0; j < metrics.size(); j++) {
            metrics[j]->updateInDegree();
        }
    }

} //AS_BOG namespace
