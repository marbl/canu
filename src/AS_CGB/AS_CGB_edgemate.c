
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
static char CM_ID[] 
= "$Id: AS_CGB_edgemate.c,v 1.4 2005-03-22 19:48:27 jason_miller Exp $";
/* *******************************************************************
 *
 * Module: AS_CGB_edgemate.c
 * 
 * Description: These routines find and access the mate directed edge
 * for a given edge of an overlap.
 *
 * Assumptions: 
 *
 * Author: Clark Mobarry
 *********************************************************************/

/*************************************************************************/
/* System include files */

/*************************************************************************/
/* Local include files */
#include "AS_CGB_all.h"

/*************************************************************************/
/* Conditional compilation */

//////////////////////////////////////////////////////////////////

void reflect_Aedge( Aedge *new_edge, Aedge *old_edge) {
/*  fragment overlaps:
       The overlapper connects two fragment-ends in an overlap
       relationship:
       
       A    ---------------->
       B          -------------->
       
       The direction mate edge preserves the which fragment-ends are
       in the overlap:
       
       B^c  <----------------
       A^c       <---------------- 
    */
  
  const IntFragment_ID old_avx = old_edge->avx;
  const int old_asx = old_edge->asx;
  const int old_ahg = old_edge->ahg;

  const IntFragment_ID old_bvx = old_edge->bvx;
  const int old_bsx = old_edge->bsx;
  const int old_bhg = old_edge->bhg;

  //set_crs_edge(edges,new_edge, old_edge);
  //set_crs_edge(edges,old_edge, new_edge);

  new_edge->avx = old_bvx;
  new_edge->asx = old_bsx;
  new_edge->ahg = old_bhg;
#ifdef STORE_OVERLAP_EXTREMES
  new_edge->amn = old_bmn;
  new_edge->amx = old_bmx;
#endif // STORE_OVERLAP_EXTREMES

  new_edge->bvx = old_avx;
  new_edge->bsx = old_asx;
  new_edge->bhg = old_ahg;
#ifdef STORE_OVERLAP_EXTREMES
  new_edge->bmn = old_amn;
  new_edge->bmx = old_amx;
#endif // STORE_OVERLAP_EXTREMES

  new_edge->nes = old_edge->nes;
  new_edge->blessed = old_edge->blessed; // Is this correct?
  new_edge->invalid = old_edge->invalid;
  new_edge->quality = old_edge->quality;
  new_edge->grangered = old_edge->grangered;
  new_edge->reflected = ! old_edge->reflected;
}   

void granger_Aedge( Aedge *new_edge, Aedge *old_edge) {
    /* To-contained edges have ahg>0
       A    ---------------->
       B          ------->...
       
       So create the granger mate edge:
       A^c  <----------------
       B^c     <-------......
    */
    /* From-contained edges bhg <0
       A    ..------->      
       B    ---------------->
       
       So create the granger mate edge:
       A^c  .......<-------
       B^c  <----------------
    */

  const IntFragment_ID old_avx = old_edge->avx;
  const int old_asx = old_edge->asx;
  const int old_ahg = old_edge->ahg;

  const IntFragment_ID old_bvx = old_edge->bvx;
  const int old_bsx = old_edge->bsx;
  const int old_bhg = old_edge->bhg;

  //set_crs_edge(edges,new_edge, old_edge);
  //set_crs_edge(edges,old_edge, new_edge);
  
  new_edge->avx =   old_avx;
  new_edge->asx = ! old_asx;
  new_edge->ahg = - old_bhg;
#ifdef STORE_OVERLAP_EXTREMES
  new_edge->amn = - old_bmx;
  new_edge->amx = - old_bmn;
#endif // STORE_OVERLAP_EXTREMES

  new_edge->bvx =   old_bvx;
  new_edge->bsx = ! old_bsx;
  new_edge->bhg = - old_ahg;
#ifdef STORE_OVERLAP_EXTREMES
  new_edge->bmn = - old_amx;
  new_edge->bmx = - old_amn;
#endif // STORE_OVERLAP_EXTREMES

  new_edge->nes = old_edge->nes;
  new_edge->blessed = old_edge->blessed; // Is this correct?
  new_edge->invalid = old_edge->invalid;
  new_edge->quality = old_edge->quality;
  new_edge->grangered = ! old_edge->grangered;
  new_edge->reflected =   old_edge->reflected;
}   

//////////////////////////////////////////////////////////////////

static QsortCompare the_compare_edge_function = NULL;

void set_compare_edge_function( QsortCompare compare_edge)
{
  the_compare_edge_function = compare_edge;
}

QsortCompare get_compare_edge_function(void)
{
  return the_compare_edge_function;
}

//////////////////////////////////////////////////////////////////

static IntEdge_ID find_overlap_edge_mate_binary_search
(/* Input Only */
 const Tfragment frags[], 
 const Tedge edges[],
 const IntEdge_ID ie0
 //const QsortCompare compare_edge
 //int (*compare_edge)(const void *,const void *)
 )
{
  // Find the other directed edge of the overlap and return the index.
  // Assert if the edge is not found.

  // This binary search version assumes that the edges are sorted with
  // respect to the compare_edge() function, and that all duplicate
  // edges with respect to the compare_edge() function have been
  // removed.

  const QsortCompare compare_edge = get_compare_edge_function();

  const IntEdge_ID nedge = GetNumEdges(edges);

  IntEdge_ID mate_edge = nedge;
  Aedge test_edge = *(GetVA_Aedge(edges,ie0));
  // For the uninitialized variable finder "man third".
  // I have cheated the accessor function rule here.

  assert(0 <= ie0);
  assert(ie0 < nedge);

  reflect_Aedge( &test_edge, &test_edge);
  
  {
    IntEdge_ID ie1 
      = get_segstart_vertex(frags,test_edge.avx,test_edge.asx);
    IntEdge_ID ie3 = ie1 
      + get_seglen_vertex(frags,test_edge.avx,test_edge.asx);
    IntEdge_ID ie2 = ie1;

    if(ie3 == ie1 ) { return AS_CGB_EDGE_NOT_FOUND;}
    // Empty edge segment.
    
    assert( ie1 < ie3);

    for(; (ie1 < ie3); ) {
      int icompare;
      
      ie2 = (ie1 + ie3)/2;

      {
        Aedge curr_edge = *(GetVA_Aedge(edges,ie2));

        icompare =
          compare_edge(&curr_edge,&test_edge);
        // This comparison function must satisfy the UNIX qsort semantics.
        
        // The comparison function will be called with two parameters that
        // point to the two elements to be compared. The comparison
        // function must return an integer less than, equal to, or greater
        // than zero, depending on whether the first element in the
        // comparison is considered less than, equal to, or greater than
        // the second element.
      }
      
      if( icompare == 0 ) {
	// Found the other half of the undirected edge?
	// const Tnes nes2 = get_nes_edge(edges,ie2);
        // found = (nes2 == nes0);
        mate_edge = ie2;
        goto escape;
      }

      assert((icompare!=0)&&(ie1 < ie3));
      if(icompare > 0) {
	ie3 = ie2;
      }
      if(icompare < 0) {
	ie1 = ie2+1;
      }
      
    }
  }
  mate_edge = AS_CGB_EDGE_NOT_FOUND;
 escape: ; // The edge is found...
  return mate_edge;
}


IntEdge_ID find_overlap_edge_mate
(/* Input Only */
 const Tfragment frags[], 
 const Tedge edges[],
 const IntEdge_ID ie0
 )
{
  return
    //find_overlap_edge_mate_linear_search
    find_overlap_edge_mate_binary_search
    ( frags, edges, ie0 );
}

void fix_overlap_edge_mate
(/* Input Only */
 const Tfragment frags[], 
       Tedge edges[],
 const IntEdge_ID ie0)
{
  const IntEdge_ID ie1 =
    find_overlap_edge_mate( frags, edges, ie0);
  // The index of the other half of the undirected edge.
  const Tnes ines0 = get_nes_edge(edges,ie0);
  Tnes ines1 = ines0; // This will be the "fix".

#if 0
  assert( is_a_dvt_edge(edges,ie0) );
#endif
  
  if(AS_CGB_EDGE_NOT_FOUND == ie1) {
    const IntFragment_ID avx = get_avx_edge(edges,ie0);
    const IntFragment_ID bvx = get_bvx_edge(edges,ie0);
    const IntFragment_ID aid = get_iid_fragment(frags,avx);
    const IntFragment_ID bid = get_iid_fragment(frags,bvx);
    fprintf(stderr,"ERROR: AS_CGB_EDGE_NOT_FOUND == ie1\n");
    fprintf(stderr,"ie0=" F_IID "\n",ie0);
    fprintf(stderr,"aid=" F_IID "\n",aid);
    fprintf(stderr,"avx=" F_IID "\n",avx);
    fprintf(stderr,"asx=%d\n",get_asx_edge(edges,ie0));
    fprintf(stderr,"bid=" F_IID "\n",bid);
    fprintf(stderr,"bvx=" F_IID "\n",bvx);
    fprintf(stderr,"bsx=%d\n",get_bsx_edge(edges,ie0));
    fprintf(stderr,"nes=%d\n",get_nes_edge(edges,ie0));
  }
  assert( AS_CGB_EDGE_NOT_FOUND != ie1 );
  // The return code for mate edge not found.
  
  set_nes_edge(edges,ie1,ines1);
}

static void fill_new_edge_with_reflected_old_edge
( Tedge * edges,
  const IntEdge_ID new_edge,
  const IntEdge_ID old_edge) 
{
  Aedge * the_new_edge = GetVA_Aedge(edges,new_edge);
  Aedge * the_old_edge = GetVA_Aedge(edges,old_edge);
  reflect_Aedge( the_new_edge, the_old_edge);
}

int verify_that_the_edges_are_in_order(Tedge edges[])
{
  const QsortCompare compare_edge = get_compare_edge_function();
  IntEdge_ID ie0 = 0, nedge = GetNumEdges(edges);
  
  for( ie0=0; ie0 < nedge-1; ie0++) {
    int icompare = compare_edge(GetVA_Aedge(edges,ie0),GetVA_Aedge(edges,ie0+1));
    assert(icompare <= 0);
    // Verify that the edges are in proper order.
  }
  return 0;
}

void append_the_edge_mates
(
 Tfragment frags[],
 Tedge edges[],
 TIntEdge_ID * next_edge_obj
)
{
  // Scan the current edges to find un-mated edges.  For each un-mated
  // edge, append the mate edge to the edge array.
  // QsortCompare compare_edge = get_compare_edge_function();

  IntEdge_ID 
    ie0,
    nedge = GetNumEdges(edges),
    nedge_delta = 0;
  fprintf(stderr,"append_the_edge_mates: nedge=" F_IID "\n", nedge);

  verify_that_the_edges_are_in_order(edges);

  for( ie0=0; ie0 < nedge; ie0++) {
    const IntEdge_ID ie1 = find_overlap_edge_mate( frags, edges, ie0);
    if( AS_CGB_EDGE_NOT_FOUND == ie1 ) {
      IntEdge_ID ie2 = nedge+nedge_delta;
      // fprintf(stderr,"nedge_delta=" F_IID "\n", nedge_delta);
      EnableRangeVA_Aedge(edges,ie2+1);
      fill_new_edge_with_reflected_old_edge( edges, ie2, ie0);
      nedge_delta ++;
    }
  }

  reorder_edges( frags, edges, next_edge_obj);
}


IntEdge_ID check_symmetry_of_the_edge_mates
(
 Tfragment frags[],
 Tedge edges[],
 TIntEdge_ID * next_edge_obj
)
{
  // const QsortCompare compare_edge = get_compare_edge_function();
  IntEdge_ID 
    ie0,
    nedge = GetNumEdges(edges),
    counter = 0;

  verify_that_the_edges_are_in_order(edges);

  for( ie0=0; ie0 < nedge; ie0++) {
    const IntEdge_ID ie1 = find_overlap_edge_mate( frags, edges, ie0);
    if( AS_CGB_EDGE_NOT_FOUND == ie1 ) {
      counter ++;
    } else {
      /*
      Aedge edge0 = *(GetVA_Aedge(edges,ie0));
      Aedge edge1 = *(GetVA_Aedge(edges,ie1));

      if(compare_edge( &edge0, &edge1) != 0) { }
      */
    }
  }
  fprintf(stderr,"check_symmetry_of_the_edge_mates: nedge=" F_IID " counter=" F_IID "\n",
          nedge, counter);
  return counter;
}

