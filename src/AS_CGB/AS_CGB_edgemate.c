
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

static char *rcsid = "$Id: AS_CGB_edgemate.c,v 1.12 2008-10-08 22:02:54 brianwalenz Exp $";

//  Description: These routines find and access the mate directed edge
//  for a given edge of an overlap.

#include "AS_CGB_all.h"
#include "AS_CGB_histo.h"

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

  new_edge->bvx = old_avx;
  new_edge->bsx = old_asx;
  new_edge->bhg = old_ahg;

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

  new_edge->bvx =   old_bvx;
  new_edge->bsx = ! old_bsx;
  new_edge->bhg = - old_ahg;

  new_edge->nes = old_edge->nes;
  new_edge->blessed = old_edge->blessed; // Is this correct?
  new_edge->invalid = old_edge->invalid;
  new_edge->quality = old_edge->quality;
  new_edge->grangered = ! old_edge->grangered;
  new_edge->reflected =   old_edge->reflected;
}


IntEdge_ID find_overlap_edge_mate
(/* Input Only */
 Tfragment frags[],
 Tedge edges[],
 IntEdge_ID ie0
 )
{
  // Find the other directed edge of the overlap and return the index.
  // Assert if the edge is not found.

  // This binary search version assumes that the edges are sorted with
  // respect to the compare_edge() function, and that all duplicate
  // edges with respect to the compare_edge() function have been
  // removed.

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

      Aedge *curr_edge = GetVA_Aedge(edges,ie2);

      icompare = compare_edge_function(curr_edge, &test_edge);

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



void fix_overlap_edge_mate
(/* Input Only */
 Tfragment frags[],
 Tedge edges[],
 IntEdge_ID ie0)
{
  IntEdge_ID ie1 = find_overlap_edge_mate( frags, edges, ie0);
  // The index of the other half of the undirected edge.
  Tnes ines0 = get_nes_edge(edges,ie0);
  Tnes ines1 = ines0; // This will be the "fix".

  if(AS_CGB_EDGE_NOT_FOUND == ie1) {
    IntFragment_ID avx = get_avx_edge(edges,ie0);
    IntFragment_ID bvx = get_bvx_edge(edges,ie0);
    IntFragment_ID aid = get_iid_fragment(frags,avx);
    IntFragment_ID bid = get_iid_fragment(frags,bvx);
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



static void fill_new_edge_with_reflected_old_edge(Tedge * edges,
                                                  const IntEdge_ID new_edge,
                                                  const IntEdge_ID old_edge) {
  Aedge * the_new_edge = GetVA_Aedge(edges,new_edge);
  Aedge * the_old_edge = GetVA_Aedge(edges,old_edge);
  reflect_Aedge( the_new_edge, the_old_edge);
}







void
verify_that_the_edges_are_in_order(Tedge edges[]) {
  IntEdge_ID ie0 = 0, nedge = GetNumEdges(edges);

  for (ie0=0; ie0 < nedge-1; ie0++) {
    int icompare = compare_edge_function(GetVA_Aedge(edges,ie0),GetVA_Aedge(edges,ie0+1));
    assert(icompare <= 0);
  }
}


void append_the_edge_mates(Tfragment frags[],
                           Tedge edges[]) {

  // Scan the current edges to find un-mated edges.  For each un-mated
  // edge, append the mate edge to the edge array.

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

  reorder_edges( frags, edges);
}


IntEdge_ID check_symmetry_of_the_edge_mates(Tfragment frags[],
                                            Tedge edges[]) {

  //  This was diagnostic code, used all over the place, but never did much.
  fprintf(stderr, "check_symmetry_of_the_edge_mates()--  Disabled.  (Was a NOP anyway)\n");
  return(0);

  IntEdge_ID
    ie0,
    nedge = GetNumEdges(edges),
    counter = 0;

  verify_that_the_edges_are_in_order(edges);

  for( ie0=0; ie0 < nedge; ie0++) {
    const IntEdge_ID ie1 = find_overlap_edge_mate( frags, edges, ie0);
    if( AS_CGB_EDGE_NOT_FOUND == ie1 )
      counter ++;
  }
  fprintf(stderr,"check_symmetry_of_the_edge_mates: nedge=" F_IID " counter=" F_IID "\n",
          nedge, counter);
  return counter;
}


//  BPW -- this really doesn't have anything to do with "edgemate" but
//  it's commonly called with check_symmetry_of_the_edge_mates(), and
//  they're both disabled.
//
void count_fragment_and_edge_labels(Tfragment frags[],
                                    Tedge     edges[],
                                    char      comment[]) {
  FILE *fout = stderr;

  fprintf(stderr, "count_fragment_and_edge_labels()--  Disabled.\n");
  return;

  const int nsample=500;
  const int nbucket=500;

  IntFragment_ID nfrag = GetNumFragments(frags);
  IntFragment_ID vid;
  Histogram_t *frag_lab_histogram = create_histogram(nsample,nbucket,TRUE,FALSE);

  fprintf(fout,"*** Histogram Fragment Labels <%s> ***\n",comment);

  for(vid=0; vid<nfrag; vid++) {
    const Tlab ilab = get_lab_fragment(frags,vid);
    add_to_histogram(frag_lab_histogram, (int)ilab, NULL);
  }

  fprintf(fout,"Histogram of the fragment label \n");
  print_histogram(fout,frag_lab_histogram, 0, 1);
  free_histogram(frag_lab_histogram);

  IntEdge_ID ie;
  IntEdge_ID nedge = GetNumEdges(edges);
  Histogram_t *inter_chunk_edge_nes_histogram = create_histogram(nsample,nbucket,TRUE,FALSE);
  Histogram_t *intra_chunk_edge_nes_histogram = create_histogram(nsample,nbucket,TRUE,FALSE);

  fprintf(fout,
          "*** Histogram Edge Labels (2 edges/overlap) <%s> ***\n",
          comment);

  for(ie=0; ie<nedge; ie++) {
    const Tnes nes = get_nes_edge(edges,ie);
    const IntFragment_ID avx = get_avx_edge(edges,ie);
    const IntFragment_ID bvx = get_bvx_edge(edges,ie);
    const IntChunk_ID a_cid = get_cid_fragment(frags,avx);
    const IntChunk_ID b_cid = get_cid_fragment(frags,bvx);
    if( a_cid == b_cid ) {
      add_to_histogram(intra_chunk_edge_nes_histogram, (int)nes, NULL);
    } else {
      add_to_histogram(inter_chunk_edge_nes_histogram, (int)nes, NULL);
    }
  }

  fprintf(fout,"Histogram of the inter-chunk overlap labels \n");
  print_histogram(fout,inter_chunk_edge_nes_histogram, 0, 1);
  free_histogram(inter_chunk_edge_nes_histogram);

  fprintf(fout,"Histogram of the intra-chunk overlap labels \n");
  print_histogram(fout,intra_chunk_edge_nes_histogram, 0, 1);
  free_histogram(intra_chunk_edge_nes_histogram);
}
