
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
= "$Id: AS_CGB_count_fragment_and_edge_labels.c,v 1.2 2004-09-23 20:25:01 mcschatz Exp $";
/* *******************************************************************
 *
 * Module: AS_CGB_count_fragment_and_edge_labels
 * 
 * Description: 
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
#define ZEBA

void count_fragment_and_edge_labels
(
 Tfragment frags[],
 Tedge     edges[],
 char      comment[])
{
  FILE *fout = stderr;

  {
    IntFragment_ID nfrag = GetNumFragments(frags);
    IntFragment_ID vid;
    const int nsample=500;
    const int nbucket=500;
    Histogram_t 
      *frag_lab_histogram 
      = create_histogram(nsample,nbucket,TRUE,FALSE);
    
    fprintf(fout,"*** Histogram Fragment Labels <%s> ***\n",comment);
      // #pragma omp parallel for
    for(vid=0; vid<nfrag; vid++) {
      const Tlab ilab = get_lab_fragment(frags,vid);
      add_to_histogram(frag_lab_histogram, (int)ilab, NULL);
    }

    fprintf(fout,"Histogram of the fragment label \n");
    print_histogram(fout,frag_lab_histogram, 0, 1);
    free_histogram(frag_lab_histogram);
  }
  
#ifndef ZEBA
  {
    IntEdge_ID ie;
    IntEdge_ID nedge = GetNumEdges(edges);
    const int nsample=500;
    const int nbucket=500;
    Histogram_t 
      *edge_nes_histogram 
      = create_histogram(nsample,nbucket,TRUE,FALSE);
    
    fprintf(fout,
	    "*** Histogram Edge Labels (2 edges/overlap) <%s> ***\n",
	    comment);
      // #pragma omp parallel for
    for(ie=0; ie<nedge; ie++) {
      const Tnes nes = get_nes_edge(edges,ie);
      add_to_histogram(edge_nes_histogram, (int)nes, NULL);
    }
    fprintf(fout,"Histogram of the overlap labels \n");
    print_histogram(fout,edge_nes_histogram, 0, 1);
    free_histogram(edge_nes_histogram);
  }
#else // ZEBA
  {
    IntEdge_ID ie;
    IntEdge_ID nedge = GetNumEdges(edges);
    const int nsample=500;
    const int nbucket=500;
    Histogram_t 
      *inter_chunk_edge_nes_histogram 
      = create_histogram(nsample,nbucket,TRUE,FALSE),
      *intra_chunk_edge_nes_histogram 
      = create_histogram(nsample,nbucket,TRUE,FALSE);
    
    fprintf(fout,
	    "*** Histogram Edge Labels (2 edges/overlap) <%s> ***\n",
	    comment);
      // #pragma omp parallel for
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
#endif // ZEBA
}
