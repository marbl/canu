
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
= "$Id: AS_FGB_contained.c,v 1.2 2004-09-23 20:25:01 mcschatz Exp $";
/*********************************************************************
 *
 * Module: AS_FGB_contained.c
 * Description: Determines the contained status of the fragments.
 * Assumptions:
 * 
 * Author: Clark Mobarry
 *********************************************************************/

/*************************************************************************/
/* System include files */

/*************************************************************************/
/* Local include files */
#include "AS_CGB_all.h"

/* Conditional compilation */

/*************************************************************************/
/* Constant definitions; Macro definitions; type definitions */

/*************************************************************************/
/* Static Globals */

/*************************************************************************/
/* Function prototypes for internal static functions */

/**************************************************************/
/* Utility functions */

void check_containment_edges
(
 Tfragment frags[],
 Tedge edges[])
{
  /* 
   */
  const IntEdge_ID nedge = GetNumEdges(edges);

  { IntEdge_ID ie; for(ie=0; ie<nedge; ie++) {
    const IntFragment_ID iavx = get_avx_edge(edges,ie);
    const IntFragment_ID ibvx = get_bvx_edge(edges,ie);
    const IntFragment_ID aid = get_iid_fragment(frags,iavx);
    const IntFragment_ID bid = get_iid_fragment(frags,ibvx);
    const int iahg = get_ahg_edge(edges,ie);
    const int ibhg = get_bhg_edge(edges,ie);
    const Tnes ines = get_nes_edge(edges,ie);

    if(!(/* The containment overlaps... */
         (AS_CGB_CONTAINED_EDGE == ines) ||
         /* The dovetail overlaps... */
         (AS_CGB_DOVETAIL_EDGE == ines)
         )
       ) {
      fprintf(stderr,"ie=" F_IID " nes=%d\n", ie,ines);
    }
    assert(/* The containment overlaps... */
	   (AS_CGB_CONTAINED_EDGE == ines) ||
	   /* The dovetail overlaps... */
           (AS_CGB_DOVETAIL_EDGE == ines) ||
	   /* The other overlaps ... */
	   (AS_CGB_REMOVED_BY_DUPLICATE_DVT == ines) ||
	   (AS_CGB_REMOVED_BY_DUPLICATE_CON == ines));
    
    if((AS_CGB_REMOVED_BY_DUPLICATE_DVT != ines) &&
       (AS_CGB_REMOVED_BY_DUPLICATE_CON != ines)
       ) {
      assert(!((iahg<0)&&(ibhg<0)));
      assert((!(iahg<0))||(ibhg>=0)); /* if A then B  is (!A)||(B). */
      assert((!(ibhg<0))||(iahg>=0)); /* if A then B  is (!A)||(B). */
      
      if(!( is_a_dvt_simple(iahg,ibhg) ||
            ((is_a_toc_simple(iahg,ibhg) || is_a_frc_simple(iahg,ibhg))
             && (ines == AS_CGB_CONTAINED_EDGE))
           )) {
	fprintf(stderr,"XX " F_IID " " F_IID " " F_IID " " F_IID " " F_IID " %d %d %d\n", 
		ie, aid, bid,
		iavx, ibvx, iahg, ibhg, ines);
	fprintf(stderr,"Break the containment tie!!\n");
	assert(FALSE);
      }
    }
  }}
}

void contained_fragment_marking_frc
(
 Tfragment frags[],
 Tedge edges[])
{
  const IntFragment_ID nfrag = GetNumEdges(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);

  { IntFragment_ID iv; for(iv=0; iv<nfrag; iv++) {
    set_con_fragment(frags,iv,FALSE);
  }}
  
  { IntEdge_ID ie; for(ie=0; ie<nedge; ie++) {
    const IntFragment_ID iavx = get_avx_edge(edges,ie);
    const IntFragment_ID ibvx = get_bvx_edge(edges,ie);
#if 0
    const IntFragment_ID aid = get_iid_fragment(frags,iavx);
    const IntFragment_ID bid = get_iid_fragment(frags,ibvx);
    const Tlab alabel = get_lab_fragment(frags,iavx);
    const Tlab blabel = get_lab_fragment(frags,ibvx);
#endif
    //const int iahg = get_ahg_edge(edges,ie);
    //const int ibhg = get_bhg_edge(edges,ie);
    const Tnes ines = get_nes_edge(edges,ie);

    switch(ines) {
    case AS_CGB_CONTAINED_EDGE:
      if(
         (is_a_frc_edge(edges,ie) ||
          (is_a_dgn_edge(edges,ie)&&(get_avx_edge(edges,ie) > get_bvx_edge(edges,ie))))
         &&
         (!(AS_CGB_HANGING_CRAPPY_FRAG == get_lab_fragment(frags,ibvx))) &&
         (FALSE == get_spur_fragment(frags,ibvx)) &&
         (! get_del_fragment(frags,iavx)) &&
         (! get_del_fragment(frags,ibvx))
         ) {
        // All containment relationships touching spur fragments
        // are ignored.
        set_con_fragment(frags,iavx,TRUE);
        set_lab_fragment(frags,iavx,AS_CGB_UNPLACEDCONT_FRAG);
      }
      break;
    case AS_CGB_DOVETAIL_EDGE:
    case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
    case AS_CGB_REMOVED_BY_DUPLICATE_CON:
      // do nothing ...
      break;
    default:
      fprintf(stderr,"Unexpected overlap edge type %d\n", ines);
      assert(FALSE);
      break;
    }
  }}
}
