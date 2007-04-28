
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
#include "AS_FGB_buildFragmentHash.h"


FragmentHashObject * build_FragmentHash
(
 const Tfragment * const frags,
 const IntFragment_ID as_cgb_max_frag_iid
)
{
  FragmentHashObject * afr_to_avx = NULL;
  IntFragment_ID max_frag_iid = 0;
  IntFragment_ID min_frag_iid = ~((IntFragment_ID)0);

  /* The current mapping from IIDs to fragment index. */
  if( frags != NULL ) {
    const IntFragment_ID nfrag = GetNumFragments(frags); 
    if(nfrag > 0){
      IntFragment_ID max_iid = get_iid_fragment(frags,0);
      IntFragment_ID min_iid = get_iid_fragment(frags,0);
      IntFragment_ID iv;
      for(iv=1;iv<nfrag;iv++) {
        const IntFragment_ID iid = get_iid_fragment(frags,iv);
        /* reset afr_to_avx[] */
        max_iid = MAX(max_iid,iid);
        min_iid = MIN(min_iid,iid);
      }
      max_frag_iid = max_iid;
      min_frag_iid = min_iid;
    }
  }

  /* Set the translation mappings converting array location and
     assembler internal fragment ids.  */

  max_frag_iid = MAX(as_cgb_max_frag_iid, max_frag_iid);
  afr_to_avx   = create_FragmentHash(max_frag_iid+1);

  /* The current mapping from IIDs to fragment index. */
  if( frags != NULL ) {
    const IntFragment_ID nfrag = GetNumFragments(frags); 
    if(nfrag > 0){
      IntFragment_ID iv;
      for(iv=0;iv<nfrag;iv++) {
        const IntFragment_ID iid = get_iid_fragment(frags,iv);
        /* reset afr_to_avx[] */
        assert(iid <= max_frag_iid);
        set_vid_FragmentHash(afr_to_avx,iid,iv);
      }
    }
  }
  return afr_to_avx;
}

