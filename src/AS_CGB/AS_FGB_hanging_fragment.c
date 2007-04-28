
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
#include "AS_CGB_all.h"
#include "AS_FGB_hanging_fragment.h"

#undef USE_REAPERS_DVT_DEGREE

void identify_early_spur_fragments
(
 Tfragment frags[],
 Tedge edges[])
{
  // The philosophy of early spur removal is to assume that a
  // hanging fragment should be ignored while the chunk backbones are
  // determined.  We must make exceptions when there is evidence to
  // consider the hanging fragment important.  In particular we need
  // to include hanging fragments at sequencing gaps in order to get a
  // tentative DNA sequence near the gap.  We also need to include a
  // hanging fragment as important if its dovetail overlap is used to
  // assign a fragment its "thru" status.  Also doubleton unitigs
  // should be assembled rather than output as two hanging fragments.

  // The chosen heuristic is if a hanging fragment "competes" in the
  // dovetail adjacency list of each neighboring fragment with a thru
  // fragment (in every location it can be placed), then aggresively
  // ignore the hanging fragment.  This is a generalization of the
  // original spur definition.  We do not assume that we know which
  // fragments are contained nor assume that we have a transitively
  // reduced overlap graph yet. The implementation uses a restricted
  // incoming dovetail degree computed from only thru fragments.

  // Determine the incoming dovetail degree at all fragment-ends
  // restricted so that only fragments-ends contribute to the degree.
  // This can be implemented by a scatter-add from the thru
  // fragment-ends using their outgoing dovetail adjacency list to the
  // edge's distal fragment-end.

  // The fragment re-labelling can be implemented by (1) first
  // re-assigning all hanging fragments to AS_CGB_HANGING_CRAPPY_FRAG,
  // and (2) for all active edges if the proximal fragment is hanging
  // (AS_CGB_HANGING_CRAPPY_FRAG or AS_CGB_HANGING_FRAG or ...) and
  // the distal fragment active (and non-contained?) and has its
  // restricted dovetail degree zero, then re-assign the proximal
  // fragment as AS_CGB_HANGING_FRAG.

  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);
  
  int * thru_fragment_prefix_dovetail_overlap_count = NULL;
  int * thru_fragment_suffix_dovetail_overlap_count = NULL;

  thru_fragment_prefix_dovetail_overlap_count = safe_malloc(sizeof(int) * nfrag);
  thru_fragment_suffix_dovetail_overlap_count = safe_malloc(sizeof(int) * nfrag);

  {
    // Initialize the incoming dovetail degree at each fragment.
    IntFragment_ID vid;
    // #pragma omp parallel for
    for(vid=0; vid<nfrag; vid++){
#if 1
      Tlab lab = get_lab_fragment(frags,vid);
      assert(
             (lab == AS_CGB_SOLO_FRAG) ||
             (lab == AS_CGB_HANGING_FRAG) ||
             (lab == AS_CGB_THRU_FRAG)
             );
#endif      
      thru_fragment_prefix_dovetail_overlap_count[vid] = 0;
      thru_fragment_suffix_dovetail_overlap_count[vid] = 0;
      set_spur_fragment(frags,vid,FALSE);
    }
  }

  {
    // Find the incoming dovetail degree from thru fragments at each
    // fragment by a scatter-add.
    IntEdge_ID ie;
    for(ie=0; ie<nedge; ie++) {
      const Tnes ines = get_nes_edge(edges,ie);
      switch(ines) {
      case AS_CGB_DOVETAIL_EDGE:
        {
          const IntFragment_ID avx = get_avx_edge(edges,ie);
          //const IntFragment_ID asx = get_asx_edge(edges,ie);
          const IntFragment_ID bvx = get_bvx_edge(edges,ie);
          const IntFragment_ID bsx = get_bsx_edge(edges,ie);
          //const IntFragment_ID aid = get_iid_fragment(frags,avx);
          //const IntFragment_ID bid = get_iid_fragment(frags,bvx);
          const Tlab alab = get_lab_fragment(frags,avx);
          //const Tlab blab = get_lab_fragment(frags,bvx);
          
          if( alab == AS_CGB_THRU_FRAG ) {
            if( bsx ) {
              thru_fragment_suffix_dovetail_overlap_count[bvx] ++;
            } else {
              thru_fragment_prefix_dovetail_overlap_count[bvx] ++;
            }
          }
        }
        break;

      case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
      case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
      case AS_CGB_CONTAINED_EDGE:
      case AS_CGB_REMOVED_BY_DUPLICATE_CON:
      case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:

        // do nothing ...
        break;
      default:
        fprintf(stderr,"Unexpected overlap edge label %d\n",ines);
        assert(FALSE);
      }
    }
  }

  {
    // All the hanging fragments will be removed unless redeemed.
    IntFragment_ID vid;
    // #pragma omp parallel for
    for(vid=0; vid<nfrag; vid++){
      const Tlab lab = get_lab_fragment(frags,vid);
      switch(lab) {
      case AS_CGB_HANGING_FRAG:
      case AS_CGB_HANGING_CRAPPY_FRAG:
      case AS_CGB_HANGING_CHUNK_FRAG:
        set_lab_fragment(frags,vid,AS_CGB_HANGING_CRAPPY_FRAG);
        set_spur_fragment(frags,vid,TRUE);
        break;
      case AS_CGB_SOLO_FRAG:
      case AS_CGB_THRU_FRAG:
      case AS_CGB_UNPLACEDCONT_FRAG:  // work-around for bubble smoothing.
        // do nothing ...
        break;
      default:
        fprintf(stderr,"Unexpected fragment label %d\n",lab);
        assert(FALSE);
      }
    }
  }

  {
    // We consider each hanging fragment individually to determine if
    // it needs to be retained for chunking.  Check every (distal)
    // fragment-end that is dovetail overlap connected with the given
    // hanging fragment.  If any distal fragment-end has zero
    // incoming-dovetail-degree-from-thru-fragments, then the hanging
    // fragment can not be labeled as a spur fragment.

    IntEdge_ID ie;
    for(ie=0; ie<nedge; ie++) {
      const Tnes ines = get_nes_edge(edges,ie);
      switch(ines) {
      case AS_CGB_DOVETAIL_EDGE:
        {
          const IntFragment_ID avx = get_avx_edge(edges,ie);
          //const IntFragment_ID asx = get_asx_edge(edges,ie);
          const IntFragment_ID bvx = get_bvx_edge(edges,ie);
          const IntFragment_ID bsx = get_bsx_edge(edges,ie);
          //const IntFragment_ID aid = get_iid_fragment(frags,avx);
          //const IntFragment_ID bid = get_iid_fragment(frags,bvx);
          const Tlab alab = get_lab_fragment(frags,avx);
          //const Tlab blab = get_lab_fragment(frags,bvx);
          if(
             (get_spur_fragment(frags,avx) == TRUE) ||
             (alab == AS_CGB_HANGING_CRAPPY_FRAG)
             ) {
            const int distal_incoming_thru_fragment_dovetail_degree = 
              ( bsx
                ? thru_fragment_suffix_dovetail_overlap_count[bvx]
                : thru_fragment_prefix_dovetail_overlap_count[bvx]);
            if(0 == distal_incoming_thru_fragment_dovetail_degree) {
	      set_lab_fragment(frags,avx,AS_CGB_HANGING_FRAG);
              set_spur_fragment(frags,avx,FALSE);
	    }
          }
        }
        break;

      case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
      case AS_CGB_CONTAINED_EDGE:
      case AS_CGB_REMOVED_BY_DUPLICATE_CON:
        // do nothing ...
        break;
      default:
        fprintf(stderr,"Unexpected overlap edge label %d\n",ines);
        assert(FALSE);
      }
    }
  }
  fprintf(stderr, "Classification of early spur hanging fragments done.\n");
  safe_free(thru_fragment_prefix_dovetail_overlap_count);
  safe_free(thru_fragment_suffix_dovetail_overlap_count);
}


void separate_fragments_as_solo_hanging_thru
(
 Tfragment frags[],
 Tedge edges[])
{
  // A "solo" fragment is has a zero outgoing dovetail overlap
  // edge degree on both fragment-ends.

  // A "hanging" fragment is has a zero outgoing dovetail overlap
  // edge degree for just one of the fragment-ends.

  // A "thru" fragment is has a positive outgoing dovetail overlap
  // edge degree on both fragment-ends.

  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);

#ifndef USE_REAPERS_DVT_DEGREE
  int * all_fragment_prefix_dovetail_overlap_count = NULL;
  int * all_fragment_suffix_dovetail_overlap_count = NULL;

  all_fragment_prefix_dovetail_overlap_count = safe_malloc(sizeof(int) * nfrag);
  all_fragment_suffix_dovetail_overlap_count = safe_malloc(sizeof(int) * nfrag);

  {
    // Initialize the outgoing dovetail degree at each fragment.
    IntFragment_ID vid;
    // #pragma omp parallel for
    for(vid=0; vid<nfrag; vid++){
      all_fragment_prefix_dovetail_overlap_count[vid] = 0;
      all_fragment_suffix_dovetail_overlap_count[vid] = 0;
    }
  }

  {
    // Find the outgoing dovetail degree at each fragment by a
    // scatter-add.
    IntEdge_ID ie;
    for(ie=0; ie<nedge; ie++) {
      const Tnes ines = get_nes_edge(edges,ie);
      switch(ines) {
      case AS_CGB_DOVETAIL_EDGE:
        {
          const IntFragment_ID avx = get_avx_edge(edges,ie);
          const IntFragment_ID asx = get_asx_edge(edges,ie);
          const IntFragment_ID bvx = get_bvx_edge(edges,ie);
          //const IntFragment_ID bsx = get_bsx_edge(edges,ie);
          //const IntFragment_ID aid = get_iid_fragment(frags,avx);
          //const IntFragment_ID bid = get_iid_fragment(frags,bvx);
          const Tlab alab = get_lab_fragment(frags,avx);
          const Tlab blab = get_lab_fragment(frags,bvx);
          
          if( (alab != AS_CGB_DELETED_FRAG) &&
              (blab != AS_CGB_DELETED_FRAG) &&
              (alab != AS_CGB_REMOVED_BREAKER_FRAG) &&
              (blab != AS_CGB_REMOVED_BREAKER_FRAG) 
              ) {
            if( asx ) {
              all_fragment_suffix_dovetail_overlap_count[avx] ++;
            } else {
              all_fragment_prefix_dovetail_overlap_count[avx] ++;
            }
          }
        }
        break;
      case AS_CGB_CONTAINED_EDGE:
        // do nothing ...
        break;
      default:
        fprintf(stderr,"Unexpected overlap edge label %d\n",ines);
        assert(FALSE);
      }
    }
  }
#endif // USE_REAPERS_DVT_DEGREE

  {
    const IntFragment_ID nfrag = GetNumFragments(frags);
    IntFragment_ID solo_fragment_count = 0;
    IntFragment_ID prefix_hanging_fragment_count = 0;
    IntFragment_ID suffix_hanging_fragment_count = 0;
    IntFragment_ID thru_fragment_count = 0;
    
    IntFragment_ID vid;
    // #pragma omp parallel for

    for(vid=0; vid<nfrag; vid++){
      //const int con = get_con_fragment(frags,vid);
      const Tlab lab = get_lab_fragment(frags,vid);

      const int raw_fragment_prefix_dovetail_overlap_count
        = get_raw_dvt_count_vertex( frags, vid, FALSE);
      const int raw_fragment_suffix_dovetail_overlap_count
        = get_raw_dvt_count_vertex( frags, vid, TRUE);
      
      const int cur_fragment_prefix_dovetail_overlap_count
#ifndef USE_REAPERS_DVT_DEGREE
        = all_fragment_prefix_dovetail_overlap_count[vid];
#else // USE_REAPERS_DVT_DEGREE
        = get_seglen_dvt_vertex(frags,vid,FALSE);
#endif // USE_REAPERS_DVT_DEGREE


      const int cur_fragment_suffix_dovetail_overlap_count
#ifndef USE_REAPERS_DVT_DEGREE
        = all_fragment_suffix_dovetail_overlap_count[vid];
#else // USE_REAPERS_DVT_DEGREE
        = get_seglen_dvt_vertex(frags,vid,TRUE);
#endif // USE_REAPERS_DVT_DEGREE

      const int fragment_prefix_dovetail_overlap_count
#if 1
        = cur_fragment_prefix_dovetail_overlap_count;
#else
        = raw_fragment_prefix_dovetail_overlap_count;
#endif
      
      const int fragment_suffix_dovetail_overlap_count
#if 1        
      = cur_fragment_suffix_dovetail_overlap_count;
#else
      = raw_fragment_suffix_dovetail_overlap_count;
#endif
      
      Tlab ilab = lab;

      assert( cur_fragment_prefix_dovetail_overlap_count <=
              raw_fragment_prefix_dovetail_overlap_count );
      assert( cur_fragment_suffix_dovetail_overlap_count <=
              raw_fragment_suffix_dovetail_overlap_count );

      {
	IntFragment_ID iid = get_iid_fragment(frags,vid);
        if(!(! (raw_fragment_prefix_dovetail_overlap_count>0)
             ||(cur_fragment_prefix_dovetail_overlap_count>0))
           ||
           !(! (raw_fragment_suffix_dovetail_overlap_count>0)
             ||(cur_fragment_suffix_dovetail_overlap_count>0))
           ) {
          fprintf(stdout,"GFC "
		  "iid=" F_IID " "
		  "raw_fragment_prefix_dovetail_overlap_count=%d "
		  "raw_fragment_suffix_dovetail_overlap_count=%d "
		  "cur_fragment_prefix_dovetail_overlap_count=%d "
		  "cur_fragment_suffix_dovetail_overlap_count=%d\n",
		  iid,
		  raw_fragment_prefix_dovetail_overlap_count,
		  raw_fragment_suffix_dovetail_overlap_count,
		  cur_fragment_prefix_dovetail_overlap_count,
		  cur_fragment_suffix_dovetail_overlap_count);
	}
      }

      if(
         (ilab != AS_CGB_DELETED_FRAG) &&
         (ilab != AS_CGB_REMOVED_BREAKER_FRAG)
         ) {

        if((fragment_prefix_dovetail_overlap_count == 0) &&
           (fragment_suffix_dovetail_overlap_count == 0)) {
          ilab = AS_CGB_SOLO_FRAG;
          solo_fragment_count ++;
        } else if ((fragment_prefix_dovetail_overlap_count == 0) ||
                   (fragment_suffix_dovetail_overlap_count == 0)) {
          ilab = AS_CGB_HANGING_FRAG;
          if (fragment_prefix_dovetail_overlap_count == 0) {
            prefix_hanging_fragment_count ++;
          }
          if (fragment_suffix_dovetail_overlap_count == 0) {
            suffix_hanging_fragment_count ++;
          }
        } else {
          ilab = AS_CGB_THRU_FRAG;
          thru_fragment_count ++;
        }
      }

      if( ilab != lab ) {
        set_lab_fragment(frags,vid,ilab);
      }
    }

    fprintf(stderr, "separate_fragments_as_solo_hanging_thru::\n");
    fprintf(stderr, "          solo_fragment_count = " F_IID "\n", solo_fragment_count);
    fprintf(stderr, "prefix_hanging_fragment_count = " F_IID "\n", prefix_hanging_fragment_count);
    fprintf(stderr, "suffix_hanging_fragment_count = " F_IID "\n", suffix_hanging_fragment_count);
    fprintf(stderr, "          thru_fragment_count = " F_IID "\n", thru_fragment_count);

  }
#ifndef USE_REAPERS_DVT_DEGREE
  safe_free(all_fragment_prefix_dovetail_overlap_count);
  safe_free(all_fragment_suffix_dovetail_overlap_count);
#endif // USE_REAPERS_DVT_DEGREE
}
