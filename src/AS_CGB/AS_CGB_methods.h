
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

//  The methods to access the vertex and edge data store.

#ifndef AS_CGB_METHODS_INCLUDE
#define AS_CGB_METHODS_INCLUDE

#include "AS_OVS_overlapStore.h"

#define VAgetaccess(Type,va,index,member)  (((Type *)GetElement_VA(va,index))->member)

/* Object access functions */

#pragma inline copy_one_fragment
static void copy_one_fragment
(Tfragment destination[],IntFragment_ID idest,Tfragment source[],IntFragment_ID isrc)
/* This is a structure assignment. */
{ *GetAfragment(destination,idest) = *GetAfragment(source,isrc);}

#pragma inline copy_one_edge
static void copy_one_edge
(Tedge destination[],IntEdge_ID idest,Tedge source[],IntEdge_ID isrc)
/* This is a structure assignment. */
{ *GetAedge(destination,idest) = *GetAedge(source,isrc);}

#pragma inline GetNumEdges
static IntEdge_ID GetNumEdges(Tedge *edges) {
  return (IntEdge_ID) GetNumVA_Aedge(edges);
}

#pragma inline GetNumFragments
static IntFragment_ID GetNumFragments(Tfragment *frags) {
  return (IntFragment_ID) GetNumVA_Afragment(frags);
}

/* We need to be very explicit about primitive type casting due to
   "features" in the Digital UNIX C compiler. */

#pragma inline set_avx_edge
static void set_avx_edge
(Tedge edges[],const IntEdge_ID i,const IntFragment_ID value)
{ VAgetaccess(Aedge,edges,i,avx) = (IntFragment_ID)value;}
#pragma inline set_bvx_edge
static void set_bvx_edge
(Tedge edges[],const IntEdge_ID i,const IntFragment_ID value)
{ VAgetaccess(Aedge,edges,i,bvx) = (IntFragment_ID)value;}
#pragma inline set_asx_edge
static void set_asx_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,asx) = (int8)value;}
#pragma inline set_bsx_edge
static void set_bsx_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,bsx) = (int8)value;}
#pragma inline set_ahg_edge
static void set_ahg_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,ahg) = (int16)value;}
#pragma inline set_bhg_edge
static void set_bhg_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,bhg) = (int16)value;}
#pragma inline set_nes_edge
static void set_nes_edge(Tedge edges[],IntEdge_ID i,Tnes value)
{ VAgetaccess(Aedge,edges,i,nes) = (int8)value;}
#pragma inline set_inv_edge
static void set_inv_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,invalid) = value;}
#pragma inline set_reflected_edge
static void set_reflected_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,reflected) = value;}
#pragma inline set_grangered_edge
static void set_grangered_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,grangered) = value;}
#pragma inline set_qua_edge
static void set_qua_edge(const Tedge * const edges,IntEdge_ID i,uint32 value)
{ VAgetaccess(Aedge,edges,i,quality) = (uint32) value;}
#pragma inline set_blessed_edge
static void set_blessed_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,blessed) = value;}

#pragma inline get_avx_edge
static IntFragment_ID get_avx_edge(const Tedge * const edges,IntEdge_ID i)
{ return (IntFragment_ID) VAgetaccess(Aedge,edges,i,avx);}
#pragma inline get_bvx_edge
static IntFragment_ID get_bvx_edge(const Tedge * const edges,IntEdge_ID i)
{ return (IntFragment_ID) VAgetaccess(Aedge,edges,i,bvx);}
#pragma inline get_asx_edge
static int get_asx_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,asx);}
#pragma inline get_bsx_edge
static int get_bsx_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,bsx);}
#pragma inline get_ahg_edge
static int get_ahg_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,ahg);}
#pragma inline get_bhg_edge
static int get_bhg_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,bhg);}
#pragma inline get_nes_edge
static Tnes get_nes_edge(const Tedge * const edges,IntEdge_ID i)
{ return (Tnes) VAgetaccess(Aedge,edges,i,nes);}
#pragma inline get_inv_edge
static int get_inv_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,invalid);}
#pragma inline get_reflected_edge
static int get_reflected_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,reflected);}
#pragma inline get_grangered_edge
static int get_grangered_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,grangered);}
#pragma inline get_qua_edge
static uint32 get_qua_edge(const Tedge * const edges,IntEdge_ID i)
{ return (uint32) VAgetaccess(Aedge,edges,i,quality);}
#pragma inline get_blessed_edge
static int get_blessed_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,blessed);}

// Unitigger overlap classification assumes that overhangs of the
// overlap are non-negative (ahg>=0)||(bhg>=0)).  Note that the
// reverse of a dvt(dgn) overlap edge is a dvt(dgn) overlap edge
// whereas the reverse of a frc(toc) overlap edge is a toc(frc)
// overlap edge.
#pragma inline is_a_dvt_simple
static int is_a_dvt_simple(const int ahg, const int bhg) {
  // A dovetail overlap edge.
  assert( (ahg >= 0) || (bhg >= 0));
  return (ahg > 0) && (bhg > 0);
}
#pragma inline is_a_frc_simple
static int is_a_frc_simple(const int ahg, const int bhg) {
  // A from-the-contained-fragment containment overlap edge.
  assert( (ahg >= 0) || (bhg >= 0));
  return (ahg < 0) || ((ahg == 0)&&(bhg > 0));
}
#pragma inline is_a_toc_simple
static int is_a_toc_simple(const int ahg, const int bhg) {
  // A to-the-contained-fragment containment overlap edge.
  assert( (ahg >= 0) || (bhg >= 0));
  return (bhg < 0) || ((bhg == 0)&&(ahg > 0));
}
#pragma inline is_a_dgn_simple
static int is_a_dgn_simple(const int ahg, const int bhg) {
  // A degenerate containment overlap edge where the fragments are
  // mutually contained.
  assert( (ahg >= 0) || (bhg >= 0));
  return (((ahg == 0) && (bhg == 0))   &&
          (! is_a_dvt_simple(ahg,bhg)) &&
          (! is_a_frc_simple(ahg,bhg)) &&
          (! is_a_toc_simple(ahg,bhg)));
}

#pragma inline is_a_dvt_edge
static int is_a_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{
  const int ahg = VAgetaccess(Aedge,edges,i,ahg);
  const int bhg = VAgetaccess(Aedge,edges,i,bhg);
  return is_a_dvt_simple( ahg, bhg);
}
// A flag is returned that indicates if the overlap is a dovetail
// overlap.

#pragma inline is_a_frc_edge
#ifndef DEGENERATE_CONTAINMENT_RESOLUTION
static int is_a_frc_edge(const Tedge * const edges,IntEdge_ID i)
{
  const int ahg = VAgetaccess(Aedge,edges,i,ahg);
  const int bhg = VAgetaccess(Aedge,edges,i,bhg);
  return is_a_frc_simple( ahg, bhg);
}
#else // DEGENERATE CONTAINMENT RESOLUTION
static int get_frc_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int)
    (! is_a_dvt_edge(edges,i)) &&

    ((VAgetaccess(Aedge,edges,i,ahg) < 0) &&
     (VAgetaccess(Aedge,edges,i,bhg) >= 0) )
    ||
    ((VAgetaccess(Aedge,edges,i,ahg) == 0) &&
     (VAgetaccess(Aedge,edges,i,bhg) > 0) )
    ||
    ((VAgetaccess(Aedge,edges,i,ahg) == 0) &&
     (VAgetaccess(Aedge,edges,i,bhg) == 0) &&
     (VAgetaccess(Aedge,edges,i,avx) < VAgetaccess(Aedge,edges,i,bvx))
     // Note that this tie breaker that is sensitive to IID->VID remapping.
     )
    ;}
// A flag is returned that indicates if the overlap is a
// from-contained overlap.
#endif // DEGENERATE CONTAINMENT RESOLUTION

#pragma inline is_a_toc_edge
static int is_a_toc_edge(const Tedge * const edges,IntEdge_ID i)
{
  const int ahg = VAgetaccess(Aedge,edges,i,ahg);
  const int bhg = VAgetaccess(Aedge,edges,i,bhg);
  return is_a_toc_simple( ahg, bhg);
}

#pragma inline is_a_dgn_edge
static int is_a_dgn_edge(const Tedge * const edges,IntEdge_ID i)
{
  const int ahg = VAgetaccess(Aedge,edges,i,ahg);
  const int bhg = VAgetaccess(Aedge,edges,i,bhg);
  return is_a_dgn_simple( ahg, bhg);
}

#pragma inline get_intrachunk_dvt_edge
static int get_intrachunk_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{ return (AS_CGB_INTRACHUNK_EDGE == get_nes_edge(edges,i)) ; }

#pragma inline get_interchunk_dvt_edge
static int get_interchunk_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{ return
    get_intrachunk_dvt_edge(edges,i) ||
    (AS_CGB_INTERCHUNK_EDGE == get_nes_edge(edges,i));
}

#pragma inline get_thickest_dvt_edge
static int get_thickest_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{ return
    get_interchunk_dvt_edge(edges,i) ||
    (AS_CGB_THICKEST_EDGE == get_nes_edge(edges,i));
}


#pragma inline get_tied_dvt_edge
static int get_tied_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{
  return FALSE;
}


#pragma inline set_iid_fragment
static void set_iid_fragment(Tfragment frags[],IntFragment_ID i,IntFragment_ID value)
{ VAgetaccess(Afragment,frags,i,iid) = (IntFragment_ID)value;}
#pragma inline set_typ_fragment
static void set_typ_fragment(Tfragment frags[],IntFragment_ID i,FragType value)
{ VAgetaccess(Afragment,frags,i,frag_type) = (FragType)value;}
#pragma inline set_lab_fragment
static void set_lab_fragment(Tfragment frags[],IntFragment_ID i,Tlab value)
{ VAgetaccess(Afragment,frags,i,label) = value;}

#pragma inline set_o3p_fragment
static void set_o3p_fragment(Tfragment frags[],IntFragment_ID i,int64  value)
{ VAgetaccess(Afragment,frags,i,offset3p) = (int64 ) value;}
#pragma inline set_o5p_fragment
static void set_o5p_fragment(Tfragment frags[],IntFragment_ID i,int64  value)
{ VAgetaccess(Afragment,frags,i,offset5p) = (int64 ) value;}

#pragma inline set_length_fragment
static void set_length_fragment(Tfragment frags[],IntFragment_ID i,int32 value)
{ VAgetaccess(Afragment,frags,i,bp_length) = (int16)value;}

#pragma inline set_cid_fragment
static void set_cid_fragment(Tfragment frags[],IntFragment_ID i,IntChunk_ID value)
{ VAgetaccess(Afragment,frags,i,cid) = (IntChunk_ID)value;}

#pragma inline set_container_fragment
static void set_container_fragment(Tfragment frags[],IntFragment_ID i,
                                   IntFragment_ID value)
{ VAgetaccess(Afragment,frags,i,container) = (IntFragment_ID)value;}

#pragma inline set_del_fragment
static void set_del_fragment(Tfragment frags[],IntFragment_ID i,int value)
{ VAgetaccess(Afragment,frags,i,deleted) = value;}
#pragma inline set_con_fragment
static void set_con_fragment(Tfragment frags[],IntFragment_ID i,int value)
{ VAgetaccess(Afragment,frags,i,contained) = value;}
#pragma inline set_spur_fragment
static void set_spur_fragment(Tfragment frags[],IntFragment_ID i,int value)
{ VAgetaccess(Afragment,frags,i,spur) = value;}


#pragma inline get_iid_fragment
#pragma inline get_typ_fragment
#pragma inline get_lab_fragment
#pragma inline get_o3p_fragment
#pragma inline get_o5p_fragment
#pragma inline get_length_fragment
#pragma inline get_cid_fragment
#pragma inline get_container_fragment
#pragma inline get_del_fragment
#pragma inline get_con_fragment
#pragma inline get_spur_fragment
#pragma inline get_forward_fragment

static IntFragment_ID get_iid_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (IntFragment_ID) VAgetaccess(Afragment,frags,i,iid);}
static FragType get_typ_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (FragType) VAgetaccess(Afragment,frags,i,frag_type);}
static Tlab get_lab_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (Tlab) VAgetaccess(Afragment,frags,i,label);}

static int64  get_o3p_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (int64 ) VAgetaccess(Afragment,frags,i,offset3p);}
static int64  get_o5p_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (int64 ) VAgetaccess(Afragment,frags,i,offset5p);}

static int get_forward_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (VAgetaccess(Afragment,frags,i,offset3p) >
	  VAgetaccess(Afragment,frags,i,offset5p) );}
// Is this fragment in a forward orientation in its unitig?

static int32 get_length_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (int32) VAgetaccess(Afragment,frags,i,bp_length);}

static IntChunk_ID get_cid_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (IntChunk_ID) VAgetaccess(Afragment,frags,i,cid);}
static IntFragment_ID get_container_fragment
(const Tfragment * const frags,IntFragment_ID i)
{ return (IntFragment_ID) VAgetaccess(Afragment,frags,i,container);}

static int get_del_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (int) VAgetaccess(Afragment,frags,i,deleted);}
static int get_con_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (int) VAgetaccess(Afragment,frags,i,contained);}
static int get_spur_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (int) VAgetaccess(Afragment,frags,i,spur);}


#pragma inline set_segstart_vertex
static void set_segstart_vertex
(Tfragment frags[],IntFragment_ID i,int flag,IntEdge_ID value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,segbgn_suffix) = (IntEdge_ID)value;
  }else{
    VAgetaccess(Afragment,frags,i,segbgn_prefix) = (IntEdge_ID)value;
  }
}
#pragma inline set_segend_vertex
static void set_segend_vertex
(Tfragment frags[],IntFragment_ID i,int flag,IntEdge_ID value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,segend_suffix) = (IntEdge_ID)value;
  }else{
    VAgetaccess(Afragment,frags,i,segend_prefix) = (IntEdge_ID)value;
  }
}
#pragma inline set_seglen_vertex
static void set_seglen_vertex
(Tfragment frags[],IntFragment_ID i,int flag,int32 value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,nsuffix_all) = (int32)value;
  }else{
    VAgetaccess(Afragment,frags,i,nprefix_all) = (int32)value;
  }
}
#pragma inline set_seglen_frc_vertex
static void set_seglen_frc_vertex
(Tfragment frags[],IntFragment_ID i,int flag,int32 value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,nsuffix_frc) = (int32)value;
  }else{
    VAgetaccess(Afragment,frags,i,nprefix_frc) = (int32)value;
  }
}
#pragma inline set_seglen_dvt_vertex
static void set_seglen_dvt_vertex
(Tfragment frags[],IntFragment_ID i,int flag,int32 value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,nsuffix_dvt) = (int32)value;
  }else{
    VAgetaccess(Afragment,frags,i,nprefix_dvt) = (int32)value;
  }
}

#pragma inline get_segstart_vertex
static IntEdge_ID get_segstart_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{ return (IntEdge_ID) (flag ?
		       VAgetaccess(Afragment,frags,i,segbgn_suffix):
		       VAgetaccess(Afragment,frags,i,segbgn_prefix));}
#pragma inline get_segend_vertex
static IntEdge_ID get_segend_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{ return (IntEdge_ID) (flag ?
		       VAgetaccess(Afragment,frags,i,segend_suffix):
		       VAgetaccess(Afragment,frags,i,segend_prefix));}
#pragma inline get_seglen_vertex
static int32 get_seglen_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{ return (int32) (flag ?
		      VAgetaccess(Afragment,frags,i,nsuffix_all) :
		      VAgetaccess(Afragment,frags,i,nprefix_all));}
static int32 get_seglen_frc_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{ return (int32) (flag ?
		      VAgetaccess(Afragment,frags,i,nsuffix_frc) :
		      VAgetaccess(Afragment,frags,i,nprefix_frc));}
static int32 get_seglen_dvt_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{ return (int32) (flag ?
		      VAgetaccess(Afragment,frags,i,nsuffix_dvt) :
		      VAgetaccess(Afragment,frags,i,nprefix_dvt) );}

#pragma inline get_blessed_vertex
static int get_blessed_vertex
(const Tfragment * const frags, IntFragment_ID vid, int suffix)
{
  return (int) ( suffix ?
                 VAgetaccess(Afragment,frags,vid,suffix_blessed) :
                 VAgetaccess(Afragment,frags,vid,prefix_blessed) ); }


#pragma inline set_blessed_vertex
static void set_blessed_vertex
(const Tfragment * const frags, IntFragment_ID vid, int suffix, int value)
{
  if(suffix) {
    VAgetaccess(Afragment,frags,vid,suffix_blessed) = value;
  } else {
    VAgetaccess(Afragment,frags,vid,prefix_blessed) = value;
  }
}

#pragma inline get_best_ovl
static int32 get_best_ovl
(const Tfragment * const frags,
 const Tedge * const edges,
 IntEdge_ID ie)
{ IntFragment_ID iavx = get_avx_edge(edges,ie);
  int32  ilen = get_length_fragment(frags,iavx);
  return (int32) (ilen - get_ahg_edge(edges,ie));}


#pragma inline get_chunk_index
static IntChunk_ID get_chunk_index
(// Input only
 const TChunkMesg thechunks[],
 const Tfragment  frags[],
 const Tedge      edges[],
 const IntEdge_ID ie
 )
{
  /* This function returns whether the chunk index of the B-chunk in
     the overlap "ie". */

  const IntFragment_ID ibvx = get_bvx_edge(edges,ie); // Get the distal fragment
  // const int ibsx = get_bsx_edge(edges,ie); // and suffix flag
  const IntChunk_ID cbvx = get_cid_fragment(frags,ibvx); // then its chunk id.
  return cbvx;
}

#pragma inline get_chunk_suffix
static int get_chunk_suffix
(// Input only
 const TChunkMesg thechunks[],
 const Tfragment  frags[],
 const Tedge      edges[],
 const IntEdge_ID ie
 )
{
  /* This function returns whether the suffix of the B-chunk is in the
     overlap "ie". */

  const IntFragment_ID  ibvx = get_bvx_edge(edges,ie); // Get the distal fragment
  const int         ibsx = get_bsx_edge(edges,ie); // and suffix flag
  const int cbsx = (ibsx ^
		    (get_o3p_fragment(frags,ibvx) <
		     get_o5p_fragment(frags,ibvx) ));
  return cbsx;
}

#pragma inline set_raw_dvt_count_vertex
static void set_raw_dvt_count_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag,int value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,raw_suffix_dvt_count) = value;
  } else {
    VAgetaccess(Afragment,frags,i,raw_prefix_dvt_count) = value;
  }
}

#pragma inline set_raw_frc_count_fragment
static void set_raw_frc_count_fragment
(const Tfragment * const frags,IntFragment_ID i,int value)
{
  VAgetaccess(Afragment,frags,i,raw_frc_count) = value;
}

#pragma inline set_raw_toc_count_fragment
static void set_raw_toc_count_fragment
(const Tfragment * const frags,IntFragment_ID i,int value)
{
  VAgetaccess(Afragment,frags,i,raw_toc_count) = value;
}

#pragma inline get_raw_dvt_count_vertex
static int get_raw_dvt_count_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{
  int value;
  if(flag) {
    value = VAgetaccess(Afragment,frags,i,raw_suffix_dvt_count);
  } else {
    value = VAgetaccess(Afragment,frags,i,raw_prefix_dvt_count);
  }
  return value;
}

#pragma inline get_raw_frc_count_fragment
static int get_raw_frc_count_fragment
(const Tfragment * const frags,IntFragment_ID i)
{
  return VAgetaccess(Afragment,frags,i,raw_frc_count);
}

#pragma inline get_raw_toc_count_fragment
static int get_raw_toc_count_fragment
(const Tfragment * const frags,IntFragment_ID i)
{
  return VAgetaccess(Afragment,frags,i,raw_toc_count);
}

#pragma inline inc_raw_dvt_count_vertex
static void inc_raw_dvt_count_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,raw_suffix_dvt_count) ++;
  } else {
    VAgetaccess(Afragment,frags,i,raw_prefix_dvt_count) ++;
  }
}

#pragma inline inc_raw_frc_count_fragment
static void inc_raw_frc_count_fragment
(const Tfragment * const frags,IntFragment_ID i)
{
  VAgetaccess(Afragment,frags,i,raw_frc_count) ++;
}

#pragma inline inc_raw_toc_count_fragment
static void inc_raw_toc_count_fragment
(const Tfragment * const frags,IntFragment_ID i)
{
  VAgetaccess(Afragment,frags,i,raw_toc_count) ++;
}

#undef VAgetaccess


#endif /* AS_CGB_METHODS_INCLUDE */
