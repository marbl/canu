
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

static const char *rcsid_AS_CGB_METHODS_INCLUDE = "$Id$";

#include "AS_OVS_overlapStore.H"

#define VAgetaccess(Type,va,index,member)  (((Type *)GetElement_VA(va,index))->member)

/* Object access functions */

static void copy_one_fragment
(Tfragment destination[],IntFragment_ID idest,Tfragment source[],IntFragment_ID isrc)
/* This is a structure assignment. */
{ *GetAfragment(destination,idest) = *GetAfragment(source,isrc);}

static void copy_one_edge
(Tedge destination[],IntEdge_ID idest,Tedge source[],IntEdge_ID isrc)
/* This is a structure assignment. */
{ *GetAedge(destination,idest) = *GetAedge(source,isrc);}

static IntEdge_ID GetNumEdges(Tedge *edges) {
  return (IntEdge_ID) GetNumVA_Aedge(edges);
}

static IntFragment_ID GetNumFragments(Tfragment *frags) {
  return (IntFragment_ID) GetNumVA_Afragment(frags);
}

/* We need to be very explicit about primitive type casting due to
   "features" in the Digital UNIX C compiler. */

static void set_avx_edge
(Tedge edges[],const IntEdge_ID i,const IntFragment_ID value)
{ VAgetaccess(Aedge,edges,i,avx) = (IntFragment_ID)value;}
static void set_bvx_edge
(Tedge edges[],const IntEdge_ID i,const IntFragment_ID value)
{ VAgetaccess(Aedge,edges,i,bvx) = (IntFragment_ID)value;}
static void set_asx_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,asx) = (int8)value;}
static void set_bsx_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,bsx) = (int8)value;}
static void set_ahg_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,ahg) = (int32)value;}
static void set_bhg_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,bhg) = (int32)value;}
static void set_nes_edge(Tedge edges[],IntEdge_ID i,Tnes value)
{ VAgetaccess(Aedge,edges,i,nes) = (int8)value;}
static void set_inv_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,invalid) = value;}
static void set_reflected_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,reflected) = value;}
static void set_grangered_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,grangered) = value;}
static void set_qua_edge(const Tedge * const edges,IntEdge_ID i,uint32 value)
{ VAgetaccess(Aedge,edges,i,quality) = (uint32) value;}
static void set_blessed_edge(Tedge edges[],IntEdge_ID i,int value)
{ VAgetaccess(Aedge,edges,i,blessed) = value;}

static IntFragment_ID get_avx_edge(const Tedge * const edges,IntEdge_ID i)
{ return (IntFragment_ID) VAgetaccess(Aedge,edges,i,avx);}
static IntFragment_ID get_bvx_edge(const Tedge * const edges,IntEdge_ID i)
{ return (IntFragment_ID) VAgetaccess(Aedge,edges,i,bvx);}
static int get_asx_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,asx);}
static int get_bsx_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,bsx);}
static int get_ahg_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,ahg);}
static int get_bhg_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,bhg);}
static Tnes get_nes_edge(const Tedge * const edges,IntEdge_ID i)
{ return (Tnes) VAgetaccess(Aedge,edges,i,nes);}
static int get_inv_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,invalid);}
static int get_reflected_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,reflected);}
static int get_grangered_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,grangered);}
static uint32 get_qua_edge(const Tedge * const edges,IntEdge_ID i)
{ return (uint32) VAgetaccess(Aedge,edges,i,quality);}
static int get_blessed_edge(const Tedge * const edges,IntEdge_ID i)
{ return (int) VAgetaccess(Aedge,edges,i,blessed);}

// Unitigger overlap classification assumes that overhangs of the
// overlap are non-negative (ahg>=0)||(bhg>=0)).  Note that the
// reverse of a dvt(dgn) overlap edge is a dvt(dgn) overlap edge
// whereas the reverse of a frc(toc) overlap edge is a toc(frc)
// overlap edge.
static int is_a_dvt_simple(const int ahg, const int bhg) {
  // A dovetail overlap edge.
  assert( (ahg >= 0) || (bhg >= 0));
  return (ahg > 0) && (bhg > 0);
}
static int is_a_frc_simple(const int ahg, const int bhg) {
  // A from-the-contained-fragment containment overlap edge.
  assert( (ahg >= 0) || (bhg >= 0));
  return (ahg < 0) || ((ahg == 0)&&(bhg > 0));
}
static int is_a_toc_simple(const int ahg, const int bhg) {
  // A to-the-contained-fragment containment overlap edge.
  assert( (ahg >= 0) || (bhg >= 0));
  return (bhg < 0) || ((bhg == 0)&&(ahg > 0));
}
static int is_a_dgn_simple(const int ahg, const int bhg) {
  // A degenerate containment overlap edge where the fragments are
  // mutually contained.
  assert( (ahg >= 0) || (bhg >= 0));
  return (((ahg == 0) && (bhg == 0))   &&
          (! is_a_dvt_simple(ahg,bhg)) &&
          (! is_a_frc_simple(ahg,bhg)) &&
          (! is_a_toc_simple(ahg,bhg)));
}

static int is_a_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{
  const int ahg = VAgetaccess(Aedge,edges,i,ahg);
  const int bhg = VAgetaccess(Aedge,edges,i,bhg);
  return is_a_dvt_simple( ahg, bhg);
}
// A flag is returned that indicates if the overlap is a dovetail
// overlap.

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

static int is_a_toc_edge(const Tedge * const edges,IntEdge_ID i)
{
  const int ahg = VAgetaccess(Aedge,edges,i,ahg);
  const int bhg = VAgetaccess(Aedge,edges,i,bhg);
  return is_a_toc_simple( ahg, bhg);
}

static int is_a_dgn_edge(const Tedge * const edges,IntEdge_ID i)
{
  const int ahg = VAgetaccess(Aedge,edges,i,ahg);
  const int bhg = VAgetaccess(Aedge,edges,i,bhg);
  return is_a_dgn_simple( ahg, bhg);
}

static int get_intrachunk_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{ return (AS_CGB_INTRACHUNK_EDGE == get_nes_edge(edges,i)) ; }

static int get_interchunk_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{ return
    get_intrachunk_dvt_edge(edges,i) ||
    (AS_CGB_INTERCHUNK_EDGE == get_nes_edge(edges,i));
}

static int get_thickest_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{ return
    get_interchunk_dvt_edge(edges,i) ||
    (AS_CGB_THICKEST_EDGE == get_nes_edge(edges,i));
}


static int get_tied_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{
  return FALSE;
}


static void set_iid_fragment(Tfragment frags[],IntFragment_ID i,IntFragment_ID value)
{ VAgetaccess(Afragment,frags,i,iid) = (IntFragment_ID)value;}
static void set_typ_fragment(Tfragment frags[],IntFragment_ID i,FragType value)
{ VAgetaccess(Afragment,frags,i,frag_type) = (FragType)value;}
static void set_lab_fragment(Tfragment frags[],IntFragment_ID i,Tlab value)
{ VAgetaccess(Afragment,frags,i,label) = value;}

static void set_o3p_fragment(Tfragment frags[],IntFragment_ID i,int64  value)
{ VAgetaccess(Afragment,frags,i,offset3p) = (int64 ) value;}
static void set_o5p_fragment(Tfragment frags[],IntFragment_ID i,int64  value)
{ VAgetaccess(Afragment,frags,i,offset5p) = (int64 ) value;}

static void set_length_fragment(Tfragment frags[],IntFragment_ID i,int32 value)
{ VAgetaccess(Afragment,frags,i,bp_length) = (int32)value;}

static void set_cid_fragment(Tfragment frags[],IntFragment_ID i,IntChunk_ID value)
{ VAgetaccess(Afragment,frags,i,cid) = (IntChunk_ID)value;}

static void set_container_fragment(Tfragment frags[],IntFragment_ID i,
                                   IntFragment_ID value)
{ VAgetaccess(Afragment,frags,i,container) = (IntFragment_ID)value;}

static void set_del_fragment(Tfragment frags[],IntFragment_ID i,int value)
{ VAgetaccess(Afragment,frags,i,deleted) = value;}
static void set_con_fragment(Tfragment frags[],IntFragment_ID i,int value)
{ VAgetaccess(Afragment,frags,i,contained) = value;}
static void set_spur_fragment(Tfragment frags[],IntFragment_ID i,int value)
{ VAgetaccess(Afragment,frags,i,spur) = value;}



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


static void set_segstart_vertex
(Tfragment frags[],IntFragment_ID i,int flag,IntEdge_ID value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,segbgn_suffix) = (IntEdge_ID)value;
  }else{
    VAgetaccess(Afragment,frags,i,segbgn_prefix) = (IntEdge_ID)value;
  }
}
static void set_segend_vertex
(Tfragment frags[],IntFragment_ID i,int flag,IntEdge_ID value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,segend_suffix) = (IntEdge_ID)value;
  }else{
    VAgetaccess(Afragment,frags,i,segend_prefix) = (IntEdge_ID)value;
  }
}
static void set_seglen_vertex
(Tfragment frags[],IntFragment_ID i,int flag,int32 value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,nsuffix_all) = (int32)value;
  }else{
    VAgetaccess(Afragment,frags,i,nprefix_all) = (int32)value;
  }
}
static void set_seglen_frc_vertex
(Tfragment frags[],IntFragment_ID i,int flag,int32 value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,nsuffix_frc) = (int32)value;
  }else{
    VAgetaccess(Afragment,frags,i,nprefix_frc) = (int32)value;
  }
}
static void set_seglen_dvt_vertex
(Tfragment frags[],IntFragment_ID i,int flag,int32 value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,nsuffix_dvt) = (int32)value;
  }else{
    VAgetaccess(Afragment,frags,i,nprefix_dvt) = (int32)value;
  }
}

static IntEdge_ID get_segstart_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{ return (IntEdge_ID) (flag ?
		       VAgetaccess(Afragment,frags,i,segbgn_suffix):
		       VAgetaccess(Afragment,frags,i,segbgn_prefix));}
static IntEdge_ID get_segend_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{ return (IntEdge_ID) (flag ?
		       VAgetaccess(Afragment,frags,i,segend_suffix):
		       VAgetaccess(Afragment,frags,i,segend_prefix));}
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

static int get_blessed_vertex
(const Tfragment * const frags, IntFragment_ID vid, int suffix)
{
  return (int) ( suffix ?
                 VAgetaccess(Afragment,frags,vid,suffix_blessed) :
                 VAgetaccess(Afragment,frags,vid,prefix_blessed) ); }


static void set_blessed_vertex
(const Tfragment * const frags, IntFragment_ID vid, int suffix, int value)
{
  if(suffix) {
    VAgetaccess(Afragment,frags,vid,suffix_blessed) = value;
  } else {
    VAgetaccess(Afragment,frags,vid,prefix_blessed) = value;
  }
}

static int32 get_best_ovl
(const Tfragment * const frags,
 const Tedge * const edges,
 IntEdge_ID ie)
{ IntFragment_ID iavx = get_avx_edge(edges,ie);
  int32  ilen = get_length_fragment(frags,iavx);
  return (int32) (ilen - get_ahg_edge(edges,ie));}


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

static void set_raw_dvt_count_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag,int value)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,raw_suffix_dvt_count) = value;
  } else {
    VAgetaccess(Afragment,frags,i,raw_prefix_dvt_count) = value;
  }
}

static void set_raw_frc_count_fragment
(const Tfragment * const frags,IntFragment_ID i,int value)
{
  VAgetaccess(Afragment,frags,i,raw_frc_count) = value;
}

static void set_raw_toc_count_fragment
(const Tfragment * const frags,IntFragment_ID i,int value)
{
  VAgetaccess(Afragment,frags,i,raw_toc_count) = value;
}

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

static int get_raw_frc_count_fragment
(const Tfragment * const frags,IntFragment_ID i)
{
  return VAgetaccess(Afragment,frags,i,raw_frc_count);
}

static int get_raw_toc_count_fragment
(const Tfragment * const frags,IntFragment_ID i)
{
  return VAgetaccess(Afragment,frags,i,raw_toc_count);
}

static void inc_raw_dvt_count_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{
  if(flag) {
    VAgetaccess(Afragment,frags,i,raw_suffix_dvt_count) ++;
  } else {
    VAgetaccess(Afragment,frags,i,raw_prefix_dvt_count) ++;
  }
}

static void inc_raw_frc_count_fragment
(const Tfragment * const frags,IntFragment_ID i)
{
  VAgetaccess(Afragment,frags,i,raw_frc_count) ++;
}

static void inc_raw_toc_count_fragment
(const Tfragment * const frags,IntFragment_ID i)
{
  VAgetaccess(Afragment,frags,i,raw_toc_count) ++;
}

#undef VAgetaccess


#endif /* AS_CGB_METHODS_INCLUDE */
