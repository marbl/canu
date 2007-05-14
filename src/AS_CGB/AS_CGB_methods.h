
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
/*********************************************************************
 * Module: AS_CGB_methods.h
 *
 * Description: The methods to access the vertex and edge data store.
 *
 * Assumptions: That the user will be "good" and not access the data
 * without using the get_*, set_*, sum_*, and copy_* methods.
 *
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_METHODS_INCLUDE
#define AS_CGB_METHODS_INCLUDE


#include "AS_OVS_overlapStore.h"


typedef struct {
  IntFragment_ID  avx,bvx;
  int16       ahg,bhg; 

  uint32    quality : 16;  // zero is a perfect overlap
  uint32    nes : 8;      // The edge labeling.
  uint32    asx : 1;
  uint32    bsx : 1;
  uint32    reflected : 1; // Is this a direction mate edge.
  uint32    grangered : 1;
  uint32    invalid : 1; // Is this edge known to be invalid?
  uint32    blessed : 1;
  uint32    bit6 : 1;
  uint32    bit7 : 1;
} Aedge;


typedef struct { 
  uint64  src; /* An index into the character array
                      "the_source_text" holding the simulator source
                      info. */
  IntFragment_ID  iid; /* The Celera Genomics Assembler Gatekeeper IID 
		      for this fragment read. */
  IntChunk_ID cid; /* The linkage to chunks. Needed due to the way
                      that interchunk edges are made. */
  IntEdge_ID  segbgn_prefix;/* The starting location of a segment of 
			   raw fragment overlaps to the prefix. */
  IntEdge_ID  segbgn_suffix;/* The starting location of a segment of 
			   fragment overlaps to the suffix. */
  IntEdge_ID  segend_prefix;/* The ending location of a segment of 
			   raw fragment overlaps to the prefix. */
  IntEdge_ID  segend_suffix;/* The ending location of a segment of 
			   fragment overlaps to the suffix. */
  int32   nprefix_all;  /* The segment length of the prefix-only 
			   fragment overlaps. */
  int32   nsuffix_all;  /* The segment length of the suffix-only 
			   fragment overlaps. */
  int32   nprefix_dvt; /* Recomputable */
  // The segment length of the prefix-only fragment dovetail edges.
  int32   nsuffix_dvt; /* Recomputable */
  // The segment length of the suffix-only fragment dovetail edges.
  int32   nprefix_frc; /* Recomputable */
  // The segment length of the prefix-only fragment containment edges.
  int32   nsuffix_frc; /* Recomputable */
  // The segment length of the suffix-only fragment containment edges.


  int32   raw_prefix_dvt_count;
  // The segment length of the raw prefix-only fragment dovetail overlaps.
  int32   raw_suffix_dvt_count;
  // The segment length of the raw suffix-only fragment dovetail overlaps.
  int32   raw_frc_count;
  int32   raw_toc_count;
  
  IntFragment_ID  container;  /* Recomputable */
  // A zero value means that the fragment is not contained by any
  // other fragment.  A positive value is the fragment IID of the
  // container fragment.  This scheme depends on the Assembler IO
  // convention that Fragment IIDs are positive integers.

  BPTYPE      offset3p;
  BPTYPE      offset5p;
  FragType    frag_type;
  // The laboratory designation of the fragment type.

  Tlab        label; /* Recomputable */
  // The chunk labeling of fragments.
  // AS_CGB_THRU_FRAG, AS_INTERCHUNK_FRAG,
  // AS_CGB_INTRACHUNK_FRAG, AS_CGB_CONTAINED_FRAG,
  // AS_CGB_DELETED_FRAG, etc ?

  int16   bp_length;
  // The length of the fragment read in bp.
  
#ifdef STORE_BRANCH_POINTS_AT_FRAGMENT
  int16       pre_br;  /* Recomputable */
  int16       suf_br;  /* Recomputable */
  int16       pre_end;  /* Recomputable */
  int16       suf_end;  /* Recomputable */
  /* URT branch point info. The non-existence of a branch point is
     signaled by (pre_br==0)&&(suf_br==0). */
#endif 

  /* FGB bit flags: */
  unsigned int deleted : 1;
  // A flag indicating if this fragment is deleted in the fragment
  // graph.
  unsigned int contained : 1;  /* Recomputable */
  // A flag indicating if this fragment is contained in the fragment
  // graph.
  unsigned int spur : 1;

  unsigned int prefix_blessed : 1;
  unsigned int suffix_blessed : 1;
  unsigned int bit05 : 1; // Unused
  unsigned int bit06 : 1;
  unsigned int bit07 : 1; // Unused
  unsigned int bit08 : 1;
  unsigned int bit09 : 1;
  unsigned int bit10 : 1;
  unsigned int bit11 : 1;
  unsigned int bit12 : 1;
  unsigned int bit13 : 1;
  unsigned int bit14 : 1;
  unsigned int bit15 : 1;

} Afragment;

/******************************************************************/

typedef struct {
  IntFragment_ID vid;
} AChunkFrag;

typedef struct {
  BPTYPE        bp_length;  
  // the length in base pairs of the chunk.
  BPTYPE        rho;        
  IntChunk_ID   iaccession;
  // An arbitrary, but dense and non-negative enumeration of the
  // unitigs.
  float       coverage_stat;
  // Gene^s coverage statistic.
  BranchType a_branch_type, b_branch_type; 
  // the branch point characterization: AS_INTO_REPEAT,
  // AS_INTO_UNIQUE, or AS_NO_BPOINT.
  int32         a_branch_point;  // a bp distance from the A chunk-end.
  int32         b_branch_point;  // a bp distance from the B chunk-end.
  // branch_point == 0 means AS_NO_BPOINT
  // branch_point >  0 means AS_INTO_UNIQUE
  // branch_point <  0 means AS_INTO_REPEAT

  IntFragment_ID    chunk_avx; /* The A and B vertices of the chunk. */
  IntFragment_ID    chunk_bvx;
  int32     chunk_asx;
  int32     chunk_bsx;
  IntFragment_ID    num_frags;
  IntFragment_ID    f_list; /* The index into a TChunkFrag array. */

  // The following is not necessary. The info available from
  // fragment-end info, chunk_avx, chunk_asx, chunk_bvx, and chunk_bsx.
  int32         a_degree_raw; 
  int32         b_degree_raw;
  IntEdge_ID        a_list_raw;
  IntEdge_ID        b_list_raw;

  size_t        seq_loc; /* An index into a character array that
			 stores a consensus sequence of the chunk. */
  size_t        seq_len; // The length of the consensus sequence of the chunk.
  size_t        qua_loc; /* An index into a character array that
			 stores a quality sequence of the chunk. */
  size_t        qua_len; // The length of the quality sequence of the chunk.
  size_t        isrc; /* An index into a character array that
			 stores an annotation string about the chunk. */
  int8 asl,bsl; 
  // Flags indicating whether the A and B chunk-ends are probably
  // touching a tandem repeat region.
} AChunkMesg;


/* The following uses AS_UTL_Var.h. */
VA_DEF(Aedge)
VA_DEF(Afragment)
typedef VA_TYPE(Aedge)   Tedge;
typedef VA_TYPE(Afragment) Tfragment;
VA_DEF(AChunkFrag)
VA_DEF(AChunkMesg)
typedef VA_TYPE(AChunkFrag) TChunkFrag;
typedef VA_TYPE(AChunkMesg) TChunkMesg;
VA_DEF(IntEdge_ID)
typedef VA_TYPE(IntEdge_ID) TIntEdge_ID;

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
static IntEdge_ID GetNumEdges(const Tedge * const edges) {
  return (IntEdge_ID) GetNumVA_Aedge(edges);
}

#pragma inline GetNumFragments
static IntFragment_ID GetNumFragments(const Tfragment * const frags) {
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
#pragma inline get_buddy_dvt_edge
static int get_buddy_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{ return
    get_intrachunk_dvt_edge(edges,i) ||
    (AS_CGB_BUDDY_EDGE == get_nes_edge(edges,i));
}
#pragma inline get_interchunk_dvt_edge
static int get_interchunk_dvt_edge(const Tedge * const edges,IntEdge_ID i)
{ return
    get_buddy_dvt_edge(edges,i) ||
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
static void set_o3p_fragment(Tfragment frags[],IntFragment_ID i,BPTYPE value)
{ VAgetaccess(Afragment,frags,i,offset3p) = (BPTYPE) value;}
#pragma inline set_o5p_fragment
static void set_o5p_fragment(Tfragment frags[],IntFragment_ID i,BPTYPE value)
{ VAgetaccess(Afragment,frags,i,offset5p) = (BPTYPE) value;}

#pragma inline set_length_fragment
static void set_length_fragment(Tfragment frags[],IntFragment_ID i,int32 value)
{ VAgetaccess(Afragment,frags,i,bp_length) = (int16)value;}

#pragma inline set_src_fragment
static void set_src_fragment(Tfragment frags[],IntFragment_ID i,size_t value)
{ VAgetaccess(Afragment,frags,i,src) = (size_t)value;}

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
#pragma inline get_src_fragment
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

static BPTYPE get_o3p_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (BPTYPE) VAgetaccess(Afragment,frags,i,offset3p);}
static BPTYPE get_o5p_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (BPTYPE) VAgetaccess(Afragment,frags,i,offset5p);}

static int get_forward_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (VAgetaccess(Afragment,frags,i,offset3p) >
	  VAgetaccess(Afragment,frags,i,offset5p) );}
// Is this fragment in a forward orientation in its unitig?

static int32 get_length_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (int32) VAgetaccess(Afragment,frags,i,bp_length);}

static size_t get_src_fragment(const Tfragment * const frags,IntFragment_ID i)
{ return (size_t) VAgetaccess(Afragment,frags,i,src);}
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


#ifdef STORE_BRANCH_POINTS_AT_FRAGMENT
#pragma inline get_bpt_vertex
static int get_bpt_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag)
{ return (int) (flag ?
		       VAgetaccess(Afragment,frags,i,suf_br):
		       VAgetaccess(Afragment,frags,i,pre_br));}

#pragma inline set_bpt_vertex
static void set_bpt_vertex
(const Tfragment * const frags,IntFragment_ID i,int flag,int value)
{
  if(flag) { 
    VAgetaccess(Afragment,frags,i,suf_br) = value;
  } else {
    VAgetaccess(Afragment,frags,i,pre_br) = value;
  }
}
#endif // STORE_BRANCH_POINTS_AT_FRAGMENT

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

static int compare_edge_weak(const void * const aa, const void * const bb) 
{
  // This comparison function is used with ANSI qsort() for sorting
  // the edges to form contiguous segments.
  
  // The lesser edge is the one we keep in the Reaper.
  int icom;
  Aedge *a = (Aedge *)aa;
  Aedge *b = (Aedge *)bb;

  icom = ((a->avx) - (b->avx));
  if( icom == 0 ) {
    icom = ((a->asx) - (b->asx));
    //if( icom == 0 ) {
    //icom = (b->blessed - a->blessed);
      if( icom == 0 ) {
        icom = ((a->ahg) - (b->ahg));
        // Favor the minimum ahg.
        if( icom == 0 ) {
          icom = ((b->bhg) - (a->bhg));
          // Favor the maximum bhg.
          
          if( icom == 0 ) {
            // The following is unnecessary, but useful for the binary
            // search in the adjaceny lists and regression output.
            icom = ((a->bvx) - (b->bvx));
            if( icom == 0 )
              icom = ((a->bsx) - (b->bsx));
          }
        } // End of regression stuff.
      }
      //}
  } 
  return icom ;
}


static int compare_edge_strong(const void * const aa, const void * const bb) 
{
  // This comparison function is used with ANSI qsort() for sorting
  // the edges to form contiguous segments.
  
  // The lesser edge is the one we keep in the Reaper.
  int icom;
  Aedge *a = (Aedge *)aa;
  Aedge *b = (Aedge *)bb;

  icom = ((a->avx) - (b->avx));
  if( icom == 0 ) {
    icom = ((a->asx) - (b->asx));
    if( icom == 0 ) {
      icom = (b->blessed - a->blessed);
      if( icom == 0 ) {
        icom = ((a->ahg) - (b->ahg));
        // Favor the minimum ahg.
        if( icom == 0 ) {
          icom = ((b->bhg) - (a->bhg));
          // Favor the maximum bhg.
          
          if( icom == 0 ) {
            // The following is unnecessary, but useful for the binary
            // search in the adjaceny lists and regression output.
            icom = ((a->bvx) - (b->bvx));
            if( icom == 0 ) {
              icom = ((a->bsx) - (b->bsx));
              if( icom == 0 )
                icom = (a->reflected - b->reflected);
            }
          }
        } // End of regression stuff.
      }
    }
  }
  return icom ;
}

#endif /* AS_CGB_METHODS_INCLUDE */
