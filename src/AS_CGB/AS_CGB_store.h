
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
 * $Id: AS_CGB_store.h,v 1.5 2007-02-12 22:16:55 brianwalenz Exp $
 *
 * Module: AS_CGB_store.h
 * Description: Header file for the code that reads and writes the 
 * check point data.
 * Assumptions:
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_STORE_INCLUDE
#define AS_CGB_STORE_INCLUDE

#define SIMINFO
#undef SIMINFO

#undef GENINFO
#define GENINFO
#if 0
typedef cds_int64 BPTYPE;
#endif

typedef struct {
  /* fragment name */
#ifdef GENINFO
  BPTYPE  genbgn,genend; /* Coordinates in genome */
#endif
#ifdef SIMINFO
  int16 pre_brp,suf_brp,pre_end,suf_end,pre_ins,suf_ins;
  char  pre_let,suf_let;
# endif
} Afraginfo;

/* The following text uses AS_UTL_Var.h. */
VA_DEF(Afraginfo)
typedef VA_TYPE(Afraginfo) Tfraginfo;

#define VAgetaccess(Type,object,index,member) \
 (Get ## Type (object,index)->member )

/* Object access functions */

#ifdef GENINFO
static void set_genbgn_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,BPTYPE value)
{ VAgetaccess(Afraginfo,fraginfo,i,genbgn) = (BPTYPE)value;}
static void set_genend_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,BPTYPE value)
{ VAgetaccess(Afraginfo,fraginfo,i,genend) = (BPTYPE)value;}
#endif
#ifdef SIMINFO
static void set_pre_let_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,char value)
{ VAgetaccess(Afraginfo,fraginfo,i,pre_let) = (char)value;}
static void set_pre_ins_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,int16 value)
{ VAgetaccess(Afraginfo,fraginfo,i,pre_ins) = (int16)value;}
static void set_pre_brp_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,int16 value)
{ VAgetaccess(Afraginfo,fraginfo,i,pre_brp) = (int16)value;}
static void set_pre_end_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,int16 value)
{ VAgetaccess(Afraginfo,fraginfo,i,pre_end) = (int16)value;}
static void set_suf_let_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,char value)
{ VAgetaccess(Afraginfo,fraginfo,i,suf_let) = (char)value;}
static void set_suf_ins_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,int16 value)
{ VAgetaccess(Afraginfo,fraginfo,i,suf_ins) = (int16)value;}
static void set_suf_brp_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,int16 value)
{ VAgetaccess(Afraginfo,fraginfo,i,suf_brp) = (int16)value;}
static void set_suf_end_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,int16 value)
{ VAgetaccess(Afraginfo,fraginfo,i,suf_end) = (int16)value;}
#endif /*SIMINFO*/

#ifdef GENINFO
static BPTYPE get_genbgn_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (BPTYPE) VAgetaccess(Afraginfo,fraginfo,i,genbgn);}
static BPTYPE get_genend_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (BPTYPE) VAgetaccess(Afraginfo,fraginfo,i,genend);}
#endif /*GENINFO*/
#ifdef SIMINFO
static char get_pre_let_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (char) VAgetaccess(Afraginfo,fraginfo,i,pre_let);}
static int16 get_pre_ins_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (int16) VAgetaccess(Afraginfo,fraginfo,i,pre_ins);}
static int16 get_pre_brp_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (int16) VAgetaccess(Afraginfo,fraginfo,i,pre_brp);}
static int16 get_pre_end_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (int16) VAgetaccess(Afraginfo,fraginfo,i,pre_end);}
static char get_suf_let_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (char) VAgetaccess(Afraginfo,fraginfo,i,suf_let);}
static int16 get_suf_ins_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (int16) VAgetaccess(Afraginfo,fraginfo,i,suf_ins);}
static int16 get_suf_brp_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (int16) VAgetaccess(Afraginfo,fraginfo,i,suf_brp);}
static int16 get_suf_end_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (int16) VAgetaccess(Afraginfo,fraginfo,i,suf_end);}
#endif /*SIMINFO*/

#undef VAgetaccess

// End of fraginfo stuff... 


typedef struct { 
  cds_int32       store_version;
  cds_int32       state_of_the_store;
  cds_int32       unused;

  unsigned int    edges_sorted_by_fragment_end : 1;
  unsigned int    unmated_edges : 1;
  unsigned int    sorted_by_adjacency_lists : 1;
  unsigned int    next_edge_array_valid : 1;
  unsigned int    dechorded : 1;
  unsigned int    transitively_marked : 1;
  unsigned int    using_to_contained_edges : 1;
  unsigned int    bit07 : 1;
  unsigned int    bit08 : 1;
  unsigned int    bit09 : 1;
  unsigned int    bit10 : 1;
  unsigned int    bit11 : 1;
  unsigned int    bit12 : 1;
  unsigned int    bit13 : 1;
  unsigned int    bit14 : 1;
  unsigned int    bit15 : 1;
  unsigned int    bit16 : 1;
  unsigned int    bit17 : 1;
  unsigned int    bit18 : 1;
  unsigned int    bit19 : 1;
  unsigned int    bit20 : 1;
  unsigned int    bit21 : 1;
  unsigned int    bit22 : 1;
  unsigned int    bit23 : 1;
  unsigned int    bit24 : 1;
  unsigned int    bit25 : 1;
  unsigned int    bit26 : 1;
  unsigned int    bit27 : 1;
  unsigned int    bit28 : 1;
  unsigned int    bit29 : 1;
  unsigned int    bit30 : 1;
  unsigned int    bit31 : 1;

  BPTYPE          nbase_in_genome;
  IntFragment_ID  nfrag_randomly_sampled_in_genome;
  float           global_fragment_arrival_rate;
  // The estimated length of the genome in base pairs.
  IntFragment_ID  min_frag_iid;
  IntFragment_ID  max_frag_iid;

} TStateGlobals;


typedef struct { 

  Tfragment     *frags; /* The current fragment array handle. */
  Tedge         *edges;   /* The current edge array handle. */
  TIntEdge_ID   *next_edge_obj;
  Tfraginfo     *fraginfo;

  TChunkFrag    *chunkfrags;
  TChunkMesg    *thechunks;
  VA_TYPE(char) *chunkseqs;
  VA_TYPE(char) *chunkquas;

  VA_TYPE(char) *frag_annotations;
  VA_TYPE(char) *chunksrcs;

  //#ifdef MINIUNITIGGER
  VA_TYPE(OFGMesg)     * the_ofg_messages;
  VA_TYPE(OverlapMesg) * the_ovl_messages;
  VA_TYPE(char)        * the_ofg_source;
  VA_TYPE(char)        * the_ovl_source;
  VA_TYPE(IntMultiPos)   * the_imp_messages;
  VA_TYPE(IntUnitigMesg) * the_ium_messages;
  VA_TYPE(char)        * the_imp_source;
  VA_TYPE(char)        * the_ium_source;
  //#endif // MINIUNITIGGER
  
} THeapGlobals;

static void ReportHeapUsage_CGB
(
 TStateGlobals * gstate,
 THeapGlobals  * heapva,
 FILE          * fout 
)
{
  size_t total = 0;
  total += ReportMemorySize_VA( heapva->frags, "frags", fout);
  total += ReportMemorySize_VA( heapva->edges, "edges", fout);
  total += ReportMemorySize_VA( heapva->next_edge_obj, "next_edge_obj", fout);
  total += ReportMemorySize_VA( heapva->fraginfo, "fraginfo", fout);
  total += ReportMemorySize_VA( heapva->chunkfrags, "chunkfrags", fout);
  total += ReportMemorySize_VA( heapva->thechunks, "thechunks", fout);
  total += ReportMemorySize_VA( heapva->chunkseqs, "chunkseqs", fout);
  total += ReportMemorySize_VA( heapva->chunkquas, "chunkquas", fout);
  total += ReportMemorySize_VA( heapva->frag_annotations, "frag_annotations", fout);
  total += ReportMemorySize_VA( heapva->chunksrcs, "chunksrcs", fout);
  total += ReportMemorySize_VA( heapva->the_ofg_messages, "the_ofg_messages", fout);
  total += ReportMemorySize_VA( heapva->the_ovl_messages, "the_ovl_messages", fout);
  total += ReportMemorySize_VA( heapva->the_ofg_source, "the_ofg_source", fout);
  total += ReportMemorySize_VA( heapva->the_ovl_source, "the_ovl_source", fout);
  total += ReportMemorySize_VA( heapva->the_imp_messages, "the_imp_messages", fout);
  total += ReportMemorySize_VA( heapva->the_ium_messages, "the_ium_messages", fout);
  total += ReportMemorySize_VA( heapva->the_imp_source, "the_imp_source", fout);
  total += ReportMemorySize_VA( heapva->the_ium_source, "the_ium_source", fout);

  fprintf( fout,"ReportHeapUsage: " F_SIZE_T " bytes\n", total);
}


void open_fgb_store
( TStateGlobals * state,
  THeapGlobals * heapva);

void close_fgb_store
( TStateGlobals * state,
  THeapGlobals * heapva);

void read_fgb_store
( const char    * const theStorePath,
  TStateGlobals * gstate,
  THeapGlobals  * heapva,
  const size_t new_additional_number_of_frags,
  const size_t new_additional_number_of_edges,
  const size_t new_additional_amount_of_text);

void write_fgb_store
( const char * const theStorePath,
  const TStateGlobals * const gstate,
  const THeapGlobals  * const heapva
  );

void check_fgb_store
( const char * const theStorePath,
  const TStateGlobals * const gstate,
  const THeapGlobals  * const heapva
  );

void clear_fgb_store
(
  TStateGlobals * state,
  THeapGlobals  * heapva,
  int clear_the_frags,
  int clear_the_edges,
  int clear_the_frag_annotations,
  int clear_the_next_edge);

#endif /* AS_CGB_STORE_INCLUDE */
