#ifdef BRANCHPOINTS

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
= "$Id: AS_CGB_branchpts.c,v 1.9 2007-02-12 22:16:55 brianwalenz Exp $";
/* *******************************************************************
 *
 * Module: AS_CGB_branchpts.c
 * 
 * Description: 
 *
 * Assumptions: 
 *
 * Authors: Clark Mobarry and Art Delcher
 *
 *********************************************************************/

/*************************************************************************/
/* System include files */
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/*************************************************************************/
/* Local include files */
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_ALN_aligners.h"
#include "AS_CGB_all.h"
#include "AS_CGB_cgb.h"
#include "AS_CGB_branchpts.h"
#include "AS_CGB_histo.h"

/*************************************************************************/
/*************************************************************************/

/*************************************************************************/
/* Conditional compilation */

#define nsample 1000
#define nbucket 1000

#undef USE_MOTHER_FRAGMENT_GRAPH
#define USE_FGB_PATH

#undef DEBUGGING

#ifdef DEBUGGING
#define DEBUG_WITH_HISTO
#undef DONT_FIND_BRANCH_POINTS
#define DEBUG16
#define DEBUG31
#define DEBUG32
#undef DEBUG33
#undef DEBUG34
#undef DEBUG35
#undef DEBUG37
#define DEBUG71
#define DEBUG77
#endif /*DEBUGGING*/


#define  AS_CGB_BPT_SEQ_LENGTH (AS_READ_MAX_LEN + 2)

#define  NUM_BRANCH_COUNTS    10
    // Number of counters in following branch-point information structures
#define  MAX_BRANCH_CT     SHRT_MAX
    // Maximum possible values in branch-point counter

typedef  struct
{
  int16  place;     // extremal right-hand branch position
  int16  ct [NUM_BRANCH_COUNTS];
  // ct [NUM_BRANCH_COUNTS - 1] is # votes for branch at  place
  // ct [NUM_BRANCH_COUNTS - 2] is # votes for branch at  place - 1
  // ct [NUM_BRANCH_COUNTS - 3] is # votes for branch at  place - 2
  //   ...
}  BranchPtInfo_t;

typedef  struct {
  IntChunk_ID  cid; /* the chunk cid */
  IntFragment_ID vid; /* the fragment vid */
  int16  minhang;
  int16  maxhang;
  unsigned int csx : 1; /* the chunk prefix/suffix flag */
  unsigned int vsx : 1; /* the fragment prefix/suffix flag */

  unsigned int use_me : 1; /* TRUE if it passed the "useful" filter */
}  BranchOlap_t;

ReadStructp Branch_Pt_Read;    // From AS_PER_ReadStruct.h

//#define  ERROR_RATE        0.02
#define  ERROR_RATE        0.04
#define  BPT_E_RATE        2.0 * ERROR_RATE   
  // Set ERROR_RATE to expected difference
  // between fragments, e.g.,  0.04
#define  BPT_PROB_THOLD    1e-6
#define  HEADROOM_FACTOR   1.2
  // Extra memory factor when realloc'ing so don't realloc EVERY time.
#define  MAX_BP_CALLS      225
  // Most calls allowed to  BPnt_Seq_Comp_AS  in a single call to
  //  Find_Branch_Points .  If more fragments than this are given
  // to that function, then randomly select this many to call.

VA_DEF(BranchOlap_t)


/*************************************************************************/
/* Global Defines */

/*************************************************************************/

#ifndef USE_FGB_PATH
typedef IntChunk_ID IntGang_ID;
#else // USE_FGB_PATH
typedef IntFragment_ID IntGang_ID;
#endif // USE_FGB_PATH

typedef struct {
  IntGang_ID next;
  int min_offset, max_offset;
  int use_me;
} GangType;

VA_DEF(GangType)

static void gang_initialize(VA_TYPE(GangType) gangs[],IntGang_ID maxgangs) {
  IntGang_ID ig;
  for(ig=0;ig<maxgangs;ig++) {
    GangType agang;
    agang.next = ig; /* This denotes a gang with one member. */
    agang.min_offset = 0;
    agang.max_offset = 0;
    agang.use_me = FALSE;
    SetGangType(gangs,ig,&agang);
  }
}

static void gang_print(VA_TYPE(GangType) gangs[]) {
  const IntGang_ID ngangs = (IntGang_ID)GetNumVA_GangType(gangs);
  IntGang_ID ig;
  printf("ngangs=" F_IID "\n", ngangs);
  for(ig=0;ig<ngangs;ig++) {
    GangType *agang;
    agang = GetVA_GangType(gangs,ig);
    printf("% 10" F_IIDP " % 10" F_IIDP " % 6d % 6d % 3d\n",
	   ig,
	   agang->next,
	   agang->min_offset,
	   agang->max_offset,
	   agang->use_me);
  }
}


static void gang_cross(GangType * const agang,GangType * const bgang) {
  // gang_cross will join two disjoint circular lists,
  // or split a circular list into two circular lists.
  // For example, if A and B should be the "ends" of two 
  // circular lists, then gang_cross will perform the two 
  // cuts and heal both lists.
  IntGang_ID it;
  it = agang->next;
  agang->next = bgang->next;
  bgang->next = it;
}

static void gang_shift(VA_TYPE(GangType) gangs[],
		       const IntGang_ID ib,
		       const int min_shift,
		       const int max_shift) 
{
  GangType *bgang, *tgang;
  IntGang_ID it;
  bgang = GetVA_GangType(gangs,ib);
  bgang->min_offset += min_shift;
  bgang->max_offset += max_shift;
  for(it = bgang->next; (it != ib); it = tgang->next) {
    tgang = GetVA_GangType(gangs,it);
    tgang->min_offset += min_shift;
    tgang->max_offset += max_shift;
  }
}

static int gang_join
(VA_TYPE(GangType) gangs[],
 const IntGang_ID ia,
 const IntGang_ID ib,
 const int new_offset,
 const int use_me) {
  /* 

     This routine does a union-find while keeping track of the
     relative (linear only) coordinates of each vertex in the graph.
     The edge is A to B, where ia and ib are the vertex indicies.
     Here new_offset is an offset from A to B.  When a redundant edge between A
     and B is asserted, it is not assumed that the offset is the same
     as before. Instead, a range of offsets supplied for each vertex.
     
     Return TRUE  if there was a new joining.
     Return FALSE for an old joining.

     Adapted from Dr Dobbs Journal, January 1991. 
     The design is for small sets and to use a minimal amount of memory.
     Programmed by Clark Mobarry in April 1999.  

*/

  int a_not_seen, b_not_seen, iret = TRUE;
  IntGang_ID maxgangs;

  GangType *agang = NULL, *bgang = NULL, *tgang = NULL;
  assert(gangs != NULL);
  maxgangs = (IntGang_ID)GetNumVA_GangType(gangs);
  assert( maxgangs > 0);
  assert( ia >= 0 );
  assert( ia < maxgangs);
  assert( ib >= 0 );
  assert( ib < maxgangs);
  agang = GetVA_GangType(gangs,ia);
  bgang = GetVA_GangType(gangs,ib);

  bgang->use_me = use_me;

  if( ia == ib ) {
    bgang->min_offset = MIN(bgang->min_offset,new_offset);
    bgang->max_offset = MAX(bgang->max_offset,new_offset);
    return 0;
  }

  a_not_seen = (agang->next == ia);
  b_not_seen = (bgang->next == ib);

  /* Consider four distinct choices. */
  if( a_not_seen && b_not_seen) {
    /* Create a new set with just A and B. */
    agang->min_offset = 0; agang->max_offset = 0; 
    bgang->min_offset = new_offset; bgang->max_offset = new_offset; 
    gang_cross(agang,bgang);
  }
  if( (!a_not_seen) && b_not_seen) {
    /* Append B to A^s set. */
    bgang->min_offset = agang->min_offset + new_offset; 
    bgang->max_offset = agang->max_offset + new_offset; 
    gang_cross(agang,bgang);
  }
  if( (!b_not_seen) && a_not_seen) {
    /* Append A to B^s set. */
    int min_shift,max_shift;
    agang->min_offset = 0; agang->max_offset = 0; 
    min_shift = new_offset + agang->min_offset - bgang->max_offset;
    max_shift = new_offset + agang->max_offset - bgang->min_offset;
    /* First move the offsets of the gang of B relative to the gang of A. */
    gang_shift(gangs,ib,min_shift,max_shift);
    gang_cross(agang,bgang);
  }
  if( (!a_not_seen) && (!b_not_seen)) {
    /* Determine if A and B are already in the same set. */
    IntGang_ID it;
    int same_set; 
    for(it = agang->next; 
	((it != ia)&&(it != ib));
	it = GetVA_GangType(gangs,it)->next);
    same_set = (it == ib);

    if(same_set) { /* Just adjust the offset spread of the pair. */
      int y_min,y_max;
      y_min = new_offset + agang->min_offset;
      y_max = new_offset + agang->max_offset;

      bgang->min_offset = MIN(bgang->min_offset,y_min);
      bgang->max_offset = MAX(bgang->max_offset,y_max);

      iret = 0;
    } else { /* Merge the two sets */
      int min_shift,max_shift;
      min_shift = new_offset + agang->min_offset - bgang->max_offset;
      max_shift = new_offset + agang->max_offset - bgang->min_offset;
      /* First move the offsets of the gang of B relative to the gang of A. */
      gang_shift(gangs,ib,min_shift,max_shift);
      /* Second join the two circular lists. */
      gang_cross(agang,bgang);
    }
  }
  return iret;
}


static void gang_clear(VA_TYPE(GangType) gangs[], IntGang_ID ia) {
  /* Clear the gang identity of the gang that includes key ia. */
  IntGang_ID it, next= ~((IntGang_ID)0), maxgangs=0;
  GangType *agang = NULL;

  assert( gangs != NULL);
  maxgangs = (IntGang_ID)GetNumVA_GangType(gangs);
  assert( maxgangs > 0);
  assert( ia >= 0 );
  assert( ia < maxgangs);
  agang = GetVA_GangType(gangs,ia);
  assert(agang != NULL);
  for(it = agang->next; (it != ia); it = next) {
    GangType *tgang = GetVA_GangType(gangs,it);
    next = tgang->next;
    tgang->next = it;
    tgang->min_offset = 0;
    tgang->max_offset = 0;
    tgang->use_me = FALSE;
  }
  agang = GetVA_GangType(gangs,ia);
  agang->next = ia;
  agang->min_offset = 0;
  agang->max_offset = 0;
  agang->use_me = FALSE;
}



/*****************************************************************/

static char Complement(char Ch)
/*  Return the Watson-Crick complement of  Ch . */
{
  switch((int) Ch) {
  case  'A' :
    return  'T';
  case  'C' :
    return  'G';
  case  'G' :
    return  'C';
  case  'T' :
    return  'A';
  case  'N' :
    return  'N';
  case  'X' :
    return  'X';
  default :
    fprintf(stderr, __FILE__ "ERROR(complement):"
	    "Unexpected DNA character `%c\'\n", Ch);
    assert(FALSE);
  }
  return (char)-1;
}

static void  Rev_Complement
    (char * S, int Len)
/* Set  S [0 .. Len - 1]  to it^s Watson-Crick (reverse) complement. */
{
  char  Ch;
  int  i, j;
  
  for  (i = 0, j = Len - 1;  i < j;  i ++, j --) {
    Ch = Complement (S [i]);
    S [i] = Complement (S [j]);
    S [j] = Ch;
  }
  
  if(i == j){
    S [i] = Complement (S [i]);
  }
  return;
}


/*****************************************************************/


static void  Add_Branch_Point_To_Bin
    (int pt, BranchPtInfo_t * b)
     
     // Add branch point at position  pt  to branch information in b.
     // In b is NUM_BRANCH_COUNTS bins that keeps a distribution of 
     // observed branch points.
  {
   int  i, diff;

   if  (b->place == 0     // indicates empty BranchPtInfo
	|| pt >= b->place + NUM_BRANCH_COUNTS)
       {
        b->place = pt;
        memset( b->ct, 0, NUM_BRANCH_COUNTS * sizeof (int16));
        b->ct[0] = 1;
        return;
       }

   if  (pt <= b->place - NUM_BRANCH_COUNTS)
       return;

   if  (pt <= b->place)
       {
        int16  * p = b->ct + (b->place - pt);

        if  (* p < MAX_BRANCH_CT)
            (* p) ++;
        return;
       }

   diff = pt - b->place;
   for(i = NUM_BRANCH_COUNTS - 1;  i >= diff;  i --)
     { b->ct[i] = b->ct[i - diff];}
   for( ;  i > 0;  i --)
     { b -> ct [i] = 0;}
   b->ct[0] = 1;
   b->place = pt;

   return;
  }

static void  Get_Chunk_End_Sequence
(char * const s, int max_len,
 BranchOlap_t * f, 
 TChunkMesg * thechunks, 
 // the handle to the chunks
 VA_TYPE(char) *chunkseqs,
 // the store for the consensus sequence of the chunks.
 const Tfragment * frags,
 FragStoreHandle frag_store,
 const int use_consensus
)
{
  const IntChunk_ID chunk_index = f->cid;
  const int chunk_suffix = f->csx;
#ifndef USE_FGB_PATH
  AChunkMesg * mychunk = GetVA_AChunkMesg(thechunks,chunk_index);
  assert(mychunk != NULL);
#else
  AChunkMesg * mychunk = NULL;
  assert(FALSE == use_consensus);
#endif
  
  if(use_consensus) {  // A consensus orientated branch point finder.
    const size_t seq_loc = mychunk->seq_loc;
    const size_t seq_len = mychunk->seq_len;
    const char * const sequence = GetVA_char(chunkseqs,seq_loc);
    int out_len;

#ifdef DEBUG34
    printf("chunk_index=" F_IID ",strlen(sequence)=" F_SIZE_T ",seq_len=" F_SIZE_T "\n",
	   chunk_index,strlen(sequence),seq_len);
#endif
    assert( strlen(sequence) == seq_len);
    out_len = MIN(seq_len,max_len);
    
    if(chunk_suffix) {
      strncpy(s, &(sequence[seq_len-out_len]), out_len);
      s[out_len] = '\0';
      Rev_Complement(s, out_len);
    } else {
      strncpy(s, sequence, out_len);
      s[out_len] = '\0';
    }
  } else { // A fragment orientated branch point finder.
    char  sequence[AS_CGB_BPT_SEQ_LENGTH+1];
    char  quality[AS_CGB_BPT_SEQ_LENGTH+1];
    IntFragment_ID frag_vid, frag_iid;
    int frag_suf, iret;
    size_t seq_len;
    uint clear_start, clear_end;
#ifndef USE_FGB_PATH
    if(chunk_suffix == 0) {
      frag_vid = mychunk->chunk_avx; frag_suf = mychunk->chunk_asx;
    } else {
      frag_vid = mychunk->chunk_bvx; frag_suf = mychunk->chunk_bsx;
    }
#else
    frag_vid = f->vid; frag_suf = f->vsx;
#endif
    frag_iid = get_iid_fragment(frags,frag_vid);
    
    iret = getFragStore( frag_store, frag_iid, FRAG_S_SEQUENCE, Branch_Pt_Read);
    assert(iret == 0);
    iret = getSequence_ReadStruct
      (Branch_Pt_Read, sequence, quality, AS_CGB_BPT_SEQ_LENGTH);
    assert(iret == 0); // Will fail if buffer is to short.
    getClearRegion_ReadStruct( Branch_Pt_Read, & clear_start, & clear_end,
			       READSTRUCT_LATEST);
    seq_len = clear_end - clear_start;
#ifdef DEBUG34
    printf("frag_index=" F_IID ",frag_suf=%d,strlen(sequence)=" F_SIZE_T ",seq_len=" F_SIZE_T "\n",
	   get_iid_fragment(frags,frag_vid),
	   frag_suf,strlen(sequence),seq_len);
#endif
    assert( seq_len <= max_len );
    // memmove(s, sequence + clear_start, seq_len);
    strncpy(s, &(sequence[clear_start]), seq_len);
    s[seq_len] = '\0';
    
    if(frag_suf) { Rev_Complement(s, seq_len);}
  }
#ifdef DEBUG35
    printf("chunk_index=" F_IID ",chunk_suffix=%d,strlen(s)=" F_SIZE_T ",s=\n%s\n",
	   chunk_index,chunk_suffix,strlen(s),s);
#endif
    
  return;
}

static void  Get_Fragment_End_Sequence
(char * const s, int max_len,
 BranchOlap_t * f, 
 const Tfragment * frags,
 FragStoreHandle frag_store
)
{
  const IntChunk_ID chunk_index = f->cid;
  const int chunk_suffix = f->csx;

 // A fragment orientated branch point finder.
  char  sequence[AS_CGB_BPT_SEQ_LENGTH+1];
  char  quality[AS_CGB_BPT_SEQ_LENGTH+1];
  IntFragment_ID frag_vid = f->vid;
  int frag_suf = f->vsx;
  int iret;
  size_t seq_len;
  uint clear_start, clear_end;
  IntFragment_ID frag_iid = get_iid_fragment(frags,frag_vid);
  
  iret = getFragStore( frag_store, frag_iid, FRAG_S_SEQUENCE, Branch_Pt_Read);
  assert(iret == 0);
  iret = getSequence_ReadStruct
    (Branch_Pt_Read, sequence, quality, AS_CGB_BPT_SEQ_LENGTH);
  assert(iret == 0); // Will fail if buffer is to short.
  getClearRegion_ReadStruct( Branch_Pt_Read, & clear_start, & clear_end,
			     READSTRUCT_LATEST);
  seq_len = clear_end - clear_start;
#ifdef DEBUG34
  printf("Before Watson-Crick complementation:\n"
	 "chunk_index=" F_IID " chunk_suffix=%d\n"
	 "frag_index=" F_IID ",frag_suf=%d,strlen(sequence)=" F_SIZE_T ",seq_len=" F_SIZE_T "\n",
	 chunk_index,chunk_suffix,
	 get_iid_fragment(frags,frag_vid),
	 frag_suf,strlen(sequence),seq_len);
#endif
  assert( seq_len <= max_len );
  // memmove(s, sequence + clear_start, seq_len);
  strncpy(s, &(sequence[clear_start]), seq_len);
  s[seq_len] = '\0';
  
  if(frag_suf) { Rev_Complement(s, seq_len);}

#ifdef DEBUG35
  printf("After Watson-Crick complementation:\n"
	 "chunk_index=" F_IID " chunk_suffix=%d\n"
	 "frag_index=" F_IID ",frag_suf=%d,strlen(sequence)=" F_SIZE_T ",seq_len=" F_SIZE_T "\n",
	 chunk_index,chunk_suffix,
	 get_iid_fragment(frags,frag_vid),
	 frag_suf,strlen(sequence),seq_len);
#endif
    
  return;
}

//  Routines borrowed from OVL to get distributions of counts

typedef  struct Distribution
{
  int  N;
  int  * Ct;
  double  * Thold, Max, Sum;
}  Distrib_t;

//  Holds  N+1  counts.  For  0 <= i <= N ,  Ct [i]  is the
//  number of values that were  <= Thold [i] .


static void  Init_Distrib  (Distrib_t * D, int N)

//  Constructor to allocate memory for a new  D  with  D -> N = N .

{
  int  i;
  
  D -> Max = - FLT_MAX;
  D -> Sum = 0.0;
  D -> N = N;
  D->Ct    = safe_calloc(sizeof(int), N + 1);
  D->Thold = safe_calloc(sizeof(double), N + 1);
  assert (D -> Thold != NULL);
  
  return;
}


static void  Incr_Distrib  (Distrib_t * D, double Val)
//  Increment the count in  D  that corresponds to  Val .
{
  int  i;
  if  (Val > D -> Max)
    { D -> Max = Val;}
  D -> Sum += Val;
  for  (i = 0;  i < D -> N;  i ++) {
    if  (Val <= D -> Thold [i]) {
      D -> Ct [i] ++;
      return;
    }
  }
  D -> Ct [D -> N] ++;
  return;
}


static void  Print_Distrib  (Distrib_t D, char * Heading)
//  Display the values in  D  on the standard error stream.
{
  int64  Total = 0;
  int  i;
  
  fprintf (stderr, "\n%s\n", Heading);
  fprintf (stderr, "     <=     Count\n");
  for  (i = 0;  i < D . N;  i ++) {
    fprintf (stderr, "%7.0f  %8d\n", D . Thold [i], D . Ct [i]);
    Total += D . Ct [i];
  }
  fprintf (stderr, "  Above   %7d\n", D . Ct [D . N]);
  Total += D . Ct [D . N];
  fprintf (stderr, "  Total %9" F_S64P "\n", Total);
  if  (Total > 0) {
    fprintf (stderr, "  Max   %9.0f\n", D . Max);
    fprintf (stderr, "  Avg   %11.1f\n", D . Sum / Total);
  }
  return;
}


static Distrib_t  BP_List_Len_Distrib;
// static FILE  * BP_Log_File;

static void  Find_Branch_Points
(int n, 
 BranchOlap_t * brncholaps, 
 BranchPtInfo_t * branch_pt,
 TChunkMesg * thechunks,
 VA_TYPE(char) * chunkseqs,
 const Tfragment * frags,
 FragStoreHandle frag_store,
 const int use_consensus,
 FILE *fbpts )
     
//  Add to  (* branchpt)  the branch-point positions for  brncholaps[0]
//  that are found by overlapping it with fragments  brncholaps[1 .. (n-1)] .
//  Fragment sequences are obtained from frag_store.
//  Assumes that  branch_pt  is already initialized.  If empty then the
//    place  field should be 0 (ct entries don't matter).  This routine
//    only adds branch points to the existing record.
//  Note:  if  brncholaps[0]  is reversed, all calculations and coordinates
//    are done after reverse-complementing it, i.e., in the reverse-complement
//     coordinate frame.

{
  int  i;
  char aseq[AS_CGB_BPT_SEQ_LENGTH+1];
  size_t  alen;
  static int  * avail = NULL;
  static int  avail_len = 0;
  static char  * chosen = NULL;
  
  if(avail_len == 0) {
    avail_len = (int)(MAX_BP_CALLS * HEADROOM_FACTOR);
    avail  = safe_malloc(sizeof(int) * avail_len);
    chosen = safe_malloc(sizeof(char) * avail_len);
    
    Init_Distrib (& BP_List_Len_Distrib, 38);
    for  (i = 0;  i < 10;  i ++)
      BP_List_Len_Distrib . Thold [i] = (i + 1);
    for  (;  i < 19;  i ++)
      BP_List_Len_Distrib . Thold [i] = 10.0 * (i - 8);
    for  ( ;  i < 28;  i ++)
      BP_List_Len_Distrib . Thold [i] = 100.0 * (i - 17);
    for  ( ;  i < 37;  i ++)
      BP_List_Len_Distrib . Thold [i] = 1000.0 * (i - 26);
  }
  Incr_Distrib (& BP_List_Len_Distrib, n);
  
  aseq[0] = '@'; // Gene^s routines do not use location zero.
  Get_Chunk_End_Sequence(aseq + 1, AS_CGB_BPT_SEQ_LENGTH, &(brncholaps[0]), 
			 thechunks, chunkseqs, frags, frag_store, use_consensus);
  alen = strlen(aseq + 1);
#ifdef DEBUG37
  printf("alen=" F_SIZE_T ",aseq=\n%s\n",alen,aseq);
#endif
  if(n > avail_len) {
    avail_len = (int)(n * HEADROOM_FACTOR);
    avail = safe_realloc(sizeof(int) avail_len);
    chosen = safe_realloc(sizeof(char) avail_len);
  }
  
  assert (brncholaps[0].use_me);
  if(n <= MAX_BP_CALLS + 1) {
    memset (chosen, '\1', n);
  } else {
    int last = 0;
    
    for(i = 1;  i < n;  i ++) {
      if(brncholaps[i].use_me) {
	avail [last ++] = i;}}
    
    if(last <= MAX_BP_CALLS) {
      memset (chosen, '\1', n);
    } else {
      memset (chosen, '\0', n);
      for(i = 0;  i < MAX_BP_CALLS;  i ++) {
	const int sub = (int) (lrand48 () % last);
	chosen [avail [sub]] = 1;
	avail [sub] = avail [-- last];
      }
    }  
  }
  
  for(i = 1;  i < n;  i ++)
    if(brncholaps[i].use_me && chosen [i]) {
      
      BranchPointResult  * branch_result = NULL;
      char  bseq [AS_CGB_BPT_SEQ_LENGTH+1];
      size_t  blen;
      
      bseq[0] = '@'; // Gene^s routines do not use location zero.
      Get_Chunk_End_Sequence(bseq + 1,AS_CGB_BPT_SEQ_LENGTH, &(brncholaps[i]), 
			     thechunks, chunkseqs, frags, frag_store, use_consensus);
      blen = strlen(bseq + 1); 
#ifdef DEBUG37
      printf("blen=" F_SIZE_T ",bseq=\n%s\n",blen,bseq);
#endif
      // Gene^s dynamic programming uses place 0 for the first
      // interface. thus the first base is at subscript one.
      
      branch_result = BPnt_Seq_Comp_AS
	(aseq, alen, bseq, blen, brncholaps[i].minhang,
	 brncholaps[i].maxhang, BPT_E_RATE, BPT_PROB_THOLD,
	 BPT_MIN_PREFIX, BPT_MIN_SUFFIX);
      
      // The return value is -1, if a branch point could not be found.
      if  (branch_result == NULL)
	{
	  fprintf (stderr, "ERROR: "
		   "Insufficient memory for branch point\n");
	  exit (EXIT_FAILURE);
	}
      
#if 0
      {
	fprintf 
	  (fbpts,
           "DIAG: " F_IID " %d " F_IID " %d " F_IID " %d " F_IID " %d "
	   "Branch pt at place=(%d,%d), ascent=%f, descent=%f.\n",
	   brncholaps[0].cid, brncholaps[0].csx,
	   brncholaps[i].cid, brncholaps[i].csx,
	   get_iid_fragment(frags,brncholaps[0].vid), brncholaps[0].vsx,
	   get_iid_fragment(frags,brncholaps[i].vid), brncholaps[i].vsx,
	   branch_result->apnt,
	   branch_result->bpnt,
	   branch_result->ascent,
	   branch_result->descent);
	// fflush(fbpts);
      }
#endif

      if((branch_result -> apnt > 0)
	 && (branch_result -> ascent > 0.0)
	 && (branch_result -> descent > 0.1) 
	 )
	{
	  Add_Branch_Point_To_Bin(branch_result->apnt, branch_pt);
	}
    }
  return;
}



static void  Dump_BP_List_Len_Distrib
(void)

//  Print distribution of list lengths  n  sent to preceding
//  function.

{
  Print_Distrib (BP_List_Len_Distrib, "Partner List Lengths:");
  return;
}



static void  Init_Branch_Compute
(void)

//  Allocate memory for structure to hold what's read from the frag_store.
//  Should be done once and re-used for all frag_store reads.
//  This memory is not freed.
    
{
  Branch_Pt_Read = new_ReadStruct (); // From AS_PER_ReadStruct.h
  //   BP_Log_File = fopen ("cgb-bpinfo.log", "w");
}

/*********************************************************************/

static int Find_Branch_Partners
( /* input only */ 
 const IntChunk_ID seed_chunk_index,
 const int seed_chunk_suffix,
 const Tfragment frags[],
 const Tedge edges[],
 const TChunkMesg thechunks[],
 /* clear and output */
 VA_TYPE(GangType) gangs[],
 /* output only */
 int * const branch_mate_ct )
{
  /* The return value is the number of branch partners including the
     seed chunk-end. */

  const AChunkMesg * the_seed_chunk 
    = GetVA_AChunkMesg(thechunks,seed_chunk_index);
  
  /* From the seed chunk in the gang, walk the adjacent vertices. */
  const IntEdge_ID  seed_chunkedgelist
    = ( seed_chunk_suffix == 0 ? 
	the_seed_chunk->a_list_raw :
	the_seed_chunk->b_list_raw );
  const int seed_chunkedgedegree 
    = ( seed_chunk_suffix == 0 ? 
	the_seed_chunk->a_degree_raw :
	the_seed_chunk->b_degree_raw );

#ifndef USE_FGB_PATH
  const IntGang_ID seed_gang_key 
    = 2*seed_chunk_index + seed_chunk_suffix;
#else // USE_FGB_PATH
  const IntGang_ID seed_gang_key 
    = 2*( seed_chunk_suffix == 0 ? 
	  the_seed_chunk->chunk_avx :
	  the_seed_chunk->chunk_bvx )
    +   ( seed_chunk_suffix == 0 ? 
	  the_seed_chunk->chunk_asx :
	  the_seed_chunk->chunk_bsx );
#endif // USE_FGB_PATH

  int imate_index, imate_used;
  int num_branch_partners=0;
  *branch_mate_ct = seed_chunkedgedegree; 
#ifdef DEBUG16
  printf("Find_Branch_Partners:seed_chunk_index=" F_IID ",seed_chunk_suffix=%d\n",
	 seed_chunk_index,seed_chunk_suffix);
#endif
  num_branch_partners++;

#ifdef DEBUG71
  { // assert that the adjacency list is already sorted by increasing ahg.
    for(imate_index=1; imate_index < seed_chunkedgedegree; imate_index++) {
      assert(
	     get_best_ovl(frags,edges,seed_chunkedgelist + imate_index-1) >=
	     get_best_ovl(frags,edges,seed_chunkedgelist + imate_index) );
    }
  }
#endif /*DEBUG71*/

  /* From the seed chunk-end find its mate chunk-ends. */
  for(imate_index=0, imate_used=0; 
      (imate_index < seed_chunkedgedegree) && (imate_used < MAX_NUM_OF_MATES); 
      imate_index++) {
    // The adjacency list is already sorted by overhang when the original
    // fragment adjacency list, that is from
    // smallest to largest. We are taking just the MAX_NUM_OF_MATES
    // 
    const IntEdge_ID mate_edge = seed_chunkedgelist + imate_index;
    const Tnes       mate_edge_nes = get_nes_edge(edges,mate_edge);
    const int use_me = TRUE;
    if( (AS_CGB_INTERCHUNK_EDGE == mate_edge_nes)
	) {
      const IntChunk_ID mate_chunk_index 
	= get_chunk_index(thechunks,frags,edges,mate_edge);
      const int mate_chunk_suffix 
	= get_chunk_suffix(thechunks,frags,edges,mate_edge);

      /* From the mate chunk, walk the adjacent vertices. */
      IntEdge_ID mate_chunkedgelist = 
	( mate_chunk_suffix == FALSE ? 
	  GetVA_AChunkMesg(thechunks,mate_chunk_index)->a_list_raw :
	  GetVA_AChunkMesg(thechunks,mate_chunk_index)->b_list_raw );
      int mate_chunkedgedegree = 
	( mate_chunk_suffix == FALSE ? 
	  GetVA_AChunkMesg(thechunks,mate_chunk_index)->a_degree_raw :
	  GetVA_AChunkMesg(thechunks,mate_chunk_index)->b_degree_raw );

      int ikin_index=0, ikin_used=0;

#ifdef DEBUG16
      printf("Find_Branch_Partners:mate_chunk_index=" F_IID ",mate_chunk_suffix=%d\n",
	     mate_chunk_index,mate_chunk_suffix);
#endif
      
      imate_used++;
      /* From the mate chunk-end find the seed^s kin chunk-ends. */
      for(ikin_index=0, ikin_used=0;
	  (ikin_index < mate_chunkedgedegree) && (ikin_used < MAX_NUM_OF_KIN);
	  ikin_index++) 
	if(mate_chunkedgedegree > 1){
	  int bp_offset;
	  
	  const IntEdge_ID kin_edge = mate_chunkedgelist + ikin_index;
	  const Tnes   kin_edge_nes = get_nes_edge(edges,kin_edge);
	  if( (AS_CGB_INTERCHUNK_EDGE == kin_edge_nes)
	      ) {
	    const IntChunk_ID kin_chunk_index 
	      = get_chunk_index(thechunks,frags,edges,kin_edge);
	    const int kin_chunk_suffix 
	      = get_chunk_suffix(thechunks,frags,edges,kin_edge);
#ifndef USE_FGB_PATH
	    const IntGang_ID kin_gang_key 
	      = 2*kin_chunk_index + kin_chunk_suffix;
#else // USE_FGB_PATH
	    const IntGang_ID kin_gang_key =
#if 0
	      2*( kin_chunk_suffix == 0 ? 
		  the_kin_chunk->chunk_avx :
		  the_kin_chunk->chunk_bvx )
	      +   ( seed_chunk_suffix == 0 ? 
		    the_kin_chunk->chunk_asx :
		    the_kin_chunk->chunk_bsx )
#else
	      2 * get_bvx_edge(edges,kin_edge) + get_bsx_edge(edges,kin_edge);
#endif
#endif // USE_FGB_PATH

	    ikin_used++;
	    
#ifdef DEBUG16
	    printf("Find_Branch_Partners:"
		   "kin_chunk_index=" F_IID ",kin_chunk_suffix=%d\n",
		   kin_chunk_index,kin_chunk_suffix);
#endif      
	    bp_offset 
	      = get_best_ovl(frags,edges,mate_edge)
	      - get_best_ovl(frags,edges,kin_edge);
	    //  The aligner code needs to take that fragment-end of
	    //  the seed and kin fragment that have common sequence,
	    //  that is in the repetitve region.  Thus a positive
	    //  bp_offset of the kin fragment-end with repect to the
	    //  seed fragment-end looks like this.

	    //  seed     RRRRRRRRRSSSSSSSSSSSS 
	    //  mate  RRRRRRRRRR  
	    //  kin        RRRRRRRKKKKKKKKK     

	    
#ifdef DEBUG16
	    printf("mate_edge aid,asx,ahg=" F_IID " %d %d\n",
		   get_iid_fragment(frags,get_avx_edge(edges,mate_edge)),
		   get_asx_edge(edges,mate_edge),
		   get_ahg_edge(edges,mate_edge));
	    printf("mate_edge bid,bsx,bhg=" F_IID " %d %d\n",
		   get_iid_fragment(frags,get_bvx_edge(edges,mate_edge)),
		   get_bsx_edge(edges,mate_edge),
		   get_bhg_edge(edges,mate_edge));

	    printf("kin_edge aid,asx,ahg=" F_IID " %d %d \n",
		   get_iid_fragment(frags,get_avx_edge(edges,kin_edge)),
		   get_asx_edge(edges,kin_edge),
		   get_ahg_edge(edges,kin_edge));
	    printf("kin_edge bid,bsx,bhg=" F_IID " %d %d\n",
		   get_iid_fragment(frags,get_bvx_edge(edges,kin_edge)),
		   get_bsx_edge(edges,kin_edge),
		   get_bhg_edge(edges,kin_edge));

	    printf("get_best_ovl(frags,edges,mate_edge)=%d\n",
		   get_best_ovl(frags,edges,mate_edge));
	    printf("get_best_ovl(frags,edges,kin_edge)=%d\n",
		   get_best_ovl(frags,edges,kin_edge));
	    printf("bp_offset=%d\n",bp_offset);
#endif
	    /* A positive bp_offset corresponds to the seed chunk end
	       overhanging the kin chunk end. */
	    
	    // use_me = (num_branch_partners < MAX_NUM_OF_PARTNERS);
	    if( (seed_gang_key != kin_gang_key) &&
		(num_branch_partners < MAX_NUM_OF_PARTNERS)) {
	      num_branch_partners +=
#ifndef USE_MOTHER_FRAGMENT_GRAPH
		gang_join(gangs, seed_gang_key, kin_gang_key, 
			  bp_offset, use_me);
#else // USE_MOTHER_FRAGMENT_GRAPH
#endif // USE_MOTHER_FRAGMENT_GRAPH

	    }
	  }
	}
    }
  }
  
  return num_branch_partners;
}


void find_the_branch_points
(
 /* input only */
 const int           use_consensus,
 const float         cgb_unique_cutoff,
 const float         global_fragment_arrival_rate,
 FragStoreHandle     frag_store,
 const Tfragment     frags[],
 const Tedge         edges[],
 TChunkFrag          chunkfrags[],
 VA_TYPE(char)       chunkseqs[],
 /* modify */
 TChunkMesg          thechunks[],
 FILE *fbpts
)
{
  VA_TYPE(BranchOlap_t) *theBranchOlaps = NULL;
  VA_TYPE(GangType)     *gangs = NULL;
  BranchPtInfo_t branch_pt;
  const IntChunk_ID nchunks = GetNumVA_AChunkMesg(thechunks);
  IntChunk_ID seed_chunk_index;

  const IntFragment_ID nfrag = GetNumFragments(frags);
#ifndef USE_FGB_PATH
  const IntGang_ID maxgangs = 2*(nchunks+1);
#else // USE_FGB_PATH
  const IntGang_ID maxgangs = 2*(nfrag+1);
#endif // USE_FGB_PATH

#ifdef DEBUG_WITH_HISTO
  FILE * fout = stdout;
  HISTOGRAM 
    *histo_branch_mates = NULL,
    *histo_branch_olaps = NULL,
    *histo_branch_mate_times_olaps = NULL;
  
  histo_branch_mates = create_histogram(nsample,nbucket,0,TRUE);
  histo_branch_olaps = create_histogram(nsample,nbucket,0,TRUE);
  histo_branch_mate_times_olaps = create_histogram(nsample,nbucket,0,TRUE);
#endif /*DEBUG_WITH_HISTO*/

  branch_pt.place = 0;
  Init_Branch_Compute();

  theBranchOlaps = CreateVA_BranchOlap_t(0);
  gangs = CreateVA_GangType(maxgangs);
  gang_initialize(gangs,maxgangs);

  fprintf(stderr,
          "CGB: find_the_branch_points() cgb_unique_cutoff=%f\n",
          cgb_unique_cutoff);

  /* Loop over chunk ends as vertices of the chunk graph. */
  for(seed_chunk_index=0;
      seed_chunk_index < nchunks; seed_chunk_index++) {
    int seed_chunk_suffix;
    
    const int number_of_randomly_sampled_fragments_in_chunk
      = count_the_randomly_sampled_fragments_in_a_chunk
      ( frags, chunkfrags, thechunks, seed_chunk_index);
    const BPTYPE rho = GetVA_AChunkMesg(thechunks,seed_chunk_index)->rho;
    const float coverage_stat = compute_coverage_statistic
      ( rho,
        number_of_randomly_sampled_fragments_in_chunk,
        global_fragment_arrival_rate );
    
        
    if( coverage_stat >= cgb_unique_cutoff)
      for(seed_chunk_suffix=0; seed_chunk_suffix<2; seed_chunk_suffix++) {
        IntGang_ID seed_gang_key;
        int branch_olap_ct = 0;
        int branch_mate_ct = 0;
        
#ifndef USE_FGB_PATH
        seed_gang_key = 2*seed_chunk_index+seed_chunk_suffix;
#else // USE_FGB_PATH
	{
	  AChunkMesg * mychunk = GetVA_AChunkMesg(thechunks,seed_chunk_index);
	  assert(mychunk != NULL);
	  if(seed_chunk_suffix == 0) {
	    seed_gang_key = 2*(mychunk->chunk_avx) + (mychunk->chunk_asx);
	  } else {
	    seed_gang_key = 2*(mychunk->chunk_bvx) + (mychunk->chunk_bsx);
	  }
	}

#endif // USE_FGB_PATH


#ifdef STORE_BRANCH_POINTS_AT_FRAGMENT
        {
          int known_bpt = get_bpt_vertex
            (frags,
             ( seed_chunk_suffix == 0 
               ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_avx
               : GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_bvx ),
             ( seed_chunk_suffix == 0 
               ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_asx
               : GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_bsx )
             );
          if(known_bpt > 0) continue;
          // For better or worse, do not re-compute a branch-point.
        }
#endif // STORE_BRANCH_POINTS_AT_FRAGMENT

        
#ifdef DEBUG32
        printf("before Find_Branch_Partners\n"
               "seed_chunk_index=" F_IID ",seed_chunk_suffix=%d\n",
               seed_chunk_index,seed_chunk_suffix);
#endif
        branch_olap_ct =
          Find_Branch_Partners
          ( /* input only */
           seed_chunk_index, seed_chunk_suffix,
           frags,
           edges,
           thechunks,
           /* clear and output */
           gangs,
           /* output only */
           &branch_mate_ct);
        
#ifdef DEBUG32
        printf("after Find_Branch_Partners\n"
               "branch_mate_ct=%d,branch_olap_ct=%d\n",
               branch_mate_ct,branch_olap_ct);
#endif
        
        /* At this point the local branch partners are known in terms of 
           their chunk-end keys and bp offset in the gangs data structure.
           Now read the gang and convert the chunk end keys to 
           chunk cid/suffix pairs. */
        
        EnableRangeVA_BranchOlap_t(theBranchOlaps,branch_olap_ct);
        
        {
          IntGang_ID kin_gang_key;
	  int ii;
          
          GetBranchOlap_t(theBranchOlaps,0)->cid = seed_chunk_index;
          GetBranchOlap_t(theBranchOlaps,0)->csx = seed_chunk_suffix;
#ifdef USE_FGB_PATH
          GetBranchOlap_t(theBranchOlaps,0)->vid = (seed_gang_key >> 1);
          GetBranchOlap_t(theBranchOlaps,0)->vsx = (seed_gang_key &  1);
#endif // USE_FGB_PATH

          GetBranchOlap_t(theBranchOlaps,0)->minhang =
            GetGangType(gangs,seed_gang_key)->min_offset;
          GetBranchOlap_t(theBranchOlaps,0)->maxhang =
            GetGangType(gangs,seed_gang_key)->max_offset;
          GetBranchOlap_t(theBranchOlaps,0)->use_me = TRUE;
          
          for(ii=1,kin_gang_key = GetGangType(gangs,seed_gang_key)->next; 
              (ii<branch_olap_ct) && (kin_gang_key != seed_gang_key); 
              ii++, kin_gang_key = GetGangType(gangs,kin_gang_key)->next) {
            
#ifndef USE_FGB_PATH
            GetBranchOlap_t(theBranchOlaps,ii)->cid = (kin_gang_key >> 1);
            GetBranchOlap_t(theBranchOlaps,ii)->csx = (kin_gang_key &  1);
#else // USE_FGB_PATH
	    {
	      const IntFragment_ID vid = (kin_gang_key >> 1);
	      const int vsx = (kin_gang_key &  1);
	      GetBranchOlap_t(theBranchOlaps,ii)->vid = vid;
	      GetBranchOlap_t(theBranchOlaps,ii)->vsx = vsx;
	      GetBranchOlap_t(theBranchOlaps,ii)->cid = get_cid_fragment(frags,vid);
	      GetBranchOlap_t(theBranchOlaps,ii)->csx = 
		(vsx ^ 
		 (get_o3p_fragment(frags,vid) <
		  get_o5p_fragment(frags,vid) ));
	      // if vsx == FALSE and o5p < 03p, then csx == FALSE.
	    }
#endif // USE_FGB_PATH
            GetBranchOlap_t(theBranchOlaps,ii)->minhang =
              GetGangType(gangs,kin_gang_key)->min_offset;
            GetBranchOlap_t(theBranchOlaps,ii)->maxhang =
              GetGangType(gangs,kin_gang_key)->max_offset;
            GetBranchOlap_t(theBranchOlaps,ii)->use_me = 
              GetGangType(gangs,kin_gang_key)->use_me;
          }
        }
        
        /* Clear the information to make ready for the next seed. */
        gang_clear(gangs, seed_gang_key);
        
        // histogram the number of mates and branch partners 
        // (branch_mate_ct,branch_olap_ct).
#ifdef DEBUG_WITH_HISTO
        add_to_histogram(histo_branch_mates,
                         branch_mate_ct,NULL);
        add_to_histogram(histo_branch_olaps,
                         branch_olap_ct,NULL);
        add_to_histogram(histo_branch_mate_times_olaps,
                         branch_mate_ct*branch_olap_ct,NULL);
#endif /*DEBUG_WITH_HISTO*/
#ifdef DEBUG31    
        // Display contents of array, then initialize, then call routine,
        // then display the results.
        fprintf(stdout, "\nBranch Overlaps for this chunk-end \n");
        { int iter; 
        for(iter = 0;  iter < branch_olap_ct;  iter ++)
          fprintf
            (stdout, "%d " F_IID " %d " F_IID " %d"
	     "%d %d\n",
             ((BranchOlap_t *)theBranchOlaps->Elements)[iter].use_me, 
             ((BranchOlap_t *)theBranchOlaps->Elements)[iter].cid, 
             ((BranchOlap_t *)theBranchOlaps->Elements)[iter].csx,
	     get_iid_fragment( frags,((BranchOlap_t *)theBranchOlaps->Elements)[iter].vid),
             ((BranchOlap_t *)theBranchOlaps->Elements)[iter].vsx,
             ((BranchOlap_t *)theBranchOlaps->Elements)[iter].minhang, 
             ((BranchOlap_t *)theBranchOlaps->Elements)[iter].maxhang);
        }
#endif /*DEBUG31*/

        branch_pt.place = 0; // Must be intialized.
        memset((char *)branch_pt.ct,0,NUM_BRANCH_COUNTS*sizeof(int16));
        // This just initializes the counts array to zero.
        
#ifndef DONT_FIND_BRANCH_POINTS
        Find_Branch_Points
          ( branch_olap_ct,
            (BranchOlap_t *)(theBranchOlaps->Elements),
            &branch_pt,
            thechunks, chunkseqs,
            frags, frag_store, use_consensus,
            fbpts);
#endif
        
        //#ifndef STORE_BRANCH_POINTS_AT_FRAGMENT
        /* Put the branch point information into the chunk end. */
        if(seed_chunk_suffix == 0) {
          GetVA_AChunkMesg(thechunks,seed_chunk_index)->a_branch_type = 
            ( branch_pt.place == 0 ? AS_NO_BPOINT : AS_INTO_UNIQUE);
          GetVA_AChunkMesg(thechunks,seed_chunk_index)->a_branch_point =
            branch_pt.place;
        } else {
          GetVA_AChunkMesg(thechunks,seed_chunk_index)->b_branch_type = 
            ( branch_pt.place == 0 ? AS_NO_BPOINT : AS_INTO_UNIQUE);
          GetVA_AChunkMesg(thechunks,seed_chunk_index)->b_branch_point =
            branch_pt.place;
        }
        //#else // STORE_BRANCH_POINTS_AT_FRAGMENT
#ifdef STORE_BRANCH_POINTS_AT_FRAGMENT
        set_bpt_vertex
          (frags,
           ( seed_chunk_suffix == 0 
             ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_avx
             : GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_bvx ),
           ( seed_chunk_suffix == 0 
             ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_asx
             : GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_bsx ),
           branch_pt.place);
#endif // STORE_BRANCH_POINTS_AT_FRAGMENT
        
        if( branch_pt.place != 0) {
          // Append to a file the known branch points.
          // "frag_iid,frag_suf,"
          // "branch_type,branch_point \n");
          fprintf
            ( fbpts,
              F_IID " %d "
              F_IID " %d %d %d\n",
              seed_chunk_index, seed_chunk_suffix,
              get_iid_fragment
              (frags,
               ( seed_chunk_suffix == 0 
                 ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_avx
                 : GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_bvx )),
              ( seed_chunk_suffix == 0 
                ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_asx
                : GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_bsx ),
              ( seed_chunk_suffix == 0 
                ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->a_branch_type
                : GetVA_AChunkMesg(thechunks,seed_chunk_index)->b_branch_type ),
              ( seed_chunk_suffix == 0 
                ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->a_branch_point
                : GetVA_AChunkMesg(thechunks,seed_chunk_index)->b_branch_point )
              );
        }
          
#ifdef DEBUG32
        if( branch_pt.place == 0) {
          fprintf(stdout, "\nNo Branch points found\n");
        } else {
          fprintf(stdout, "\nYes Branch points found\n");
          { int iter;
          for(iter = NUM_BRANCH_COUNTS - 1;  iter >= 0;  iter --)
            fprintf(stdout, " %3d", branch_pt.ct[iter]);
          }
          fprintf(stdout, "  : %3d\n", branch_pt.place);
        }
#endif /*DEBUG32*/
        

#ifdef DEBUG33
        printf
          ("seed_chunk_index=" F_IID ","
           "seed_chunk_suffix=%d,"
           "frag_iid=" F_IID ",frag_suf=%d,"
           "branch_type=%d,branch_point=%d\n",
           seed_chunk_index, seed_chunk_suffix,
           get_iid_fragment
           (frags,
            ( seed_chunk_suffix == 0 
              ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_avx
              : GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_bvx )),
           ( seed_chunk_suffix == 0 
             ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_asx
             : GetVA_AChunkMesg(thechunks,seed_chunk_index)->chunk_bsx ),
           ( seed_chunk_suffix == 0 
             ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->a_branch_type
             : GetVA_AChunkMesg(thechunks,seed_chunk_index)->b_branch_type ),
           ( seed_chunk_suffix == 0 
             ? GetVA_AChunkMesg(thechunks,seed_chunk_index)->a_branch_point
             : GetVA_AChunkMesg(thechunks,seed_chunk_index)->b_branch_point )
           );
#endif
      }
  }
  DeleteVA_BranchOlap_t(theBranchOlaps);
  DeleteVA_GangType(gangs);
#ifdef DEBUG_WITH_HISTO
  fprintf(fout,"\n\nHistogram of "
          "branch_mate_ct per chunk end.\n");
  print_histogram(fout,histo_branch_mates, 0, 1);
  
  fprintf(fout,"\n\nHistogram of "
          "branch_olap_ct per chunk end.\n");
  print_histogram(fout,histo_branch_olaps, 0, 1);
  
  fprintf(fout,"\n\nHistogram of "
          "branch_mate_ct*branch_olap_ct per chunk end.\n");
  print_histogram(fout,histo_branch_mate_times_olaps, 0, 1);
  
  free_histogram(histo_branch_mates);
  free_histogram(histo_branch_olaps);
  free_histogram(histo_branch_mate_times_olaps);
#endif /*DEBUG_WITH_HISTO*/

#ifdef  DEBUGGING
  // Dump_BP_List_Len_Distrib ();
#endif
//  fclose (BP_Log_File);
}



#endif // BRANCHPOINTS
