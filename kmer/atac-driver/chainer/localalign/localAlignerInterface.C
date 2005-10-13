// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Clark Mobarry
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "GF_ALN_global.h"
#include "GF_ALN_local.h"
#include "GF_ALN_aligners.h"

#define DEBUG_LOCALOVL 1

static Local_Overlap *desc = NULL;

void syntenicSegments
(
 FILE * outfile,
 char const * const Aseq, int const Astart, int const Astop,
 char const * const Bseq, int const Bstart, int const Bstop,
 float const erate
 )
{
  // All outputs are to stdout.

  desc = NULL; // In case an early exit happens!
  
  //fprintf(stderr,"Hello Alice1\n");

  /* Key data types ("Local_Segment" and "Local_Overlap" are defined in
     "CA_ALN_local.h" -- of which there is a copy in AS_ALN!)*/

  int bpReportableMatchInA = 0;
  int bpReportableMatchInB = 0;
  int bpChainMatchInA = 0;
  int bpChainMatchInB = 0;

  // Step 1: get local segments:
  char const * const Ausable = Aseq + Astart;
  char const * const Busable = Bseq + Bstart;
  assert(Astop >= Astart);
  assert(Bstop >= Bstart);
  int const alen = Astop - Astart;
  int const blen = Bstop - Bstart;
  /* bool */ int const reverse_bstring = 0;
  int const minlen_reportable_match = 16; /* minimum length of a reportable match */
  int const minlen_overlap = 20;
  int NumSegs = 0; /* number of local matches returned */

  Local_Segment *local_results = NULL;
  Local_Overlap *Ov = NULL;

  assert(reverse_bstring == 0); // CMM
  
  //fprintf(stderr,"Hello Alice2\n");
  
  local_results = Find_Local_Segments
    (
     Ausable /* sequence A */,
     alen,
     Busable, 
     blen,
     (reverse_bstring ? LOCAL_REVR : LOCAL_FORW),
     /* whether to compute a forward search , reverse, or both */
     minlen_reportable_match /* minimum length of a reportable match */, 
     erate /* maximum error for a match to be returned */, 
     &NumSegs /* number of local matches returned */);
  
#if DEBUG_LOCALOVL > 2
  if(NumSegs==0) {
    fprintf(outfile, "#LA_FLS no matches found\n");
  } else {
    fprintf(outfile, "#LA_FLS %d matches found\n",NumSegs);
  }
#endif
  if(NumSegs==0) {
    
    return;
  }
  
#if DEBUG_LOCALOVL > 2
  {
    int i;
    printf("@@#n\n@@#x 0 %d\n@@#y 0 %d\n@@#m 1\n",alen,blen);
    for(i=0;i<NumSegs;i++){
      bpReportableMatchInA += local_results[i].aepos - local_results[i].abpos;
      bpReportableMatchInB += local_results[i].bepos - local_results[i].bbpos;
      printf("@@#k Segment %d: (%d,%d)[----------](%d,%d), error %e\n"
             "@@%d %d\n@@%d %d\n"
             ,i
             ,local_results[i].abpos
             ,local_results[i].bbpos
             ,local_results[i].aepos
             ,local_results[i].bepos
             ,local_results[i].error
             ,local_results[i].abpos
             ,local_results[i].bbpos
             ,local_results[i].aepos
             ,local_results[i].bepos
             );
    }
    fprintf(outfile,"#NOISE2  bpReportableMatchInA= %d bpReportableMatchInB= %d\n",
            bpReportableMatchInA, bpReportableMatchInB);

  }
#endif


  // Step 2: get a chain of local segments:

  Ov = Find_Local_Overlap
    (
     alen              /* length of sequence A */,
     blen              /* length of sequence B */,
     reverse_bstring   /* comp==0 -> fwd orientation */,
     0                 /* nextbest==0 -> find best overlap*/,
     local_results     /* the input set of local segments */,
     NumSegs           /* number of input local segments */,
     minlen_overlap-6  /* shortest "overlap" to report" */, 
     1.                /* fraction of overlap not in a match -- 
                          needs to be large to allow
                          substantial mismatches */
     );

#if DEBUG_LOCALOVL > 2
  if(Ov == NULL) {
    fprintf(outfile, "#LA_FLO no decent chain found\n");
  } else {
    fprintf(outfile, "#LA_FLO a decent chain found using num_pieces= %d from NumSegs= %d\n",
	    Ov->num_pieces, NumSegs );
  }
#endif

  if(Ov == NULL){
    return;
  }
  
#if DEBUG_LOCALOVL > 2
  if(NULL != Ov) {
    int i;
    printf("#@@#m 11\n");
    for(i=0;i<Ov->num_pieces;i++){
      bpChainMatchInA += Ov->chain[i].piece.aepos - Ov->chain[i].piece.abpos;
      bpChainMatchInB += Ov->chain[i].piece.bepos - Ov->chain[i].piece.bbpos;
      
      printf("#@@#k Overlap %d: (%d,%d)[----------](%d,%d), error %e\n"
             "#@@%d %d\n@@%d %d\n",i,
             Ov->chain[i].piece.abpos,
             Ov->chain[i].piece.bbpos,
             Ov->chain[i].piece.aepos,
             Ov->chain[i].piece.bepos,
             Ov->chain[i].piece.error,
             Ov->chain[i].piece.abpos,
             Ov->chain[i].piece.bbpos,
             Ov->chain[i].piece.aepos,
             Ov->chain[i].piece.bepos
             );
    }
    fprintf(outfile,"#NOISE3  bpChainMatchInA= %d bpChainMatchInB= %d\n",
            bpChainMatchInA, bpChainMatchInB);
  }
#endif

#if 0
  fprintf(outfile,
          "#NOISE4 bpChainMatchInA/bpReportableMatchInA= %.3f bpChainMatchInB/bpReportableMatchInB= %.3f\n",
          (float)bpChainMatchInA/((float)bpReportableMatchInA + 1e-6f),
          (float)bpChainMatchInB/((float)bpReportableMatchInB + 1e-6f));
#endif
          
  // Step 3 (optional): 
  //
  // a) fix the chain of segments so that the segments don't overlap.
  // It must be a 1-1 mapping. (can either trim or delete segments--or
  // leave them completely alone)
  //
  // b) construct an alignment "trace" 
  //
  // The "trace" is the standard "AS" encoding of an alignment.

  if(Ov != NULL) {

    // coordinate munge between Gene's local aligner and
    // DP_Compare()-related routines coordinates from Find_Local
    // routines will be one off from those expected by the trace
    // routines, so adjust them!

    for(int i=0;i<=Ov->num_pieces;i++){
      if(i<Ov->num_pieces){
        Ov->chain[i].piece.abpos++;
        Ov->chain[i].piece.bbpos++;
        Ov->chain[i].piece.aepos++;
        Ov->chain[i].piece.bepos++;
      }
    }

    //  AS_Local_Trace assumes string pointer one before start of string!

    if(Ov != NULL) {
      int *trace = AS_Local_Trace(Ov,Ausable-1,Busable-1);
      if(trace == NULL) {
        fprintf(stderr,"EXCEPTION Ov=%p trace=%p\n", Ov, trace);
      }
    }

    for(int i=0;i<=Ov->num_pieces;i++){
      if(i<Ov->num_pieces){
        Ov->chain[i].piece.abpos--;
        Ov->chain[i].piece.bbpos--;
        Ov->chain[i].piece.aepos--;
        Ov->chain[i].piece.bepos--;
      }
    }
  }

  if(Ov != NULL) { Ov->next = 0;}
  desc = Ov;
}



static void Print_Local_Overlap2
( FILE *outfile,
  Local_Overlap *desc,
  char const * selfid,
  char const * parentid,
  char const * a_uid,
  char const * b_uid,
  int const a_scaf_pos, // Lowest position in string "A" for the alignment.
  int const b_scaf_pos, // Lowest position in string "B" for the alignment.
  int const alen,       // Length of the alignment box for string "A".
  int const blen,       // Length of the alignment box for string "B".
  int const sequence_reversed_match, /* bool */
  int const a_lft_seed_len,
  int const b_lft_seed_len,
  int const a_rht_seed_len,
  int const b_rht_seed_len
  )
{ 
  /* Adapted from Print_Local_Overlap(). */
  assert(NULL != outfile);
  assert(NULL != desc);
  {
    Local_Chain *chain = desc->chain;
    assert(NULL != desc->chain);

#if DEBUG_LOCALOVL > 2
    {
      fprintf(outfile,"#LA %d segments ", desc->num_pieces);
      if(desc->num_pieces > 1){
        fprintf(stderr,"with %d unaligned bp's between them",
                desc->score);
      }
      fprintf(outfile," diffs = %d indif = %d score = %d overlap = %d\n",
              desc->diffs, desc->indif, desc->score, desc->length);
    }
#endif
    {
      int i;
      int const ifirst=0;
      int const ilast=desc->num_pieces -1;
      for( i = ifirst; i <= ilast; i++) {
        const Local_Segment * const seg =  & (chain[i].piece);
        assert(NULL != seg);
        assert(!chain[i].reversed);
#if DEBUG_LOCALOVL > 2
        fprintf(outfile,"#LA [%d,%d] x [%d,%d] @%.2f\n",
                seg->abpos,seg->aepos,
                seg->bbpos,seg->bepos, 100.*seg->error);
#endif
        assert(seg->aepos > seg->abpos);
        assert(seg->bepos > seg->bbpos);

        if( i==ifirst ) {
          /* bool */ int const overlaps_lft = (seg->abpos < a_lft_seed_len);
          /* bool */ int const overlaps_bot = (seg->bbpos < b_lft_seed_len);
          char aflag = 'A';  // Abuts the corner.
          char bflag = 'A';  // Abuts the corner.
          if(!overlaps_lft ) {
            aflag = 'N'; // The A-string first piece does not overlap the left seed hit.
          } else if(seg->abpos != 0) {
            aflag = 'O'; // The A-string first piece does not abut the left corner.
          }
          if(!overlaps_bot ) {
            bflag = 'N'; // The B-string first piece does not overlap the left seed hit.
          } else if(seg->bbpos != 0) {
            bflag = 'O'; // The B-string first piece does not abut the left corner.
          }
#if DEBUG_LOCALOVL > 2
          fprintf(outfile,"#LA PIECE F %c %c %d %d\n", aflag, bflag, seg->abpos, seg->bbpos);
#endif
        }

        fprintf(outfile,
                "M %s %s_%d %s %s %d %d %d %s %d %d %d > %.2f\n",
                "l", selfid, i, parentid,
                a_uid, (a_scaf_pos + seg->abpos), (seg->aepos - seg->abpos), 1,
                b_uid,
                (sequence_reversed_match ? (b_scaf_pos + blen - seg->bepos) : (b_scaf_pos + seg->bbpos)),
                (seg->bepos - seg->bbpos),
                (sequence_reversed_match ? -1 : 1),
                100.*seg->error
                );
        
        if( i==ilast ) {
          /* bool */ int const overlaps_rht = (alen - seg->aepos < a_rht_seed_len);
          /* bool */ int const overlaps_top = (blen - seg->bepos < b_rht_seed_len);
          char aflag = 'A';  // Abuts the corner.
          char bflag = 'A';  // Abuts the corner.
          if(!overlaps_rht ) {
            aflag = 'N'; // The A-string last piece does not overlap the right seed.
          } else if(seg->aepos != alen) {
            aflag = 'O'; // The A-string last piece does not abut the left corner.
          }
          if(!overlaps_top ) {
            bflag = 'N'; // The B-string last piece does not overlap the right seed.
          } else if(seg->bepos != blen) {
            bflag = 'O'; // The B-string last piece does not abut the right corner.
          }
#if DEBUG_LOCALOVL > 2
          fprintf(outfile,"#LA PIECE L %c %c %d %d\n",  aflag, bflag, alen - seg->aepos, blen - seg->bepos);
#endif
        }
      }
    }
  }
}


int iterate_Local_Overlap
(
 int& seg_abpos,
 int& seg_bbpos,
 int& seg_alen,
 int& seg_blen,
 float& seg_error
 )
{
  // Using global variable "desc".
  
  /* Adapted from Print_Local_Overlap(). */

  if( NULL != desc ) {
    Local_Chain *chain = desc->chain;
    assert(NULL != desc->chain);
    for(; 0 <= desc->next && desc->next < desc->num_pieces; ) {
      int the_piece = (desc->next)++;
      
      const Local_Segment * const seg =  & (chain[the_piece].piece);
      assert(NULL != seg);
      assert(!chain[the_piece].reversed);
      
      // Set the return data:
      seg_abpos = seg->abpos;
      seg_alen = seg->aepos - seg->abpos;
      seg_bbpos = seg->bbpos;
      seg_blen = seg->bepos - seg->bbpos;
      seg_error = seg->error;

      // Skip over the "deleted in-place" segments.
      if((seg->aepos <= seg->abpos)&&(seg->bepos <= seg->bbpos)) { continue;}
      
      //assert(seg->aepos > seg->abpos);
      //assert(seg->bepos > seg->bbpos);
      
      return 1; // the data is valid
    }
  }
  return 0; // the data is invalid
}

int Print_Local_Overlap4(
 FILE *outfile,
 char const * selfid,
 char const * parentid,
 char const * a_uid,
 char const * b_uid,
 int const a_scaf_pos, // Lowest position in string "A" for the alignment.
 int const b_scaf_pos, // Lowest position in string "B" for the alignment.
 int const alen,       // Length of the alignment box for string "A".
 int const blen,       // Length of the alignment box for string "B".
 int const sequence_reversed_match, /* bool */
 int const a_lft_seed_len,
 int const b_lft_seed_len,
 int const a_rht_seed_len,
 int const b_rht_seed_len
 )
{
  // returns 1 for valid output
  // returns 0 for invalid output
  
#if DEBUG_LOCALOVL > 2
    {
      fprintf(outfile,"#LA %d segments ", desc->num_pieces);
      if(desc->num_pieces > 1){
        fprintf(stderr,"with %d unaligned bp's between them",
                desc->score);
      }
      fprintf(outfile," diffs = %d indif = %d score = %d overlap = %d\n",
              desc->diffs, desc->indif, desc->score, desc->length);
    }
#endif

    {
      int valid = 1;
      int the_piece = 0;
      while( valid == 1) {
        int seg_abpos;
        int seg_alen;
        int seg_bbpos;
        int seg_blen;
        float seg_error;
        
        valid = iterate_Local_Overlap
          (
           seg_abpos,
           seg_bbpos,
           seg_alen,
           seg_blen,
           seg_error
           );
        
        if(valid){

          fprintf(outfile,
                  "M %s %s_%d %s %s %d %d %d %s %d %d %d > %.2f\n",
                  "l", selfid, the_piece, parentid,
                  a_uid,
                  (a_scaf_pos + seg_abpos), seg_alen, 1,
                  b_uid,
                  (b_scaf_pos + (sequence_reversed_match ? (blen - seg_blen - seg_bbpos) : seg_bbpos)),
                  seg_blen,
                  (sequence_reversed_match ? -1 : 1),
                  100.*seg_error
                  );
          
        }      
        the_piece ++;
      }
    }
    return 0;
}

    
