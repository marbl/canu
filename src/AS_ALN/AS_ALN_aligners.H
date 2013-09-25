
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

#ifndef AS_ALN_ALIGNERS_H
#define AS_ALN_ALIGNERS_H

static const char *rcsid_AS_ALN_ALIGNERS_H = "$Id$";

#include "AS_global.H"
#include "AS_MSG_pmesg.H"
#include "SUBDELREZ.H"

typedef struct {
  uint32           aifrag;
  uint32           bifrag;
  int32            ahg;
  int32            bhg;
  PairOrient       orientation;
  OverlapType      overlap_type;
  float            quality;
  int32            *alignment_trace;
} ALNoverlapFull;  //  Former 'OverlapMesg'

typedef struct {
  int begpos;
  int endpos;
  int length;
  int diffs;
  int comp;
  int *trace;
} ALNoverlap;  //  Former 'Overlap'


void
PrintALNoverlap(char *message,
                char *aseq,
                char *bseq,
                ALNoverlap *align);

ALNoverlap *CopyALNoverlap(ALNoverlap *ovl);


#define AS_ANALYZE_ALL           0
#define AS_ANALYZE_DELETES       1
#define AS_ANALYZE_INSERTS       2
#define AS_ANALYZE_SUBSTITUTIONS 3


/* Just like Analyze_Overlap_AS with addition of evaluation of affine model
   Analyze the overlap between fragments a and b.
   Returns:  (previous versions returned a list of positions at which these
   errors occured, but that was buggy and not used)
     alen - # of a symbols in overlap,
     blen - # of b symbols in overlap,
     del  - # of unaligned symbols in a,
     sub  - # of substitutions,
     ins  - # of unaligned symbols in b,
     affdel - # of runs of unaligned symbols in a,
     affins - # of runs of unaligned symbols in b,
     blockdel - # of runs of size > blocksize of unaligned symbols in a,
     blockins - # of runs of size > blocksize of unaligned symbols in b,
     const blocksize - min length of an indel to count as a block.
     optional biggestBlock - the length of the largest mismatch, OK to give it NULL here

*/
void Analyze_Affine_ALNoverlapFull(InternalFragMesg *a,
                                   InternalFragMesg *b,
                                   ALNoverlapFull *align,
                                   int amode, int *alen, int *blen,
                                   int *del, int *sub, int *ins,
                                   int *affdel, int *affins,
                                   int *blockdel, int *blockins, int blocksize,
                                   int *biggestBlock);



typedef enum {
  AS_FIND_OVERLAP,
  AS_FIND_ALIGN,
  AS_FIND_ALIGN_NO_TRACE,  // Same as FIND_ALIGN, but no deltas are returned... SAK
  AS_FIND_QVALIGN,
  AS_FIND_AFFINE_OVERLAP,  //series of AFFINE equivalents to basics!
  AS_FIND_AFFINE_ALIGN,
  AS_FIND_AFFINE_ALIGN_NO_TRACE,
  AS_FIND_LOCAL_OVERLAP,  //series of LOCAL equivalents to basics!
  AS_FIND_LOCAL_ALIGN,
  AS_FIND_LOCAL_ALIGN_NO_TRACE
} CompareOptions;

typedef ALNoverlapFull *(AS_ALN_Aligner_AS)(InternalFragMesg *a, InternalFragMesg *b,
                                            int beg, int end,
                                            int opposite,
                                            double erate, double thresh, int minlen,
                                            CompareOptions what, int *where);

typedef ALNoverlap *(AS_ALN_Aligner)(char *aseq, char *bseq,
                                     int beg, int end,
                                     int ahang, int bhang,
                                     int opposite,
                                     double erate, double thresh, int minlen,
                                     CompareOptions what);

typedef int (AS_ALN_OverlapScorer)(ALNoverlap *O,
                                   int      expected_length,
                                   int      ahang_input,
                                   int      bhang_input,
                                   double   maxerate,
                                   double  *lScore_out,
                                   double  *aScore_out,
                                   double  *bScore_out);

AS_ALN_Aligner_AS  DP_Compare_AS;
AS_ALN_Aligner_AS  Affine_Overlap_AS;
AS_ALN_Aligner_AS  Local_Overlap_AS;

AS_ALN_Aligner     DP_Compare;
AS_ALN_Aligner     Local_Overlap_AS_forCNS;
AS_ALN_Aligner     Affine_Overlap_AS_forCNS;
AS_ALN_Aligner     Optimal_Overlap_AS_forCNS;

AS_ALN_OverlapScorer ScoreOverlap;

#include "CA_ALN_local.H"
int *AS_Local_Trace(Local_Overlap *local_overlap, char *aseq, char *bseq);


//  Check that the supplied trace makes sense
int IsInvalidAlign(int prefix, int suffix,
                   char *a, char *b, int *trace);

/* Print an alignment to file between a and b given in trace (unpacked).
   Prefix gives the length of the initial prefix of a that is unaligned.  */
void PrintAlign(FILE *file, int prefix, int suffix, char *a, char *b, int *trace);


/* O(kn) identity-based alignment algorithm.  Find alignment between
   a and b (of lengths alen and blen), that begins at finishing
   boundary position *spnt.  Return at *spnt the diagonal at which the
   alignment starts.                                                   */
int32 *AS_ALN_OKNAlign(char *a, int64 alen, char *b, int64 blen, int32 *spnt, int64 diff);


/* O(kn) affine gap cost alignment algorithm.  Find best alignment between
   a and b (of lengths alen and blen), within a band of width 2*diff centered
   on the diagonal containing finishing boundary position *spnt.  Return at
   *spnt the diagonal at which the alignment starts.  A quick implementation
   that is space inefficient, takes O(kn) space as opposed to the O(n) that
   is possible.                                                             */
int32 *AS_ALN_OKNAffine(char *a, int64 alen, char *b, int64 blen, int32 *bpnt, int32 *epnt, int64 diff); 


/* fix_overlapping_pieces():
   Handle Local_Overlap pieces (local segments) which overlap,
   by trimming them back until they abut,
   such that the number of mismatches is minimized */
void fix_overlapping_pieces(char *aseq, char *bseq, Local_Overlap *O,int piece0, int piece1);


/*  AS_ALN_clean_up_trace removes leading and trailing gaps */
void AS_ALN_clean_up_trace(int *trace,int alen, int blen,int *spos,int *epos);

#endif


