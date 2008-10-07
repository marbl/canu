
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

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "SUBDELREZ.h"

void Complement_Seq(char *a);

/* Complement the sequence a.  The operation does the
   complementation/reversal in place.  Calling it a second
   time on a given fragment restores it to its original state.  */

void Complement_Fragment_AS(InternalFragMesg *a);

/* Complement the sequence in fragment message a.  This include also
   revsersing the order of the quality values.  The operation does the
   complementation/reversal in place.  Calling it a second time on a
   given fragment restores it to its original state.                */

void Print_Overlap_AS(FILE *file, InternalFragMesg *a,
                                  InternalFragMesg *b, OverlapMesg *align);

/*  Print an ASCII representation of the alignment between fragments a and
    b encoded in align to the file "file".

    Within the file AS_ALN_qvaligner.c the defined constant PRINT_WIDTH
    (set to 50) controls the number of columns per line in the display
    of the alignment.                                                     */


typedef struct {
  int begpos;
  int endpos;
  int length;
  int diffs;
  int comp;
  int *trace;
} Overlap;

void Print_Overlap(FILE *file, char *aseq, char *bseq, Overlap *align);

/*  Print an ASCII representation of the alignment between fragments a and
    b encoded in align to the file "file".

    Within the file AS_ALN_qvaligner.c the defined constant PRINT_WIDTH
    (set to 50) controls the number of columns per line in the display
    of the alignment.                                                     */

Overlap *Copy_Overlap(Overlap *ovl);

/* Make a copy of overlap ovl, allocating memory for the copy. */

#define AS_ANALYZE_ALL           0
#define AS_ANALYZE_DELETES       1
#define AS_ANALYZE_INSERTS       2
#define AS_ANALYZE_SUBSTITUTIONS 3


void Analyze_Affine_Overlap_AS(InternalFragMesg *a, InternalFragMesg *b,
                               OverlapMesg *align, int amode,
                               int *alen, int *blen,
                               int *del, int *sub, int *ins,
                               int *affdel, int *affins,
                               int *blockdel, int *blockins, int blocksize,
                               int *biggestBlock);
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


Overlap *DP_Compare(char *aseq, char *bseq,
                    int beg, int end, int opposite,
                    double erate, double thresh, int minlen,
                    CompareOptions what);

/* As below for DP_Compare_AS(), but with different interface; this (DP_Compare()) is now the core routine while DP_Compare_AS() is a wrapper */

OverlapMesg *LD_DP_Compare_AS(InternalFragMesg *a, InternalFragMesg *b,
                           int beg, int end, int opposite,
                           double erate, double thresh, int minlen,
                           CompareOptions what, int *where, double my_ld_ratio);

/* Temporarily modify ld_ratio, run DP_Compare_AS, then change ld_ratio back */

OverlapMesg *DP_Compare_AS(InternalFragMesg *a, InternalFragMesg *b,
                           int beg, int end, int opposite,
                           double erate, double thresh, int minlen,
                           CompareOptions what, int *where);

/* Given fragments a and b, find the best overlap between them subject
   to the following parameters/thresholds.  The overlap must start on
   one of the diagonals of the d.p. matrix in the interval [beg,end].
   For example if one gives the interval [-10,20], then the overlap
   either has less than the first 20bp of a unaligned or less than the
   first 10bp of b unaligned.  If the boolean variable `opposite' is nonzero
   then the fragments are to be aligned in the opposite orientation.  One
   is assuming an error rate of `erate' in the sequences, and is guaranteed
   to find only alignments for which the number of differences d in each
   prefix of length n is such that
           Sum_k=d^n (n choose k) erate^k (1-erate)^(n-k) < thresh.
   One should note carefully, that alignments not satisfying this property
   may be found, the point is that ones that don't may be missed.
   In addition, the alignment must involve at least `minlen' symbols of the
   prefix-sequence in the overlap.  The option `what' specifies what kind
   of comparison is to be performed as follows:

   AS_FIND_OVERLAP:
      Just find a good alignment to the boundary of the d.p. matrix.  From
      this extrapolate a rough overlap relationship without an alignment
      and return the result (if there is one).
   AS_FIND_ALIGN:
      For this option, further go to the trouble of computing the alignment
      and store it in the overlap message.
   AS_FIND_ALIGN_NO_TRACE:
      For this option, further go to the trouble of computing the alignment
      and store it in the overlap message.  Don't compute or return the
      encoded alignment.  A Hack by SAK.
   AS_FIND_AFFINE_ALIGN
      For this option, find the best alignment using an affine gap cost
      measure.  Substitutions cost SUBCOST and gaps of length len cost
      GAPCOST + len, where SUBCOST and GAPCOST are defined constants within
      AS_ALN_dpaligner.c
   AS_FIND_QVALIGN:
      NOT YET IMPLEMENTED.  Will ultimately use quality values to compute
      the best possible alignment.

   As for all other routines, the space for the overlap message is owned
   by the routine, and must be copied if it is to be retained beyond the
   given call.  The routine also returns in the integer pointed at by
   where, the starting diagonal of the alignment it finds.  While technically
   this can be inferred from the overlap message, it is easier to take this
   position directly and see if the returned alignment actually starts
   in the interval [beg,end].  If it does not, then it calls the initial
   beg,end range into question.                                            */




/*****************************************************************/
/* New prototypes for bubble overlap detectors                   */
/*****************************************************************/

OverlapMesg *AS_ALN_affine_overlap(InternalFragMesg *a, InternalFragMesg *b,
                           int beg, int end, int opposite,
                           double erate, double thresh, int minlen,
                           CompareOptions what, int *where) ;

/* Bubble smoothing overlap detector based on
   affine dp_compare with subsequent filtering.

   Usage/arguments as for DP_Compare_AS, except that some of the parameters
   are hijacked:
   - "erate" is used as a filter after the fact, but not passed to DP_C*_AS
   - "what" is overridden
   - "minlen" is passed to DP_Compare_AS, but it is also used as a filter
     to require that minlen matches were found

   Procedure:
   - call DP_Compare_AS, using maximal erate and what=AS_FIND_ALIGN;
            (should use affine alignment option when Gene checks it in).
   - use Analyze_Affine_Overlap_AS() to evaluate the resulting alignment.
   - if the error rate in the alignment, adjusted for an affine scoring scheme,
     is better than user's erate, return the OverlapMesg*; else return NULL;
     optionally, also test number of large indels.

   Assumptions/caveats:
   - as with DP_Compare_AS, the returned message must be copied if it is to
   be retained.
   - this version does not test for sanity with respect to placement of gaps;
     in principle, a gap of several hundred bases at the end (most likely
     indicating a true branchpoint) would be accepted if proposed; this
     seems to be safe enough since DP_Compare_AS doesn't seem to find overlaps
     above about 12% simple (non-affine) error rate.  However, more minor
     versions of this could cause (false) overlaps of shallow branchpoints
     (true branchpoints occurring near the ends of fragments)

       Branchpoint:

               .........+++++++          ("." matches; +:# mismatch)
	       .........#######

          could be treated as overlapping with an affine gap:

               .........+++++++
	       .........-------######    ("-" a gap)

     Equally, a (short) bad fragment end has a better change of being
     overlapped if affine gaps allow it to find the best match within a
     modest window of uncertainty:

         Insufficient trimming of the second fragment:

               .....ACAGTAGACGAGATAGGATAGATAGAGTAGACAGATAGTTGACTAAC
	            ||||||||||||||||||||||
               .....ACAGTAGACGAGATAGGATAGACAGTTA

	 could be interpreted as:

               .....ACAGTAGACGAGATAGGATAGATAGAGTAGACAGATAGTTGACTAAC
	            ||||||||||||||||||||||         ||| ||
               .....ACAGTAGACGAGATAGGATAGA---------CAGTTA

*/
/*end comments for AS_ALN_affine_overlap */



#include "CA_ALN_local.h"
int *AS_Local_Trace(Local_Overlap *local_overlap, char *aseq, char *bseq);

/* Create a trace to be interpreted as with DP_Compare_AS, but based
   on a Local_Overlap record.  A Local_Segment within the overlap
   will be aligned using OKNAlign(), generating a subtrace.  Subtraces,
   with their indices appropriately adjusted, will be spliced together
   by an encoding of the gaps between segments; for now, we'll simply insert
   gaps as follows:

      A "gap" with x bases in A and y bases in B will become a section
      of the alignment x+y positions long, with the A fragment first
      and the B fragment second:

               AAAAAAAAAA--------------
               ----------BBBBBBBBBBBBBB

   Obviously, a more compact treatment is possible!

   Assumptions: both sequences should be in the forward orientation
   and all segments are forward.

*/

OverlapMesg *Local_Overlap_AS(InternalFragMesg *a, InternalFragMesg *b,
                           int beg, int end, int opposite,
                           double erate, double thresh, int minlen,
                           CompareOptions what, int *where);


/* Given fragments a and b, find the best overlap between them using
   local alignment code.  The purpose is to permit "bubbles" in the fragment
   graph due to multi-base polymorphisms to be smoothed out, such that a
   unitig can be constructed.

   The method relies on Myers' local overlap code.  MORE DETAILS????

   The local overlap code returns the best overlap available as a chain of
   matching segments.

   The current function is concerned with three things:

   1) Providing a wrapper around the local alignment code so that users
   can carry out operations on the level of the "message" (IFM, OVL, etc).

   2) Determine whether the best possible overlap is good enough to constitute
   an overlap for the purposes of unitigging.  A preliminary set of criteria
   are given below.

   3) Combine bits of existing code to provide a full alignment specification
   for the local overlapper; one can debate whether large gaps should be
   encoded in an edit trace, but for the purposes of integration with the
   existing code base, the current function will aim to specify an
   alignment in terms of an ahang, a bhang and a trace, which will get
   converted into the core elements of an overlap message just as in
   DP_Compare_AS.

   PRELIMINARY DECISION CRITERIA:

     - if an overall %-match of less than erate with length MinLen[=40]
       exists, then the overlap should be allowed regardless of the
       following additional criteria
     - at least one segment of at least MinLen bases
     - if we treat an affine event as a single error (and the length involved
       to likewise be a single base) for the purpose of calculating an
       error rate, the error rate should be no more than erate.
     - no more than MaxGaps [=2?] gaps in the overlap chain;
       this also entails no more than MaxGaps+1 matching segments

     Additional criteria that we might want to use:

     - no more than MaxEndGap [=20?] unaligned bases at either end; this
       is to prevent a classic branchpoint out of a repeat into distinct
       regions from being overlapped; on the other hand, the requirement
       that all fragments in a bubble go back together might mean that
       this filter is unnecessary.
     - only one end may end in a gap >= BigGapSize; having both ends
       be in polymorphic regions should be quite unlikely.

     One thing to bear in mind: it is not the intent of this function
     (or bubble smoothing in general) to handle fragments which are
     generally noisy or have ends which were incompletely quality
     trimmed.  At least the latter is a spur and handled by other
     means.


   AS_FIND_LOCAL_OVERLAP:
      Just find a good alignment to the boundary of the d.p. matrix.  From
      this extrapolate a rough overlap relationship without an alignment
      and return the result (if there is one).
   AS_FIND_LOCAL_ALIGN:
      For this option, further go to the trouble of computing the alignment
      and store it in the overlap message.
   AS_FIND_LOCAL_ALIGN_NO_TRACE
      A hack by SAK...don't compute delta encoding.

   As with DP_Compare_AS, the resulting overlap message occupies
   memory which is persistently owned (reused) by the function; it
   must be copied if it is to be retained beyond the given call.

*/

/* Given fragments a and b, find the best chain of local alignments, trimmed
   to non-overlapping.  Related to Local_Overlap_AS, but returns really bad
   alignments if that's the best there is to be had.
*/
OverlapMesg *BoxFill_AS(InternalFragMesg *a, InternalFragMesg *b,
			int beg, int end, int opposite,
			double erate, double thresh, int minlen,
			CompareOptions what, int *where);



void PrintAlign(FILE *file, int prefix, int suffix,
                       char *a, char *b, int *trace);

/* Print an alignment to file between a and b given in trace (unpacked).
   Prefix gives the length of the initial prefix of a that is unaligned.  */

int *AS_ALN_OKNAlign(char *a, int alen, char *b, int blen, int *spnt, int diff);

/* O(kn) identity-based alignment algorithm.  Find alignment between
   a and b (of lengths alen and blen), that begins at finishing
   boundary position *spnt.  Return at *spnt the diagonal at which the
   alignment starts.                                                   */

int *AS_ALN_OKNAffine(char *a, int alen, char *b, int blen,
                      int *bpnt, int *epnt, int diff);

/* O(kn) affine gap cost alignment algorithm.  Find best alignment between
   a and b (of lengths alen and blen), within a band of width 2*diff centered
   on the diagonal containing finishing boundary position *spnt.  Return at
   *spnt the diagonal at which the alignment starts.  A quick implementation
   that is space inefficient, takes O(kn) space as opposed to the O(n) that
   is possible.                                                             */


// From AS_ALN_pieceOlap.c

/* fix_overlapping_pieces():
   Handle Local_Overlap pieces (local segments) which overlap,
   by trimming them back until they abut,
   such that the number of mismatches is minimized */

void fix_overlapping_pieces(char *aseq, char *bseq,
			    Local_Overlap *O,int piece0, int piece1);


/*  Compute_Olap_Version converts and OverlapMesg into the orientation and hangs that would be reported as an "olap" (e.g. by dump-olap-store) */

void Compute_Olap_Version(InternalFragMesg* a,InternalFragMesg *b,OverlapMesg *O,int *ahang,int *bhang, char *ori);

/*  AS_ALN_clean_up_trace removes leading and trailing gaps */

void AS_ALN_clean_up_trace(int *trace,int alen, int blen,int *spos,int *epos);

#endif


