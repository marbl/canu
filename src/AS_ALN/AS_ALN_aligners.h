
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

int *Unpack_Alignment_AS(OverlapMesg *align);

/* Convert overlap message align's delta into an array of ints that more
   directly encode the alignment.  For an unpacked trace < i1, i2, ... in, 0>
   a negative number j indicates that a dash should be placed before A[-j]
   and a positive number k indicates that a dash should be placed before
   B[k], where A and B are the two sequences of the overlap.  These indels
   occur in order along the alignment.

   A pointer to an array containing the unpacked trace is returned.  This
   array is owned by the routine and is reused by it with each subsequent
   call.  If the unpacked trace is needed beyond a subsequent call, the
   caller must copy its contents to a memory area they have allocated.   */

signed char *Pack_Alignment_AS(int *trace, int prefix);

/*  Produce an overlap delta for an unpacked trace between two sequences,
    say A and B, where the first prefix symbols of A are unaligned with B
    (prefix must be > 0).

    A pointer to the delta is returned.  This delta is owned by the routine
    and reused by it with each subsequence call.                           */

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

#if 0

//  2005-oct-19, BPW, dead code.

int *Analyze_Overlap_AS(InternalFragMesg *a, InternalFragMesg *b,
                        OverlapMesg *align, int amode,
                        int *alen, int *blen, int *del, int *sub, int *ins);

/* Analyze the overlap between fragments a and b.
   Returns a -1 terminated list of the positions in the a sequences at
   which an error of type "amode" occurs, as well as:
     alen - # of a symbols in overlap,
     blen - # of b symbols in overlap,
     del  - # of unaligned symbols in a,
     sub  - # of substitutions,
     ins  - # of unaligned symbols in b.                                 */

#endif

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




OverlapMesg *QV_ReAligner_AS(InternalFragMesg *a, InternalFragMesg *b,
                             OverlapMesg *align);

/* Given fragments a and b and overlap align between them (computed by an
   identity-based method), refine the overlap using quality values for the
   fragments and return a pointer to a buffered overlap message containing
   the new result.  The realignment takes place within a band of width
   BAND_WIDTH (defined constant in AS_ALN_qvaligner.c).  Like other routines
   in the module, the space for the overlap message is owned by the routine
   and is reused by it on each subsequence call.  One must make a copy of
   a result if the object is to be retained beyond a single call.

   One must be aware that both Unpack_Alignment_AS and Pack_Alignment_AS
   are used by this routine, so the most recent returns from these routines
   (if called by the user) are destroyed by calling this routine.

   THIS ROUTINE STILL NEEDS TO HAVE QV-BASED SCORING ADDED TO IT, CURRENTLY
   IT REALIGNS UNDER THE IDENTITY METRIC.                                  */

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



typedef struct {
  int apnt, bpnt;  /* A branchpoint occurs at matrix position (apnt,bpnt)
                      in a comparision between two sequences A and B.  That is,
                      theres is a b.p. between the apnt'th and (apnt+1)'st
                      symbol of A and the bpnt'th and (bpnt+1)'st symbo of B.
                      To within first order the length of the matching part of
                      the overlap is MIN(apnt,bpnt) and the length of the
                      non-matching part is MIN(|A|-apnt,|B|-bpnt).         */
  float ascent;    /* The ratio the score of the matching part of the overlap
                      to its length.  If the two sequences match at (1-e)%
                      then this ratio should be about BP_RATIO - e where
                      BP_RATIO is a defined constant in AS_ALN_dpaligners.c
                      currently set to .25.  Thus this number gives one a
                      sense of the fidelity of the match between the two
                      sequences before the branch point.                   */
  float descent;   /* The ratio of the score of the longest extension of the
                      matching part of the overlap into the non-matching part
                      divided by its length.  For random sequences this ratio
                      should be about .5-BP_RATIO.  This number gives one an
                      idea of how divergent the sequences are after the
                      branchpoint, if they don't diverge enough it may be
                      that this is just a more highly variable region within
                      a repeat.                                             */
} BranchPointResult;

BranchPointResult *BPnt_Compare_AS(InternalFragMesg *a, InternalFragMesg *b,
                                   int beg, int end, int opposite,
                                   double erate, double thresh,
                                   int prefix, int suffix);

/* Given fragments a and b, determine if there is a branchpoint between
   them subject to the parameters beg, end, opposite, erate, and thresh as
   for DP_Compare_AS above.  Typically one should use twice the prevailing
   sequencing error rate for erate.  Moreover the branchpoint must occur
   at least minprefix symbols into the overlap and involve at least
   minsuffix symbols beyond the branch.  If there is insufficient memory
   then a NULL pointer is returned, otherwise a branchpoint record is
   returned.  The position of the branchpoint is returned as (-1,-1) if
   a branchpoint satisfying the given constraints was not found.  Otherwise
   a complete record as described with the typedef of the b.p. record
   is returned.  *This record is reused every time the routine is called.* */


BranchPointResult *  BPnt_Seq_Comp_AS
    (char * aseq, int alen,  char * bseq, int blen,
     int beg, int end, double erate, double thresh,
     int minprefix, int minsuffix);

/* Identical to previous except that sequences  aseq  and  bseq  
*  have already been extracted, complemented if necessary,
*  and shifted so that the first characters
*  of each are  aseq [1]  and  bseq [1] , respectively.
*  Their lengths also have been determined to be  alen  and  blen . */

#define MICROMX  6  /* Largest length of a micro-sat to be detected. */

int MicroFinder_AS(char *seq, int ispref, void (*handler)(char *));

/* Find micro-sats of length up to MICROMX(6) that matches a prefix
   (iff "ispref" != 0) or suffix of length MICROTESTLEN(40) of "seq" with
   not more than MICROERRORS(6) errors and call "handler" with each one
   (if any).  Return the first level filter cutoff value (good for
   testing only)
*/




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

	 

 
*/ /*end comments for AS_ALN_affine_overlap */


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


