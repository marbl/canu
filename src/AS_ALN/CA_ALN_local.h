
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

#ifndef CA_ALN_local_h
#define CA_ALN_local_h

static const char *rcsid_CA_ALN_local_h = "$Id: CA_ALN_local.h,v 1.3 2008-10-08 22:02:54 brianwalenz Exp $";

/* Local alignment record:
     Coordinates are in terms of the d.p. matrix that go from (0,0)
     to (|A|,|B|), where A is the A-sequence argument and B is the
     B-sequence argument.  A coordinate is a position between chars
     of the sequence.  For example, an alignment from (3,5) to (6,9)
     aligns characters 4-5 of A with characters 7-9 of B.

     If the B start coordinate is greater than the B end coordinate
     then the alignment is of A versus the complement of B, i.e. the
     alignment runs along an anti-diagonal, not a diagonal.

     Diagonal k is the set of coordinates (a,b) s.t. a-b = k.
     Anti-diagonal k is the set of coordinates (a,b) s.t. a+b = k.
*/

typedef struct {
  int abpos, bbpos;  /* Start coordinate of local alignment */
  int aepos, bepos;  /* End coordinate of local alignment */
  int ldiag, hdiag;  /* Alignment is between (anti)diagonals ldiag & hdiag */
  int   score;       /* Score of alignment where match = 1, difference = -3 */
  float error;       /* Lower bound on error rate of match */
} Local_Segment;

#define LOCAL_FORW 0  /* Compare A to B only */
#define LOCAL_REVR 1  /* Compare A to complement(B) only */
#define LOCAL_BOTH 2  /* Compare A to both B and its complement */

/* Find_Local_Segments compares sequence A of length Alen against sequence B
   of length Blen, in the same, opposite, or both orientations depending
   on the setting of Action, and returns a pointer to an array of local
   alignment records, where the number of such records is in the integer
   pointed at by NumSegs.  Find_Local_Segments only reports alignments that are
   longer than MinLen bps. and contain less than MaxDiff errors as a fraction
   of the alignment's length.  Find_Local_Segments reuses the storage for the
   array of local alignment segments with each successive call, a user should
   copy the array if they which for it to persist beyond a given invocation.

   Find_Local_Segments finds alignments that contain at least 36bp that match
   at 95% or better.  The encompassing local alignment has to match at
   about 75% or better as matching chars are scored SAMECOST (1) and
   differences are scored -DIFFCOST (3).  Extension of a local alignment
   in a given direction ends at a cumulative maximum from which all
   extensions drop by BLOCKCOST (15 = DIFFCOST(3)*MAXIGAP(5)) in score,
   the equivalent of MAXIGAP(5) consecutive differences.

   Find_Local_Segments is most efficiently applied to large sequences.  Any
   application that is doing an all-against-all of smaller fragment
   sequences, would best utilize Find_Local_Segments by concatenating the
   fragments, applying Find_Local_Segments to the large concatenations, and
   then mapping the local alignments back to the fragments and
   coordinates to which they pertain.  By separating fragments in the
   concatentaion by MAXIGAP+1(6) N's, one guarantees that local alignments
   do not span the boundaries between fragment sequences.

   Find_Local_Segments builds an index of the A sequence as part of its
   acceleration method.  If successive calls to Find_Local_Segments involve
   the same A-sequence, this table is built only once, improving
   efficiency.
*/

Local_Segment *Find_Local_Segments
                  (char *A, int Alen, char *B, int Blen,
                   int Action, int MinLen, float MaxDiff, int *NumSegs);

/* Local_Overlap Record:
     A local overlap is a chain of Local_Segments computed by
     Find_Local_Segments that when strung together form an overlap between
     the two sequences involved.  The Local_Overlap record contains the
     the number of segments, a pointer to an ordered array of the segments
     in the chain, and the following parameters:
       score: The score of a chain is the total number of indels required
              to build an overlap out of the chained elements.
       begpos,endpos: As for DP_Compare, the diagonals on which the overlap
              begins and ends.
       diffs: The number of substitutions and indels required to build an
              overlap out of the chained elements.
       length: (|A-seg|+|B-seg|)/2.
     The field chain points to an array of num_pieces+1 Local_Chain records.
     Records 0..num_pieces encode the nature of the gap between segments so
     that record 0 gives the gap to the start border (if any), record
     num_pieces gives the gap to the finish border (if any), and record i
     gives the gap between the segment of record i-1 and record i.  Record
     num_pieces does not contain a segment description.  The parameters
     agap and bgap give the delta in the A- and B-coordinates if of the
     end of the previous segment and the start of the next one.  The
     coordinates can be negative and both can be zero only for the boundary
     gaps (first and last).  Each segment is identified by its position in
     the array of segments passed to the routine, and if the segment was
     complemented in order to form part of the chain, then and only then is
     the field reversed set to a non-zero value.  The type field gives an
     indication of the type of the gap as follows:

       LOCAL_BOUNDARY -- if both agap and bgap are zero at a boundary gap
                         then this indication is given.
       LOCAL_MINOR -- if the a- and b-gaps are less than a user-supplied
                      limit "MinorThresh", then the gap is considered a
                      minor break between two segments of similarity.
       LOCAL_INDEL -- if the gap in one sequence is minor, but major, positive,
                      and at least 4 times as large in the other, then the
                      gap is considered an indel.
       LOCAL_REPEAT -- if the gap in one sequence is minor or negative, but
                       major and negative then the gap is considered a repeat
                       gap in the sense that a tandem repeat must occur in
                       one or both of the sequences around the junction between
                       the two adjacent segments.
       LOCAL_REPnDEL -- if the gap in one sequence is major and negative, and
                        the other is major and positive, then there is a
                        repeated element on both sides of the sequence with
                        the inserted sequence.
       LOCAL_DISAGREE -- anything else, i.e. both gap deltas are positive,
                         at least one is major, and if the other is minor
                         then the ration is less than 1 to 4.
*/

#define LOCAL_BOUNDARY  0x0  /* No gap, at boundary */
#define LOCAL_MINOR     0x1  /* Small break in alignment */
#define LOCAL_DISAGREE  0x2  /* The two sequences significantly disagree */
#define LOCAL_INDEL     0x3  /* One sequence has missing/added sequence */
#define LOCAL_REPEAT    0x4  /* A tandem repeat occurs at the junction */
#define LOCAL_REPnDEL   0x5  /* Both a tandem repeat and an indel */

typedef struct {
  int agap, bgap;      /* A- and B-seq deltas from last segment to this one */
  short type;          /* Type of gap as given by the defined cons. above */
  short reversed;      /* Is segment reversed for inclusion in chain */
  Local_Segment piece; /* Segment in the chain */
} Local_Chain;

typedef struct {
  int begpos;     /* Entry diagonal of boundary point (a,b) on which
                       overlap starts, where diagonal = a - b.          */
  int endpos;     /* Exit diagonal of boundary point (a,b) on which
                       overlap ends, where diagoanl = (|B|-b) - (|A|-a) */
  int length;     /* Length of overlap (|A|+|B|)/2                      */
  int diffs;      /* Estimated number of differences in overlap         */
  int comp;       /* B sequence was complemented for this comparison    */
  int indif;      /* Estimated number of diffs in segments of overlap   */
  int score;      /* Sum of all gap lengths */
  int num_pieces; /* # of segments in overlap chain */
  Local_Chain *chain; /* chain[0..num_pieces] describe each gap between
                         local segments in the overlap chain            */
} Local_Overlap;

/* Find_Local_Overlap takes an array of local alignments as returned by
   Compare_Local and finds the best scoring local overlap between the
   underlying sequences.  One must pass in the length of the two sequences
   from which Compare_Local produced the local alignments as well as the
   number of local alignments in the array.  If the parameter comp is nonzero
   then the comparison will effectively be between A and the complement of B.
   Normally, the parameter nextbest is zero -- after such a call, a second
   alternate overlap, third alternate, and so on can be generated by
   subsequent calls with nextbest set to a nonzero value.  The alternates
   are the best scoring overlaps that starts with a segment not in any
   previous overlap.  The best overlap is returned as a pointer to local
   overlap structure described above.  Unlike many of my routines, the
   reclamation of the storage for this data structure is the responsibility
   of the caller and requires simply calling free on it, as the entire
   structure, including the chain array, is in a single memory block.  The
   parameter MinorThresh determines whether a gap delta is consider minor or
   major (see the description above on gap types).  An overlap is returned
   only if the ratio of the difference to the length of the overlap is less
   than GapThresh, otherwise NULL is returned.
*/

Local_Overlap *Find_Local_Overlap(int Alen, int Blen, int comp, int nextbest,
                                  Local_Segment *Segs, int NumSegs,
                                  int MinorThresh, float GapThresh);

/* Print out a representation of the local overlap, desc, on the file ofile.
   Indent the input by "indent" spaces.                                     */

void Print_Local_Overlap(FILE *ofile, Local_Overlap *desc, int indent);
void Print_Local_Overlap_Picture(FILE *file, Local_Overlap *align, int indent);
void Print_Local_Overlap_withAlign(FILE *file, Local_Overlap *desc,char *a,char *b);

#endif
