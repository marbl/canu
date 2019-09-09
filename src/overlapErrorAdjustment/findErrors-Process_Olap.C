
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2016-JUN-27
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "findErrors.H"

#include "sequence.H"
#include "edlib.H"


//  EDLIB
void
Analyze_Alignment(Thread_Work_Area_t *wa,
                  EdlibAlignResult   &result,
                  char *aseq, int32 abgn, int32 aend,
                  char *bseq, int32 bbgn, int32 bend,
                  int32   sub);

//  ORIGINAL
void
Analyze_Alignment(Thread_Work_Area_t *wa,
                  char   *a_part, int32 a_len, int32 a_offset,
                  char   *b_part, int32 b_len,
                  int32   sub);



//  Find the alignment referred to in  olap , where the  a_iid
//  fragment is in  Frag  and the  b_iid  sequence is in  b_seq .
//  Use the alignment to increment the appropriate vote fields
//  for the a fragment.   shredded  is true iff the b fragment
//  is from shredded data, in which case the overlap will be
//  ignored if the a fragment is also shredded.
//  rev_seq  is a buffer to hold the reverse complement of  b_seq
//  if needed.  (* rev_id) is used to keep track of whether
//  rev_seq  has been created yet.  (* wa) is the work-area
//  containing space for the process to use in case of multi-threading.


void
Process_Olap(Olap_Info_t        *olap,
             char               *b_seq_in,
             bool                shredded,
             Thread_Work_Area_t *wa) {

#if 0
  fprintf(stderr, "Process_Olap:  %8d %8d %5d %5d  %c\n",
          olap->a_iid, olap->b_iid,
          olap->a_hang, olap->b_hang,
          olap->innie == true ? 'I' : 'N');
#endif

  int32  ri = olap->a_iid - wa->G->bgnID;
  int32  rj = olap->b_iid - wa->G->bgnID;

  if ((shredded == true) && (wa->G->reads[ri].shredded == true))
    return;

  char  *a_seq   = wa->G->reads[ri].sequence;
  char  *b_seq   = (olap->normal == true) ? b_seq_in : wa->rev_seq;

  //  If innie, reverse-complement the B sequence.

  if ((olap->innie == true) && (wa->rev_id != olap->b_iid)) {
    strcpy(b_seq, b_seq_in);
    reverseComplementSequence(b_seq, 0);
    wa->rev_id = olap->b_iid;
  }

  //  Count degree - just how many times we cover the end of the read?

  if ((olap->a_hang <= 0) && (wa->G->reads[ri].left_degree < MAX_DEGREE))
    wa->G->reads[ri].left_degree++;

  if ((olap->b_hang >= 0) && (wa->G->reads[ri].right_degree < MAX_DEGREE))
    wa->G->reads[ri].right_degree++;

  //  Extract the read

  int32    abgn = 0;
  int32    aend = wa->G->reads[ri].clear_len;

  int32    bbgn = 0;
  int32    bend = wa->G->reads[rj].clear_len;


  if (olap->a_hang > 0)   abgn +=  olap->a_hang;
  if (olap->a_hang < 0)   bbgn += -olap->a_hang;

  if (olap->b_hang > 0)   bend -=  olap->b_hang;
  if (olap->b_hang < 0)   aend -= -olap->b_hang;

  double maxAlignErate = 0.06;   //  6% is a good default for corrected reads, probably too high.

  EdlibAlignResult result = edlibAlign(a_seq + abgn, aend - abgn,
                                       b_seq + bbgn, bend - bbgn,
                                       edlibNewAlignConfig((int32)ceil(1.1 * maxAlignErate * ((aend - abgn) + (bend - bbgn)) / 2.0),
                                                           EDLIB_MODE_NW,
                                                           EDLIB_TASK_PATH));

#define WITH_DELTA

#ifdef WITH_DELTA
  int32   *delta    = NULL;
  uint32   deltaLen = 0;
  uint32   deltaMax = 0;

  edlibAlignmentToDelta(result.alignment, result.alignmentLength,
                        delta, deltaLen, deltaMax);
#endif

  //  The original code was testing only if (errors <= wa->G->Error_Bound[olap_len]) to decide
  //  if the alignment was any good.  It'd be nice to get rid of that array.

  int32 olap_len = min(aend - abgn, bend - bbgn);

  if ((result.numLocations > 0) &&
      (result.editDistance < wa->G->Error_Bound[olap_len])) {
    wa->passedOlaps++;

    //fprintf(stderr, "numLocations %d editDistance %d limit %d\n",
    //        result.numLocations, result.editDistance, wa->G->Error_Bound[olap_len]);

    Analyze_Alignment(wa,
                      result,
                      a_seq, abgn, aend,
                      b_seq, bbgn, bend,
                      ri);

#ifdef WITH_DELTA
    Analyze_Alignment(wa,
                      a_seq, aend - abgn, abgn,
                      b_seq, bend - bbgn,
                      ri);
#endif
  }

  else {
    //fprintf(stderr, "numLocations %d editDistance %d limit %d\n",
    //        result.numLocations, result.editDistance, wa->G->Error_Bound[olap_len]);
    wa->failedOlaps++;
  }

#ifdef WITH_DELTA
  delete [] delta;
#endif

  edlibFreeAlignResult(result);
}
