
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "findErrors.H"

#include "AS_UTL_reverseComplement.H"


int32
Prefix_Edit_Dist(char   *A, int m,
                 char   *T, int n,
                 int     Error_Limit,
                 int32  &A_End,
                 int32  &T_End,
                 bool   &Match_To_End,
                 pedWorkArea_t * wa);

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
             char               *b_seq,
             bool                shredded,
             Thread_Work_Area_t *wa) {

#if 0
  fprintf(stderr, "Process_Olap:  %8d %8d %5d %5d  %c\n",
          olap->a_iid, olap->b_iid,
          olap->a_hang, olap->b_hang,
          olap->innie == true ? 'I' : 'N');
#endif

  int32  ri = olap->a_iid - wa->G->bgnID;

  if ((shredded == true) && (wa->G->reads[ri].shredded == true))
    return;

  char  *a_part   = wa->G->reads[ri].sequence;
  int32  a_offset = 0;

  char  *b_part   = (olap->normal == true) ? b_seq : wa->rev_seq;
  int32  b_offset = 0;

  //  If innie, reverse-complement the B sequence.

  if ((olap->innie == true) && (wa->rev_id != olap->b_iid)) {
    strcpy(b_part, b_seq);
    reverseComplementSequence(b_part, 0);
    wa->rev_id = olap->b_iid;
  }

  //  Adjust for hangs.

  if (olap->a_hang > 0) {
    a_offset  = olap->a_hang;
    a_part   += a_offset;
  }

  if (olap->a_hang < 0) {
    b_offset  = -olap->a_hang;
    b_part   +=  b_offset;
  }

  //  Count degree - just how many times we cover the end of the read?

  if ((olap->a_hang <= 0) && (wa->G->reads[ri].left_degree < MAX_DEGREE))
    wa->G->reads[ri].left_degree++;

  if ((olap->b_hang >= 0) && (wa->G->reads[ri].right_degree < MAX_DEGREE))
    wa->G->reads[ri].right_degree++;

  // Get the alignment

  uint32   a_part_len = strlen(a_part);
  uint32   b_part_len = strlen(b_part);

  int32    a_end = 0;
  int32    b_end = 0;

  uint32   olap_len = min(a_part_len, b_part_len);

  bool     match_to_end = false;

  //fprintf(stderr, "A: offset %d length %d\n", a_offset, a_part_len);
  //fprintf(stderr, "B: offset %d length %d\n", b_offset, b_part_len);

  int32    errors = Prefix_Edit_Dist(a_part, a_part_len,
                                     b_part, b_part_len,
                                     wa->G->Error_Bound[olap_len],
                                     a_end,
                                     b_end,
                                     match_to_end,
                                     &wa->ped);

  if ((a_end < 0) || (a_end > a_part_len) || (b_end < 0) || (b_end > b_part_len)) {
    fprintf (stderr, "ERROR:  Bad edit distance.\n");
    fprintf (stderr, "  errors = %d  a_end = %d  b_end = %d\n", errors, a_end, b_end);
    fprintf (stderr, "  a_part_len = %d  b_part_len = %d\n", a_part_len, b_part_len);
    fprintf (stderr, "  a_iid = %d  b_iid = %d  match_to_end = %c\n", olap->a_iid, olap->b_iid, match_to_end ? 'T' : 'F');
  }
  assert(a_end >= 0);
  assert(a_end <= a_part_len);
  assert(b_end >= 0);
  assert(b_end <= b_part_len);

  //printf("  errors = %d  delta_len = %d\n", errors, wa->ped.deltaLen);
  //printf("  a_align = %d/%d  b_align = %d/%d\n", a_end, a_part_len, b_end, b_part_len);
  //Display_Alignment(a_part, a_end, b_part, b_end, wa->delta, wa->deltaLen, wa->G->reads[ri].clear_len - a_offset);

  if ((match_to_end == false) && (a_end + a_offset >= wa->G->reads[ri].clear_len - 1)) {
    olap_len = min(a_end, b_end);
    match_to_end = true;
  }

  if ((errors <= wa->G->Error_Bound[olap_len]) && (match_to_end == true))
    Analyze_Alignment(wa,
                      a_part, a_end, a_offset,
                      b_part, b_end,
                      ri);
  else
    wa->failedOlaps++;
}
