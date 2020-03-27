
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "findErrors.H"

#include "sequence.H"

#define  DISPLAY_WIDTH   250

//  Show (to  stdout ) the alignment encoded in  delta [0 .. (deltaLen - 1)]
//  between strings  a [0 .. (a_len - 1)]  and  b [0 .. (b_len - 1)] .

static
void
Display_Alignment(char    *a,   int32 aLen,
                  char    *b,   int32 bLen,
                  int32   *delta,
                  int32    deltaLen) {

  int32  i = 0;
  int32  j = 0;

  char  *top    = new char [32 * 1024];
  int32  topLen = 0;

  char  *bot    = new char [32 * 1024];
  int32  botLen = 0;

  for (int32 k = 0;  k < deltaLen;  k++) {
    for (int32 m = 1;  m < abs(delta[k]);  m++) {
      top[topLen++] = a[i++];
      j++;
    }

    if (delta[k] < 0) {
      top[topLen++] = '-';
      j++;
    } else {
      top[topLen++] = a[i++];
    }
  }

  while (i < aLen && j < bLen) {
    top[topLen++] = a[i++];
    j++;
  }
  top[topLen] = '\0';


  i = j = 0;

  for (int32 k = 0;  k < deltaLen;  k++) {
    for (int32 m = 1;  m < abs(delta[k]);  m++) {
      bot[botLen++] = b[j++];
      i++;
    }

    if (delta[k] > 0) {
      bot[botLen++] = '-';
      i++;
    } else {
      bot[botLen++] = b[j++];
    }
  }

  while (j < bLen && i < aLen) {
    bot[botLen++] = b[j++];
    i++;
  }
  bot[botLen] = '\0';


  for (i = 0;  i < topLen || i < botLen;  i += DISPLAY_WIDTH) {
    putc('\n', stderr);

    fprintf(stderr, "A: ");
    for (j = 0;  j < DISPLAY_WIDTH && i + j < topLen;  j++)
      putc(top[i + j], stderr);
    putc('\n', stderr);

    fprintf(stderr, "B: ");
    for (j = 0;  j < DISPLAY_WIDTH && i + j < botLen;  j++)
      putc(bot[i + j], stderr);
    putc('\n', stderr);

    fprintf(stderr, "   ");
    for (j = 0;  j < DISPLAY_WIDTH && i + j < botLen && i + j < topLen; j++)
      if (top[i + j] != ' ' && bot[i + j] != ' ' && top[i + j] != bot[i + j])
        putc('^', stderr);
      else
        putc(' ', stderr);
    putc('\n', stderr);
  }

  delete [] top;
  delete [] bot;
}


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
                  char   *b_part, //int32 b_len,
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
  fprintf(stderr, "Process_Olap:  %8d %8d %5ld %5ld  %c\n",
          olap->a_iid, olap->b_iid,
          olap->a_hang, olap->b_hang,
          olap->innie == true ? 'I' : 'N');

  //if (olap->a_iid != 39861 && olap->b_iid != 39861 && olap->a_iid != 2283 && olap->b_iid != 2283)
  //  return;
#endif

  int32  ri = olap->a_iid - wa->G->bgnID;

  if ((shredded == true) && (wa->G->reads[ri].shredded == true)) {
    //fprintf(stderr, "%8d %8d shredded\n", olap->a_iid, olap->b_iid);
    return;
  }
  //fprintf(stderr, "%8d %8d not shredded\n", olap->a_iid, olap->b_iid);

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

  //fprintf(stderr, "  errors = %d  delta_len = %d\n", errors, wa->ped.deltaLen);
  //fprintf(stderr, "  a_align = %d/%d  b_align = %d/%d\n", a_end, a_part_len, b_end, b_part_len);
  //Display_Alignment(a_part, a_end, b_part, b_end, wa->ped.delta, wa->ped.deltaLen);//, wa->G->reads[ri].clear_len - a_offset);

  if ((match_to_end == false) && (a_end + a_offset >= wa->G->reads[ri].clear_len - 1)) {
    olap_len = min(a_end, b_end);
    match_to_end = true;
  }

  if (match_to_end && errors <= wa->G->Error_Bound[olap_len]) {
    wa->passedOlaps++;
    //fprintf(stderr, "%8d %8d passed overlap\n", olap->a_iid, olap->b_iid);
    Analyze_Alignment(wa,
                      a_part, a_end, a_offset,
                      b_part, //b_end,
                      ri);
  } else {
    wa->failedOlaps++;
    //fprintf(stderr, "%8d %8d failed overlap\n", olap->a_iid, olap->b_iid);
    //fprintf(stderr, "%8d %8d match to end %c\n", olap->a_iid, olap->b_iid, match_to_end ? 'T' : 'F');
    //fprintf(stderr, "%8d %8d too many errors %c\n", olap->a_iid, olap->b_iid, (errors > wa->G->Error_Bound[olap_len]) ? 'T' : 'F');
  }
}
