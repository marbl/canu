
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
 *    Brian P. Walenz from 2015-MAY-14 to 2015-JUN-03
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-FEB-25
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "correctOverlaps.H"
#include "correctionOutput.H"

#include "sequence.H"
#include "edlib.H"



void
correctRead(uint32 curID,
            char *fseq, uint32 &fseqLen, Adjust_t *fadj, uint32 &fadjLen,
            char *oseq, uint32  oseqLen,
            Correction_Output_t  *C,
            uint64               &Cpos,
            uint64                Clen,
            uint64               *changes=NULL);











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













//  Set hanging offset values for reversed fragment in
//   rev_adj[0 .. (adj_ct - 1)]  based on corresponding forward
//  values in  fadj[0 .. (adj_ct - 1)].  frag_len  is the length
//  of the fragment.

static
void
Make_Rev_Adjust(Adjust_t    *radj,
                Adjust_t    *fadj,
                int32        adj_ct,
                int32        frag_len) {

  if (adj_ct == 0)
    return;

  int32  i = 0;
  int32  j = 0;
  int32  prev = 0;

  for (i=adj_ct-1; i>0; i--) {
    if (fadj[i].adjust == fadj[i-1].adjust + 1) {
      radj[j].adjpos = 2 + frag_len - fadj[i].adjpos;
      radj[j].adjust = prev + 1;

      prev = radj[j].adjust;
    }

    else if (fadj[i].adjust == fadj[i-1].adjust - 1) {
      radj[j].adjpos = 3 + frag_len - fadj[i].adjpos;
      radj[j].adjust = prev - 1;

      prev = radj[j].adjust;
    }

    else {
      fprintf(stderr, "ERROR:  Bad adjustment value.  i = %d  adj_ct = %d  adjust[i] = %d  adjust[i-1] = %d\n",
              i, adj_ct, fadj[i].adjust, fadj[i-1].adjust);
      assert(0);
    }

    j++;
  }

  assert(i == 0);

  if (fadj[i].adjust == 1) {
    radj[j].adjpos = 2 + frag_len - fadj[i].adjpos;
    radj[j].adjust = prev + 1;
  }

  else if (fadj[i].adjust == -1) {
    radj[j].adjpos = 3 + frag_len - fadj[i].adjpos;
    radj[j].adjust = prev - 1;
  }

  else {
    fprintf(stderr, "ERROR:  Bad adjustment value.  i = %d  adj_ct = %d  adjust[i] = %d\n",
             i, adj_ct, fadj[i].adjust);
    assert(0);
  }

  assert(j+1 == adj_ct);
}





//  Return the adjusted value of  hang  based on
//   adjust[0 .. (adjust_ct - 1)] .
static
int32
Hang_Adjust(int32     hang,
            Adjust_t *adjust,
            int32     adjust_ct) {
  int32  delta = 0;

  assert(hang >= 0);

  //  Replacing second test >= with just > didn't change anything.  Both had 14 fails.

  for  (int32 i=0; (i < adjust_ct) && (hang >= adjust[i].adjpos); i++) {
    //if (delta != adjust[i].adjust)
    //  fprintf(stderr, "hang_adjust i=%d adjust_ct=%d adjust=%d pos=%d\n", i, adjust_ct, adjust[i].adjust, adjust[i].adjpos);
    delta = adjust[i].adjust;
  }

  if (hang + delta < 0) {
    int32 i=0;

    fprintf(stderr, "\n");
    fprintf(stderr, "hang_adjust hang=%d\n", hang);

    for  (; (i < adjust_ct) && (hang >= adjust[i].adjpos); i++)
      fprintf(stderr, "hang_adjust i=%d adjust_ct=%d adjust=%d pos=%d --\n", i, adjust_ct, adjust[i].adjust, adjust[i].adjpos);

    for  (int32 j=i+10; (i < adjust_ct) && (i < j); i++)
      fprintf(stderr, "hang_adjust i=%d adjust_ct=%d adjust=%d pos=%d\n", i, adjust_ct, adjust[i].adjust, adjust[i].adjpos);

    return(0);
  }

  //fprintf(stderr, "hang adjust delta %d\n", delta);
  return(hang + delta);
}











//  Read old fragments in  seqStore  and choose the ones that
//  have overlaps with fragments in  Frag. Recompute the
//  overlaps, using fragment corrections and output the revised error.
void
Redo_Olaps(coParameters *G, sqStore *seqStore) {

  //  Figure out the range of B reads we care about.  We probably could just loop over every read in
  //  the store with minimal penalty.

  uint64     thisOvl = 0;
  uint64     lastOvl = G->olapsLen - 1;

  uint32     loBid   = G->olaps[thisOvl].b_iid;
  uint32     hiBid   = G->olaps[lastOvl].b_iid;

  //  Open all the corrections.

  memoryMappedFile     *Cfile = new memoryMappedFile(G->correctionsName);
  Correction_Output_t  *C     = (Correction_Output_t *)Cfile->get();
  uint64                Cpos  = 0;
  uint64                Clen  = Cfile->length() / sizeof(Correction_Output_t);

  //  Allocate some temporary work space for the forward and reverse corrected B reads.

  fprintf(stderr, "--Allocate " F_SIZE_T " MB for fseq and rseq.\n", (2 * sizeof(char) * 2 * (AS_MAX_READLEN + 1)) >> 20);
  char          *fseq    = new char     [AS_MAX_READLEN + 1 + AS_MAX_READLEN + 1];
  uint32         fseqLen = 0;

  char          *rseq    = new char     [AS_MAX_READLEN + 1 + AS_MAX_READLEN + 1];

  fprintf(stderr, "--Allocate " F_SIZE_T " MB for fadj and radj.\n", (2 * sizeof(Adjust_t) * (AS_MAX_READLEN + 1)) >> 20);
  Adjust_t      *fadj    = new Adjust_t [AS_MAX_READLEN + 1];
  Adjust_t      *radj    = new Adjust_t [AS_MAX_READLEN + 1];
  uint32         fadjLen  = 0;  //  radj is the same length

  fprintf(stderr, "--Allocate " F_SIZE_T " MB for pedWorkArea_t.\n", sizeof(pedWorkArea_t) >> 20);
  sqRead        *read     = new sqRead;
  pedWorkArea_t *ped      = new pedWorkArea_t;

  uint64         Total_Alignments_Ct           = 0;

  uint64         Failed_Alignments_Ct          = 0;
  uint64         Failed_Alignments_Both_Ct     = 0;
  uint64         Failed_Alignments_End_Ct      = 0;
  uint64         Failed_Alignments_Length_Ct   = 0;

  uint32         rhaFail = 0;
  uint32         rhaPass = 0;

  uint64         olapsFwd = 0;
  uint64         olapsRev = 0;

  uint64         nBetter = 0;
  uint64         nWorse  = 0;
  uint64         nSame   = 0;


  ped->initialize(G, G->errorRate);

  //  Process overlaps.  Loop over the B reads, and recompute each overlap.

  for (uint32 curID=loBid; curID<=hiBid; curID++) {
    if (((curID - loBid) % 1024) == 0)
      fprintf(stderr, "Recomputing overlaps - %9u - %9u - %9u\r", loBid, curID, hiBid);

    //  Overlaps are sorted on the B iid.  If the current B iid isn't used in
    //  an overlap, skip the loop.

    if (curID < G->olaps[thisOvl].b_iid)
      continue;

    //  Fetch the B read and apply corrections to it (also convert to lower
    //  case, reverse, etc).

    //fprintf(stderr, "Correcting B read %u at Cpos=%u Clen=%u\n", curID, Cpos, Clen);

    seqStore->sqStore_getRead(curID, read);

    fseqLen = 0;
    fadjLen = 0;

    correctRead(curID,
                fseq, fseqLen, fadj, fadjLen,
                read->sqRead_sequence(),
                read->sqRead_length(),
                C, Cpos, Clen);

    //fprintf(stderr, "Finished   B read %u at Cpos=%u Clen=%u\n", curID, Cpos, Clen);

    memcpy(rseq, fseq, sizeof(char) * (fseqLen + 1));

    reverseComplementSequence(rseq, fseqLen);

    Make_Rev_Adjust(radj, fadj, fadjLen, fseqLen);

    //  Recompute alignments for all overlaps involving the B read.

    for (; ((thisOvl <= lastOvl) &&
            (G->olaps[thisOvl].b_iid == curID)); thisOvl++) {
      Olap_Info_t  *olap = G->olaps + thisOvl;

      //fprintf(stderr, "processing overlap %u - %u\n", olap->a_iid, olap->b_iid);

      //  Find the A segment.  It's always forward.  It's already been corrected.

      char    *a_seq = G->reads[olap->a_iid - G->bgnID].bases;
      char    *b_seq = (olap->normal == true) ? fseq : rseq;

      int32    abgn = 0;
      int32    aend = G->reads[ olap->a_iid - G->bgnID ].basesLen;

      int32    bbgn = 0;
      int32    bend = fseqLen;

      //  Adjust the bgn/end points based on corrections made to the read.
      //
      //  EXCEPT we can't adjust the B read because we don't have those corrections
      //  anymore.
      //
      //  rha was set if olap->a_hang < 0 -- if we adjusted the begin of the B read
      //  from zero to something else.

#if 0
      if (olap->a_hang > 0)   abgn +=Hang_Adjust(olap->a_hang,
                                                 G->reads[olap->a_iid - G->bgnID].adjusts,
                                                 G->reads[olap->a_iid - G->bgnID].adjustsLen);
      if (olap->a_hang < 0)   bbgn += ((olap->normal == true) ? Hang_Adjust(-olap->a_hang, fadj, fadjLen) :
                                                                Hang_Adjust(-olap->a_hang, radj, fadjLen));

      if (olap->b_hang > 0)   bend -=  olap->b_hang;
      if (olap->b_hang < 0)   aend -= -olap->b_hang;
#endif

      if (olap->a_hang > 0)   abgn +=  olap->a_hang;
      if (olap->a_hang < 0)   bbgn += -olap->a_hang;

      if (olap->b_hang > 0)   bend -=  olap->b_hang;
      if (olap->b_hang < 0)   aend -= -olap->b_hang;


      if (olap->normal == true)
        olapsFwd++;
      else
        olapsRev++;

      Total_Alignments_Ct++;


      double maxAlignErate = 0.06;   //  6% is a good default for corrected reads, probably too high.

      EdlibAlignResult result = edlibAlign(a_seq + abgn, aend - abgn,
                                           b_seq + bbgn, bend - bbgn,
                                           edlibNewAlignConfig((int32)ceil(1.1 * maxAlignErate * ((aend - abgn) + (bend - bbgn)) / 2.0),
                                                               EDLIB_MODE_NW,
                                                               EDLIB_TASK_PATH));


      int32 olap_len = min(aend - abgn, bend - bbgn);

#if 0
      if ((result.numLocations > 0) &&
          (result.editDistance < wa->G->Error_Bound[olap_len])) {
        wa->passedOlaps++;

      }

      else {
#endif


      if        ((double)result.editDistance / olap_len < G->olaps[thisOvl].evalue)
        nBetter++;
      else if ((double)result.editDistance / olap_len > G->olaps[thisOvl].evalue)
        nWorse++;
      else
        nSame++;


      G->olaps[thisOvl].evalue = AS_OVS_encodeEvalue((double)result.editDistance / olap_len);

      edlibFreeAlignResult(result);

      //fprintf(stderr, "REDO - errors = %u / olapLep = %u -- %f\n", errors, olapLen, AS_OVS_decodeEvalue(G->olaps[thisOvl].evalue));
    }
  }

  fprintf(stderr, "\n");

  delete    ped;
  delete    read;
  delete [] radj;
  delete [] fadj;
  delete [] rseq;
  delete [] fseq;
  delete    Cfile;

  fprintf(stderr, "--  Release bases, adjusts and reads.\n");

  delete [] G->bases;     G->bases   = NULL;
  delete [] G->adjusts;   G->adjusts = NULL;
  delete [] G->reads;     G->reads   = NULL;

  fprintf(stderr, "Olaps Fwd " F_U64 "\n", olapsFwd);
  fprintf(stderr, "Olaps Rev " F_U64 "\n", olapsRev);

  fprintf(stderr, "Total:  " F_U64 "\n", Total_Alignments_Ct);
  fprintf(stderr, "Failed: " F_U64 " (both)\n", Failed_Alignments_Both_Ct);
  fprintf(stderr, "Failed: " F_U64 " (either)\n", Failed_Alignments_Ct);
  fprintf(stderr, "Failed: " F_U64 " (match to end)\n", Failed_Alignments_End_Ct);
  fprintf(stderr, "Failed: " F_U64 " (negative length)\n", Failed_Alignments_Length_Ct);

  fprintf(stderr, "rhaFail %u rhaPass %u\n", rhaFail, rhaPass);

  fprintf(stderr, "Changed %lu overlaps.\n", nBetter + nWorse + nSame);
  fprintf(stderr, "Better: %lu overlaps.\n", nBetter);
  fprintf(stderr, "Worse:  %lu overlaps.\n", nWorse);
  fprintf(stderr, "Same:   %lu overlaps.\n", nSame);
}
