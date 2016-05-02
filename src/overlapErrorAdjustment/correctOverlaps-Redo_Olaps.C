
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
#include "AS_UTL_reverseComplement.H"




void
correctRead(uint32 curID,
            char *fseq, uint32 &fseqLen, Adjust_t *fadj, uint32 &fadjLen,
            char *oseq, uint32  oseqLen,
            Correction_Output_t  *C,
            uint64               &Cpos,
            uint64                Clen,
            uint64               *changes=NULL);


int32
Prefix_Edit_Dist(char    *A,  int32 m,
                 char    *T,  int32 n,
                 int32    Error_Limit,
                 int32   &A_End,
                 int32   &T_End,
                 bool    &Match_To_End,
                 pedWorkArea_t *ped);










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

  //fprintf(stderr, "hang adjust delta %d\n", delta);
  return(hang + delta);
}











//  Read old fragments in  gkpStore  and choose the ones that
//  have overlaps with fragments in  Frag. Recompute the
//  overlaps, using fragment corrections and output the revised error.
void
Redo_Olaps(coParameters *G, gkStore *gkpStore) {

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

  fprintf(stderr, "--Allocate "F_U64" MB for fseq and rseq.\n", (2 * sizeof(char) * 2 * (AS_MAX_READLEN + 1)) >> 20);
  char          *fseq    = new char     [AS_MAX_READLEN + 1 + AS_MAX_READLEN + 1];
  uint32         fseqLen = 0;

  char          *rseq    = new char     [AS_MAX_READLEN + 1 + AS_MAX_READLEN + 1];
  uint32         rseqLen = 0;

  fprintf(stderr, "--Allocate "F_U64" MB for fadj and radj.\n", (2 * sizeof(Adjust_t) * (AS_MAX_READLEN + 1)) >> 20);
  Adjust_t      *fadj    = new Adjust_t [AS_MAX_READLEN + 1];
  Adjust_t      *radj    = new Adjust_t [AS_MAX_READLEN + 1];
  uint32         fadjLen  = 0;  //  radj is the same length

  fprintf(stderr, "--Allocate "F_U64" MB for pedWorkArea_t.\n", sizeof(pedWorkArea_t) >> 20);
  gkReadData    *readData = new gkReadData;
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



  ped->initialize(G, G->errorRate);

  //  Process overlaps.  Loop over the B reads, and recompute each overlap.

  for (uint32 curID=loBid; curID<=hiBid; curID++) {
    if (((curID - loBid) % 1024) == 0)
      fprintf(stderr, "Recomputing overlaps - %9u - %9u - %9u\r", loBid, curID, hiBid);

    if (curID < G->olaps[thisOvl].b_iid)
      continue;

    gkRead *read = gkpStore->gkStore_getRead(curID);

    gkpStore->gkStore_loadReadData(read, readData);

    //  Apply corrections to the B read (also converts to lower case, reverses it, etc)

    //fprintf(stderr, "Correcting B read %u at Cpos=%u\n", curID, Cpos);

    fseqLen = 0;
    rseqLen = 0;

    fadjLen = 0;

    correctRead(curID,
                fseq, fseqLen, fadj, fadjLen,
                readData->gkReadData_getSequence(),
                read->gkRead_sequenceLength(),
                C, Cpos, Clen);

    //  Create copies of the sequence for forward and reverse.  There isn't a need for the forward copy (except that
    //  we mutate it with corrections), and the reverse copy could be deferred until it is needed.

    memcpy(rseq, fseq, sizeof(char) * (fseqLen + 1));

    reverseComplementSequence(rseq, fseqLen);

    Make_Rev_Adjust(radj, fadj, fadjLen, fseqLen);

    //  Recompute alignments for all overlaps involving the B read.

    for (; ((thisOvl <= lastOvl) &&
            (G->olaps[thisOvl].b_iid == curID)); thisOvl++) {
      Olap_Info_t  *olap = G->olaps + thisOvl;

      //fprintf(stderr, "processing overlap %u - %u\n", olap->a_iid, olap->b_iid);

      //  Find the A segment.  It's always forward.  It's already been corrected.

      char *a_part = G->reads[olap->a_iid - G->bgnID].bases;

      if (olap->a_hang > 0) {
        int32 ha = Hang_Adjust(olap->a_hang,
                               G->reads[olap->a_iid - G->bgnID].adjusts,
                               G->reads[olap->a_iid - G->bgnID].adjustsLen);
        a_part += ha;
        //fprintf(stderr, "offset a_part by ha=%d\n", ha);
      }

      //  Find the B segment.

      char *b_part = (olap->normal == true) ? fseq : rseq;

      //if (olap->normal == true)
      //  fprintf(stderr, "b_part = fseq %40.40s\n", fseq);
      //else
      //  fprintf(stderr, "b_part = rseq %40.40s\n", rseq);

      if (olap->normal == true)
        olapsFwd++;
      else
        olapsRev++;

      bool rha=false;
      if (olap->a_hang < 0) {
        int32 ha = (olap->normal == true) ? Hang_Adjust(-olap->a_hang, fadj, fadjLen) :
                                            Hang_Adjust(-olap->a_hang, radj, fadjLen);
        b_part += ha;
        //fprintf(stderr, "offset b_part by ha=%d normal=%d\n", ha, olap->normal);
        rha=true;
      }

      //  Compute the alignment.

      int32   a_part_len  = strlen(a_part);
      int32   b_part_len  = strlen(b_part);
      int32   olap_len    = min(a_part_len, b_part_len);

      int32   a_end        = 0;
      int32   b_end        = 0;
      bool    match_to_end = false;

      //fprintf(stderr, ">A\n%s\n", a_part);
      //fprintf(stderr, ">B\n%s\n", b_part);

      int32 errors = Prefix_Edit_Dist(a_part, a_part_len,
                                      b_part, b_part_len,
                                      G->Error_Bound[olap_len],
                                      a_end,
                                      b_end,
                                      match_to_end,
                                      ped);

      //  ped->delta isn't used.

      //  ??  These both occur, but the first is much much more common.

      if ((ped->deltaLen > 0) && (ped->delta[0] == 1) && (0 < G->olaps[thisOvl].a_hang)) {
        int32  stop = min(ped->deltaLen, (int32)G->olaps[thisOvl].a_hang);  //  a_hang is int32:31!
        int32  i = 0;

        for  (i=0; (i < stop) && (ped->delta[i] == 1); i++)
          ;

        //fprintf(stderr, "RESET 1 i=%d delta=%d\n", i, ped->delta[i]);
        assert((i == stop) || (ped->delta[i] != -1));

        ped->deltaLen -= i;

        memmove(ped->delta, ped->delta + i, ped->deltaLen * sizeof (int));

        a_part     += i;
        a_end      -= i;
        a_part_len -= i;
        errors     -= i;

      } else if ((ped->deltaLen > 0) && (ped->delta[0] == -1) && (G->olaps[thisOvl].a_hang < 0)) {
        int32  stop = min(ped->deltaLen, - G->olaps[thisOvl].a_hang);
        int32  i = 0;

        for  (i=0; (i < stop) && (ped->delta[i] == -1); i++)
          ;

        //fprintf(stderr, "RESET 2 i=%d delta=%d\n", i, ped->delta[i]);
        assert((i == stop) || (ped->delta[i] != 1));

        ped->deltaLen -= i;

        memmove(ped->delta, ped->delta + i, ped->deltaLen * sizeof (int));

        b_part     += i;
        b_end      -= i;
        b_part_len -= i;
        errors     -= i;
      }


      Total_Alignments_Ct++;


      int32  olapLen = min(a_end, b_end);

      if ((match_to_end == false) && (olapLen <= 0))
        Failed_Alignments_Both_Ct++;

      if (match_to_end == false)
        Failed_Alignments_End_Ct++;

      if (olapLen <= 0)
        Failed_Alignments_Length_Ct++;

      if ((match_to_end == false) || (olapLen <= 0)) {
        Failed_Alignments_Ct++;

#if 0
        //  I can't find any patterns in these errors.  I thought that it was caused by the corrections, but I
        //  found a case where no corrections were made and the alignment still failed.  Perhaps it is differences
        //  in the alignment code (the forward vs reverse prefix distance in overlapper vs only the forward here)?

        fprintf(stderr, "Redo_Olaps()--\n");
        fprintf(stderr, "Redo_Olaps()--\n");
        fprintf(stderr, "Redo_Olaps()--  Bad alignment  errors %d  a_end %d  b_end %d  match_to_end %d  olapLen %d\n",
                errors, a_end, b_end, match_to_end, olapLen);
        fprintf(stderr, "Redo_Olaps()--  Overlap        a_hang %d b_hang %d innie %d\n",
                olap->a_hang, olap->b_hang, olap->innie);
        fprintf(stderr, "Redo_Olaps()--  Reads          a_id %u a_length %d b_id %u b_length %d\n",
                G->olaps[thisOvl].a_iid,
                G->reads[ G->olaps[thisOvl].a_iid ].basesLen,
                G->olaps[thisOvl].b_iid,
                G->reads[ G->olaps[thisOvl].b_iid ].basesLen);
        fprintf(stderr, "Redo_Olaps()--  A %s\n", a_part);
        fprintf(stderr, "Redo_Olaps()--  B %s\n", b_part);

        Display_Alignment(a_part, a_part_len, b_part, b_part_len, ped->delta, ped->deltaLen);

        fprintf(stderr, "\n");
#endif

        if (rha)
          rhaFail++;

        continue;
      }

      if (rha)
        rhaPass++;

      G->olaps[thisOvl].evalue = AS_OVS_encodeEvalue((double)errors / olapLen);

      //fprintf(stderr, "REDO - errors = %u / olapLep = %u -- %f\n", errors, olapLen, AS_OVS_decodeEvalue(G->olaps[thisOvl].evalue));
    }
  }

  fprintf(stderr, "\n");

  delete    ped;
  delete    readData;
  delete [] radj;
  delete [] fadj;
  delete [] rseq;
  delete [] fseq;
  delete    Cfile;

  fprintf(stderr, "--  Release bases, adjusts and reads.\n");

  delete [] G->bases;     G->bases   = NULL;
  delete [] G->adjusts;   G->adjusts = NULL;
  delete [] G->reads;     G->reads   = NULL;

  fprintf(stderr, "Olaps Fwd "F_U64"\n", olapsFwd);
  fprintf(stderr, "Olaps Rev "F_U64"\n", olapsRev);

  fprintf(stderr, "Total:  "F_U64"\n", Total_Alignments_Ct);
  fprintf(stderr, "Failed: "F_U64" (both)\n", Failed_Alignments_Both_Ct);
  fprintf(stderr, "Failed: "F_U64" (either)\n", Failed_Alignments_Ct);
  fprintf(stderr, "Failed: "F_U64" (match to end)\n", Failed_Alignments_End_Ct);
  fprintf(stderr, "Failed: "F_U64" (negative length)\n", Failed_Alignments_Length_Ct);

  fprintf(stderr, "rhaFail %u rhaPass %u\n", rhaFail, rhaPass);
}
