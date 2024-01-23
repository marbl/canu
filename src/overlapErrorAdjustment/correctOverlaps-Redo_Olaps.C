
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

#include "correctOverlaps.H"
#include "correctionOutput.H"
#include "sequence.H"
#include "computeDiff.H"

#include <tuple>



void
correctRead(uint32 curID,
            char *fseq, uint32 &fseqLen, Adjust_t *fadj, uint32 &fadjLen,
            const char *oseq, uint32  oseqLen,
            Correction_Output_t  *C,
            uint64               &Cpos,
            uint64                Clen,
            uint64               *changes=NULL);


int32
Prefix_Edit_Dist(const char    *A,  int32 m,
                 const char    *T,  int32 n,
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
Display_Alignment(const char    *a,   int32 aLen,
                  const char    *b,   int32 bLen,
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

static
void
PrepareRead(/*const*/ sqStore *seqStore, uint32 curID,
            uint32 &fseqLen, char *fseq, char *rseq,
            uint32 &fadjLen, Adjust_t *fadj, Adjust_t *radj,
            Correction_Output_t  *C, uint64 &Cpos, uint64 Clen) {
  sqRead read;
  seqStore->sqStore_getRead(curID, &read);
  //  Apply corrections to the B read (also converts to lower case, reverses it, etc)

  //fprintf(stderr, "Correcting B read %u at Cpos=%u Clen=%u\n", curID, Cpos, Clen);

  fseqLen = 0;
  fadjLen = 0;

  //Correcting "b" read. "a" reads were corrected beforehand.
  correctRead(curID,
              fseq, fseqLen, fadj, fadjLen,
              read.sqRead_sequence(),
              read.sqRead_length(),
              C, Cpos, Clen);

  //fprintf(stderr, "Finished   B read %u at Cpos=%u Clen=%u\n", curID, Cpos, Clen);

  //  Create copies of the sequence for forward and reverse.  There isn't a need for the forward copy (except that
  //  we mutate it with corrections), and the reverse copy could be deferred until it is needed.

  memcpy(rseq, fseq, sizeof(char) * (fseqLen + 1));

  reverseComplementSequence(rseq, fseqLen);

  Make_Rev_Adjust(radj, fadj, fadjLen, fseqLen);
}

//FIXME code duplication
//returns error rate of the alignment or -1. if (!match_to_end || invalid_olap)
static
double
ProcessAlignment(const int32 a_part_len, char *a_part,
                 const int32 b_part_len, char *b_part,
                 int32 error_bound, bool check_trivial_dna,
                 pedWorkArea_t *ped, bool *match_to_end, bool *invalid_olap) {
  int32   a_end        = 0;
  int32   b_end        = 0;

  //fprintf(stderr, ">A\n%s\n", a_part);
  //fprintf(stderr, ">B\n%s\n", b_part);

  int32 all_errors = Prefix_Edit_Dist(a_part, a_part_len,
                                      b_part, b_part_len,
                                      error_bound,
                                      a_end,
                                      b_end,
                                      *match_to_end,
                                      ped);


  //Adjusting the extremities
  //TODO discuss the logic!
  //TODO discuss changes regarding a_hang usage!
  //TODO refactor out the code duplication
  all_errors -= TrimStartingIndels(a_part, a_end, ped->delta, ped->deltaLen, 1);
  all_errors -= TrimStartingIndels(b_part, b_end, ped->delta, ped->deltaLen, -1);

  //fprintf(stderr, "Showing alignment\n");
  //Display_Alignment(a_part, a_end, b_part, b_end, ped->delta, ped->deltaLen);

  *invalid_olap = (std::min(a_end, b_end) <= 0);

  if (!*match_to_end || *invalid_olap) {
    return -1.;
  }

  int32 events;
  int32 alignment_len;

  //fprintf(stderr, "Computing errors\n");
  //fprintf(stderr, "Checking for trivial DNA regions: %d\n", check_trivial_dna);
  std::tie(events, alignment_len) = ComputeErrors(a_part, b_part, ped->deltaLen, ped->delta,
                                                  a_end, b_end, check_trivial_dna);

  //fprintf(stderr, "Old errors %d new events %d\n", all_errors, events);

  //Assert can fail if DEFAULT_FLANK_IGNORE > 0
  //assert(check_trivial_dna || all_errors == events);

  assert(events >= 0 && alignment_len > 0);
  assert(all_errors > int32min);  //  Useless but makes all_errors be 'used',
  return (double) events / alignment_len;
}

//  Read old fragments in  seqStore  and choose the ones that
//  have overlaps with fragments in  Frag. Recompute the
//  overlaps, using fragment corrections and output the revised error.
void
Redo_Olaps(coParameters *G, /*const*/ sqStore *seqStore) {

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
  pedWorkArea_t *ped      = new pedWorkArea_t;

  uint64         Total_Alignments_Ct           = 0;

  uint64         Failed_Alignments_Ct          = 0;
  uint64         Failed_Alignments_Both_Ct     = 0;
  uint64         Failed_Alignments_End_Ct      = 0;
  uint64         Failed_Alignments_Length_Ct   = 0;

  uint64         olapsFwd = 0;
  uint64         olapsRev = 0;

  uint64         nBetter = 0;
  uint64         nWorse  = 0;
  uint64         nSame   = 0;


  ped->initialize(G, G->errorRate);

  //  Process overlaps.  Loop over the B reads, and recompute each overlap.
  //  Loop over the B reads ...
  for (uint32 curID = loBid; curID <= hiBid; curID++) {
    if (((curID - loBid) % 1024) == 0)
      fprintf(stderr, "Recomputing overlaps - %9u - %9u - %9u\n", loBid, curID, hiBid);

    //  ... towards the read for current overlap
    if (curID < G->olaps[thisOvl].b_iid)
      continue;

    assert(curID == G->olaps[thisOvl].b_iid);

    //  Load and correct the B read
    PrepareRead(seqStore, curID,
                fseqLen, fseq, rseq,
                fadjLen, fadj, radj,
                C, Cpos, Clen);

    //  Recompute alignments for ALL overlaps involving the B read
    for (; thisOvl <= lastOvl && G->olaps[thisOvl].b_iid == curID; thisOvl++) {
      const Olap_Info_t &olap = G->olaps[thisOvl];

      if (G->secondID != UINT32_MAX && olap.b_iid != G->secondID)
        continue;

      if (olap.normal) {
      //  fprintf(stderr, "b_part = fseq %40.40s\n", fseq);
        olapsFwd++;
      } else {
      //  fprintf(stderr, "b_part = rseq %40.40s\n", rseq);
        olapsRev++;
      }

      //  Find the A segment.  It's always forward.  It's already been corrected.
      char *a_part = G->reads[olap.a_iid - G->bgnID].bases;
      if (olap.a_hang > 0) {
        int32 ha = Hang_Adjust(olap.a_hang,
                               G->reads[olap.a_iid - G->bgnID].adjusts,
                               G->reads[olap.a_iid - G->bgnID].adjustsLen);
        a_part += ha;
        //fprintf(stderr, "offset a_part by ha=%d\n", ha);
      }

      //  Find the B segment.
      char *b_part = (olap.normal == true) ? fseq : rseq;

      if (olap.a_hang < 0) {
        int32 ha = olap.normal ? Hang_Adjust(-olap.a_hang, fadj, fadjLen) :
                                            Hang_Adjust(-olap.a_hang, radj, fadjLen);
        b_part += ha;
        //fprintf(stderr, "offset b_part by ha=%d normal=%d\n", ha, olap.normal);
      }

      //  Compute and process the alignment
      Total_Alignments_Ct++;
      //TODO discuss difference with error finding code
      //In errors finding one of the sequences is the (almost) entire read and the length of its prefix is passed
      int32   a_part_len  = strlen(a_part);
      int32   b_part_len  = strlen(b_part);

      bool    match_to_end = false;
      bool    invalid_olap = false;
      double err_rate = ProcessAlignment(a_part_len, a_part, //olap.a_hang,
                                         b_part_len, b_part,
                                         G->Error_Bound[std::min(a_part_len, b_part_len)],
                                         /*check trivial DNA*/G->checkTrivialDNA,
                                         ped, &match_to_end, &invalid_olap);

      if (err_rate >= 0.) {
        //if (err_rate > /*report_threshold*/ 0.) {
        //  fprintf(stderr, "Err rate of overlap %u - %u is %f\n", olap.a_iid, olap.b_iid, err_rate);
        //}

        const uint32 err_encoded = AS_OVS_encodeEvalue(err_rate);

        const uint32 base_encoded = G->olaps[thisOvl].evalue;
        if (err_encoded < base_encoded)
          nBetter++;
        else if (err_encoded > base_encoded)
          nWorse++;
        else
          nSame++;

        G->olaps[thisOvl].evalue = err_encoded;
        //fprintf(stderr, "REDO - err rate = %f\n", AS_OVS_decodeEvalue(G->olaps[thisOvl].evalue));
      } else {
        //fprintf(stderr, "Err rate of overlap %u - %u failed\n", olap.a_iid, olap.b_iid);

        Failed_Alignments_Ct++;

        if (!match_to_end && invalid_olap)
          Failed_Alignments_Both_Ct++;

        if (!match_to_end)
          Failed_Alignments_End_Ct++;

        if (invalid_olap)
          Failed_Alignments_Length_Ct++;

      #if 0
        //  I can't find any patterns in these errors.  I thought that it was caused by the corrections, but I
        //  found a case where no corrections were made and the alignment still failed.  Perhaps it is differences
        //  in the alignment code (the forward vs reverse prefix distance in overlapper vs only the forward here)?

        fprintf(stderr, "Redo_Olaps()--\n");
        fprintf(stderr, "Redo_Olaps()--\n");
        fprintf(stderr, "Redo_Olaps()--  Bad alignment  errors %d  a_end %d  b_end %d  match_to_end %d  olapLen %d\n",
                errors, a_end, b_end, match_to_end, olapLen);
        fprintf(stderr, "Redo_Olaps()--  Overlap        a_hang %d b_hang %d innie %d\n",
                olap.a_hang, olap.b_hang, olap.innie);
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
      }
    }
  }

  fprintf(stderr, "\n");

  delete    ped;
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

  fprintf(stderr, "Changed " F_U64 " overlaps.\n", nBetter + nWorse + nSame);
  fprintf(stderr, "Better: " F_U64 " overlaps.\n", nBetter);
  fprintf(stderr, "Worse:  " F_U64 " overlaps.\n", nWorse);
  fprintf(stderr, "Same:   " F_U64 " overlaps.\n", nSame);
}
