

#include "correctOverlaps.H"
#include "correctionOutput.H"
#include "AS_UTL_reverseComplement.H"




void
correctRead(uint32 curID,
            char *fseq, uint32 &fseqLen, Adjust_t *fadj, uint32 &fadjLen,
            char *oseq, uint32  oseqLen,
            Correction_Output_t  *C,
            uint64               &Cpos,
            uint64                Clen);

int
Prefix_Edit_Dist(char A[], int m, char T[], int n, int Error_Limit,
                     int &A_End, int &T_End, bool &Match_To_End,
                     int *Delta, int &Delta_Len) {
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
       radj[j].pos    = 2 + frag_len - fadj[i].pos;
       radj[j].adjust = prev + 1;

       prev = radj[j].adjust;
     }

     else if (fadj[i].adjust == fadj[i-1].adjust - 1) {
       radj[j].pos    = 3 + frag_len - fadj[i].pos;
       radj[j].adjust = prev - 1;

       prev = radj[j].adjust;
     }

     else {
       fprintf (stderr, "ERROR:  Bad adjustment value.  i = %d  adj_ct = %d  adjust[i] = %d  adjust[i-1] = %d\n",
                i, adj_ct, fadj[i].adjust, fadj[i-1].adjust);
       assert(0);
     }

     j++;
   }

   assert(i == 0);

   if (fadj[i].adjust == 1) {
     radj[j].pos    = 2 + frag_len - fadj[i].pos;
     radj[j].adjust = prev + 1;
   }

   else if (fadj[i].adjust == -1) {
     radj[j].pos    = 3 + frag_len - fadj[i].pos;
     radj[j].adjust = prev - 1;
   }

   else {
     fprintf (stderr, "ERROR:  Bad adjustment value.  i = %d  adj_ct = %d  adjust[i] = %d\n",
              i, adj_ct, fadj[i].adjust);
     assert(0);
   }
  }





//  Return the adjusted value of  hang  based on
//   adjust[0 .. (adjust_ct - 1)] .
static
int32
Hang_Adjust(int32     hang,
            Adjust_t *adjust,
            int       adjust_ct) {
  int  delta = 0;

  for  (int32 i=0; (i < adjust_ct) && (hang >= adjust[i].pos); i++)
    delta = adjust[i].adjust;

  return(hang + delta);
}











//  Read old fragments in  gkpStore  and choose the ones that
//  have overlaps with fragments in  Frag. Recompute the
//  overlaps, using fragment corrections and output the revised error.
void
Redo_Olaps(coParameters &G, gkStore *gkpStore) {

  //  Figure out the range of B reads we care about.  We probably could just loop over every read in
  //  the store with minimal penalty.

  uint64     thisOvl = 0;
  uint64     lastOvl = G.olapsLen - 1;

  uint32     loBid   = G.olaps[thisOvl].b_iid;
  uint32     hiBid   = G.olaps[lastOvl].b_iid;

  //  Open all the corrections.

  memoryMappedFile     *Cfile = new memoryMappedFile(G.correctionsName);
  Correction_Output_t  *C     = (Correction_Output_t *)Cfile->get();
  uint64                Cpos  = 0;
  uint64                Clen  = Cfile->length() / sizeof(Correction_Output_t);

  //  Allocate some temporary work space for the forward and reverse corrected B reads.

  char         *fseq    = new char     [AS_MAX_READLEN + AS_MAX_READLEN];
  uint32        fseqLen = 0;
  Adjust_t     *fadj    = new Adjust_t [AS_MAX_READLEN];
  uint32        fadjLen  = 0;

  char         *rseq    = new char     [AS_MAX_READLEN + AS_MAX_READLEN];
  uint32        rseqLen = 0;
  Adjust_t     *radj    = new Adjust_t [AS_MAX_READLEN];
  uint32        radjLen = 0;

  gkReadData    readData;

  uint64  Total_Alignments_Ct  = 0;
  uint64  Failed_Alignments_Ct = 0;

  //  Process overlaps.  Loop over the B reads, and recompute each overlap.

  for (uint32 curID=loBid; curID<hiBid; curID++) {

    if (curID < G.olaps[thisOvl].b_iid)
      continue;

    gkRead *read = gkpStore->gkStore_getRead(curID);

    if (read->gkRead_isDeleted() == true)
      continue;

    gkpStore->gkStore_loadReadData(read, &readData);

    //  Apply corrections to the B read (also converts to lower case, reverses it, etc)

    correctRead(curID,
                fseq, fseqLen, fadj, fadjLen,
                readData.gkReadData_getSequence(),
                read->gkRead_clearRegionLength(),
                C, Cpos, Clen);

    //  Create copies of the sequence for forward and reverse.  There isn't a need for the forward copy (except that
    //  we mutate it with corrections), and the reverse copy could be deferred until it is needed.

    memcpy(rseq, fseq, sizeof(char) * (fseqLen + 1));

    reverseComplementSequence(rseq, fseqLen);

    Make_Rev_Adjust(radj, fadj, fadjLen, fseqLen);

    //  Recompute alignments for all overlaps involving the B read.

    for (; ((thisOvl <= lastOvl) &&
            (G.olaps[thisOvl].b_iid == curID)); thisOvl++) {
      Olap_Info_t  *olap = G.olaps + thisOvl;

      //  Find the A segment.  It's always forward.  It's already been corrected.

      char *a_part = G.reads[olap->a_iid].bases;

      if (olap->a_hang > 0)
        a_part += Hang_Adjust(olap->a_hang, G.reads[curID].adjusts, G.reads[curID].adjustsLen);

      //  Find the B segment.

      char *b_part = (olap->normal == true) ? fseq : rseq;

      if (olap->a_hang < 0)
        b_part += (olap->normal == true) ? Hang_Adjust(-olap->a_hang, fadj, fadjLen) :
                                           Hang_Adjust(-olap->a_hang, radj, radjLen);

      //  Compute the alignment.

      int32   a_part_len  = strlen(a_part);
      int32   b_part_len  = strlen(b_part);
      int32   olap_len    = min(a_part_len, b_part_len);

      int32   a_end        = 0;
      int32   b_end        = 0;
      bool    match_to_end = false;
      int32  *delta        = NULL;
      int32   deltaLen     = 0;

      int32 errors = Prefix_Edit_Dist(a_part, a_part_len,
                                      b_part, b_part_len,
                                      G.Error_Bound[olap_len],
                                      a_end,
                                      b_end,
                                      match_to_end,
                                      delta, deltaLen);

      if ((deltaLen > 0) && (delta[0] == 1) && (0 < G.olaps[thisOvl].a_hang)) {
        int32  stop = min(deltaLen, (int32)G.olaps[thisOvl].a_hang);  //  a_hang is int32:31!
        int32  i = 0;

        for  (i=0; (i < stop) && (delta[i] == 1); i++)
          ;

        assert((i == stop) || (delta[i] != -1));

        deltaLen -= i;

        memmove(delta, delta + i, deltaLen * sizeof (int));

        a_part     += i;
        a_end      -= i;
        a_part_len -= i;
        errors     -= i;

      } else if ((deltaLen > 0) && (delta[0] == -1) && (G.olaps[thisOvl].a_hang < 0)) {
        int32  stop = min(deltaLen, - G.olaps[thisOvl].a_hang);
        int32  i = 0;

        for  (i=0; (i < stop) && (delta[i] == -1); i++)
          ;

        assert((i == stop) || (delta[i] != 1));

        deltaLen -= i;

        memmove(delta, delta + i, deltaLen * sizeof (int));

        b_part     += i;
        b_end      -= i;
        b_part_len -= i;
        errors     -= i;
      }

      Total_Alignments_Ct++;

      int32  olapLen = min(a_end, b_end);

      if ((match_to_end == false) || (olapLen <= 0)) {
        Failed_Alignments_Ct++;

        fprintf(stderr, "Redo_Olaps()--  Bad alignment  a_iid %u  b_iid %u  errors %d  a_end %d  b_end %d\n",
                G.olaps[thisOvl].a_iid, G.olaps[thisOvl].b_iid, errors, a_end, b_end);

        continue;
      }

      //Display_Alignment(a_part, a_part_len, b_part, b_part_len, delta, deltaLen);

      G.olaps[thisOvl].evalue = AS_OVS_encodeQuality(errors / olapLen);

      //  Ancient, would output an overlap _message_ if the fp was open.  This is the -o option, which runCA doesn't use.
      //if ((quality <= Quality_Threshold) ||
      //    (G.olaps[thisOvl].a_hang <= 0 && Frag[sub].keep_left) ||
      //    (G.olaps[thisOvl].b_hang >= 0 && Frag[sub].keep_right))
      //  Output_OVL(olap, quality);
    }
  }

  fprintf(stderr, "Total:  "F_U64"\n", Total_Alignments_Ct);
  fprintf(stderr, "Failed: "F_U64"\n", Failed_Alignments_Ct);
}



