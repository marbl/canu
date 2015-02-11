

#include "correctOverlaps.H"
#include "correctionOutput.H"


//  Shouldn't be global.
char  filter[256];


void
correctRead(uint32 curID,
            char *fseq, uint32 &fseqLen, Adjust_t *fadj, uint32 &fadjLen,
            char *oseq, uint32  oseqLen,
            Correction_Output_t  *C,
            uint64               &Cpos,
            uint64                Clen) {

  //  Find the correct corrections.

  while ((Cpos < Clen) && (C[Cpos].readID < curID))
    Cpos++;

  //  Skip any IDENT message.

  assert(C[Cpos].type == IDENT);

  //G.reads[G.readsLen].keep_left  = C[Cpos].keep_left;
  //G.reads[G.readsLen].keep_right = C[Cpos].keep_right;

  Cpos++;

  int32   adjVal = 0;

  for (uint32 i=0; i<oseqLen; i++) {

    //  No more corrections, just copy bases till the end.
    if ((Cpos == Clen) || (C[Cpos].readID != curID)) {
      fseq[fseqLen++] = filter[oseq[i]];
      continue;
    }

    assert((i == C[Cpos].pos) ||
           (i == C[Cpos].pos + 1));
      
    switch (C[Cpos].type) {

    case DELETE:  //  Delete base
      fadj[fadjLen].pos    = i + 1;
      fadj[fadjLen].adjust = --adjVal;
      fadjLen++;
      break;

    case A_SUBST:  fseq[fseqLen++] = 'a';  break;
    case C_SUBST:  fseq[fseqLen++] = 'c';  break;
    case G_SUBST:  fseq[fseqLen++] = 'g';  break;
    case T_SUBST:  fseq[fseqLen++] = 't';  break;

    case A_INSERT:
      if (i != C[Cpos].pos + 1)                             // Insert not immediately after subst
        fseq[fseqLen++] = filter[oseq[i++]];     //  i++ to undo the 'undo' below
      fseq[fseqLen++] = 'a';

      fadj[fadjLen].pos    = i + 1;
      fadj[fadjLen].adjust = ++adjVal;
      fadjLen++;
      i--;  //  Undo the automagic loop increment
      break;

    case C_INSERT:
      if (i != C[Cpos].pos + 1)
        fseq[fseqLen++] = filter[oseq[i++]];
      fseq[fseqLen++] = 'c';

      fadj[fadjLen].pos    = i + 1;
      fadj[fadjLen].adjust = ++adjVal;
      fadjLen++;
      i--;
      break;

    case G_INSERT:
      if (i != C[Cpos].pos + 1)
        fseq[fseqLen++] = filter[oseq[i++]];
      fseq[fseqLen++] = 'g';

      fadj[fadjLen].pos    = i + 1;
      fadj[fadjLen].adjust = ++adjVal;
      fadjLen++;
      i--;
      break;

    case T_INSERT:
      if (i != C[Cpos].pos + 1)
        fseq[fseqLen++] = filter[oseq[i++]];
      fseq[fseqLen++] = 't';

      fadj[fadjLen].pos    = i + 1;
      fadj[fadjLen].adjust = ++adjVal;
      fadjLen++;
      i--;
      break;

    default:
      fprintf (stderr, "ERROR:  Illegal vote type\n");
      break;
    }

    Cpos++;
  }

  //  Terminate the sequence.

  fseq[fseqLen++] = 0;
}









//  Open and read corrections from  Correct_File_Path  and
//  apply them to sequences in  Frag .

//  Load reads from gkpStore, and apply corrections.

void
Correct_Frags(coParameters &G,
              gkStore      *gkpStore) {

  //  The original converted to lowercase, and made non-acgt be 'a'.

  for (uint32 i=0; i<256; i++)
    filter[i] = 'a';

  filter['A'] = filter['a'] = 'a';
  filter['C'] = filter['c'] = 'c';
  filter['G'] = filter['g'] = 'g';
  filter['T'] = filter['t'] = 't';

  //  Open the corrections, as an array.

  memoryMappedFile     *Cfile = new memoryMappedFile(G.correctionsName);
  Correction_Output_t  *C     = (Correction_Output_t *)Cfile->get();
  uint64                Cpos  = 0;
  uint64                Clen  = Cfile->length() / sizeof(Correction_Output_t);

  uint64     firstRecord   = 0;
  uint64     currentRecord = 0;

  //  Count the number of bases, so we can do two gigantic allocations for bases and adjustments.
  //  Adjustments are always less than the number of corrections; we could also count exactly.

  G.basesLen   = 0;
  G.adjustsLen = 0;

  for (uint32 curID=G.bgnID; curID<G.endID; curID++) {
    gkRead *read = gkpStore->gkStore_getRead(curID);

    G.basesLen += read->gkRead_clearRegionLength() + 1;
  }

  for (uint64 c=0; c<Clen; c++) {
    switch (C[c].type) {
    case DELETE:
    case A_INSERT:
    case C_INSERT:
    case G_INSERT:
    case T_INSERT:
      G.adjustsLen++;
      break;
    }
  }

  G.bases        = new char          [G.basesLen];
  G.adjusts      = new Adjust_t      [G.adjustsLen];
  G.reads        = new Frag_Info_t   [G.endID - G.bgnID + 1];
  G.readsLen     = 0;

  G.basesLen   = 0;
  G.adjustsLen = 0;

  //  Load reads and apply corrections for each one.

  gkReadData  readData;

  for (uint32 curID=G.bgnID; curID<G.endID; curID++) {
    gkRead *read       = gkpStore->gkStore_getRead(curID);

    if (read->gkRead_isDeleted() == true)
      continue;

    gkpStore->gkStore_loadReadData(read, &readData);

    uint32  readLength = read->gkRead_clearRegionLength();
    char   *readBases  = readData.gkReadData_getSequence();

    //  Save pointers to the bases and adjustments.

    G.reads[G.readsLen].bases       = G.bases   + G.basesLen;
    G.reads[G.readsLen].basesLen    = 0;
    G.reads[G.readsLen].adjusts     = G.adjusts + G.adjustsLen;
    G.reads[G.readsLen].adjustsLen  = 0;

    //  Find the correct corrections.

    while ((Cpos < Clen) && (C[Cpos].readID < curID))
      Cpos++;

    //  We should be at the IDENT message.

    assert(C[Cpos].type == IDENT);

    G.reads[G.readsLen].keep_left  = C[Cpos].keep_left;
    G.reads[G.readsLen].keep_right = C[Cpos].keep_right;

    Cpos++;

    //  Now do the corrections.

    correctRead(curID, 
                G.reads[G.readsLen].bases,
                G.reads[G.readsLen].basesLen,
                G.reads[G.readsLen].adjusts,
                G.reads[G.readsLen].adjustsLen,
                readData.gkReadData_getSequence(),
                read->gkRead_clearRegionLength(),
                C,
                Cpos,
                Clen);

    //  Update the lengths in the globals.

    G.basesLen   += G.reads[G.readsLen].basesLen;
    G.adjustsLen += G.reads[G.readsLen].adjustsLen;
    G.readsLen   += 1;
  }

  delete Cfile;
}
