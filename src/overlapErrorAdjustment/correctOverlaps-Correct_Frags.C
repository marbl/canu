
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
 *    Brian P. Walenz from 2015-MAY-20 to 2015-JUN-03
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-MAY-02
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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
            uint64                Clen,
            uint64               *changes) {

#if DO_NO_CORRECTIONS
  //  for testing if the adjustments are screwed up.  yup.
  strcpy(fseq, oseq);
  fseqLen += oseqLen;
  return;
#endif

  //fprintf(stderr, "Correcting read %u\n", curID);

  //  Find the correct corrections.

  while ((Cpos < Clen) && (C[Cpos].readID < curID)) {
    //fprintf(stderr, "SKIP Cpos=%d for read %u, want read %u\n", Cpos, C[Cpos].readID, curID);
    Cpos++;
  }

  //  Skip any IDENT message.

  assert(C[Cpos].type == IDENT);

  //G.reads[G.readsLen].keep_left  = C[Cpos].keep_left;
  //G.reads[G.readsLen].keep_right = C[Cpos].keep_right;

  Cpos++;

  //fprintf(stderr, "Start at Cpos=%d position=%d type=%d id=%d\n", Cpos, C[Cpos].pos, C[Cpos].type, C[Cpos].readID);

  int32   adjVal = 0;

  for (uint32 i=0; i<oseqLen; i++) {

    //  No more corrections, or no more corrections for this read -- just copy bases till the end.
    if ((Cpos == Clen) || (C[Cpos].readID != curID)) {
      //fprintf(stderr, "no more corrections at i=%u, copy rest of read as is\n", i);
      while (i < oseqLen)
        fseq[fseqLen++] = filter[oseq[i++]];
      break;
    }

    //  Not at a correction -- copy the base.
    if (i < C[Cpos].pos) {
      fseq[fseqLen++] = filter[oseq[i]];
      continue;
    }

    if ((i != C[Cpos].pos) &&
        (i != C[Cpos].pos + 1))
      fprintf(stderr, "i="F_U32" Cpos="F_U64" C[Cpos].pos="F_U32"\n", i, Cpos, C[Cpos].pos);
    assert((i == C[Cpos].pos) ||
           (i == C[Cpos].pos + 1));

    if (changes)
      changes[C[Cpos].type]++;

    switch (C[Cpos].type) {
      case DELETE:  //  Delete base
        //fprintf(stderr, "DELETE %u pos %u adjust %d\n", fadjLen, i+1, adjVal-1);
        fadj[fadjLen].adjpos = i + 1;
        fadj[fadjLen].adjust = --adjVal;
        fadjLen++;
        break;

      case A_SUBST:  fseq[fseqLen++] = 'a';  break;
      case C_SUBST:  fseq[fseqLen++] = 'c';  break;
      case G_SUBST:  fseq[fseqLen++] = 'g';  break;
      case T_SUBST:  fseq[fseqLen++] = 't';  break;

      case A_INSERT:
        if (i != C[Cpos].pos + 1) {                // Insert not immediately after subst
          //fprintf(stderr, "A i=%d != C[%d].pos+1=%d\n", i, Cpos, C[Cpos].pos+1);
          fseq[fseqLen++] = filter[oseq[i++]];
        }
        fseq[fseqLen++] = 'a';

        fadj[fadjLen].adjpos = i + 1;
        fadj[fadjLen].adjust = ++adjVal;
        fadjLen++;
        i--;  //  Undo the automagic loop increment
        break;

      case C_INSERT:
        if (i != C[Cpos].pos + 1) {
          //fprintf(stderr, "C i=%d != C[%d].pos+1=%d\n", i, Cpos, C[Cpos].pos+1);
          fseq[fseqLen++] = filter[oseq[i++]];
        }
        fseq[fseqLen++] = 'c';

        fadj[fadjLen].adjpos = i + 1;
        fadj[fadjLen].adjust = ++adjVal;
        fadjLen++;
        i--;
        break;

      case G_INSERT:
        if (i != C[Cpos].pos + 1) {
          //fprintf(stderr, "G i=%d != C[%d].pos+1=%d\n", i, Cpos, C[Cpos].pos+1);
          fseq[fseqLen++] = filter[oseq[i++]];
        }
        fseq[fseqLen++] = 'g';

        fadj[fadjLen].adjpos = i + 1;
        fadj[fadjLen].adjust = ++adjVal;
        fadjLen++;
        i--;
        break;

      case T_INSERT:
        if (i != C[Cpos].pos + 1) {
          //fprintf(stderr, "T i=%d != C[%d].pos+1=%d\n", i, Cpos, C[Cpos].pos+1);
          fseq[fseqLen++] = filter[oseq[i++]];
        }
        fseq[fseqLen++] = 't';

        fadj[fadjLen].adjpos = i + 1;
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

  fseq[fseqLen] = 0;

  //fprintf(stdout, ">%u\n%s\n", curID, fseq);
}









//  Open and read corrections from  Correct_File_Path  and
//  apply them to sequences in  Frag .

//  Load reads from gkpStore, and apply corrections.

void
Correct_Frags(coParameters *G,
              gkStore      *gkpStore) {

  //  The original converted to lowercase, and made non-acgt be 'a'.

  for (uint32 i=0; i<256; i++)
    filter[i] = 'a';

  filter['A'] = filter['a'] = 'a';
  filter['C'] = filter['c'] = 'c';
  filter['G'] = filter['g'] = 'g';
  filter['T'] = filter['t'] = 't';

  //  Open the corrections, as an array.

  memoryMappedFile     *Cfile = new memoryMappedFile(G->correctionsName);
  Correction_Output_t  *C     = (Correction_Output_t *)Cfile->get();
  uint64                Cpos  = 0;
  uint64                Clen  = Cfile->length() / sizeof(Correction_Output_t);

  uint64     firstRecord   = 0;
  uint64     currentRecord = 0;

  fprintf(stderr, "Reading "F_U64" corrections from '%s'.\n", Clen, G->correctionsName);

  //  Count the number of bases, so we can do two gigantic allocations for bases and adjustments.
  //  Adjustments are always less than the number of corrections; we could also count exactly.

  G->basesLen   = 0;
  G->adjustsLen = 0;

  for (uint32 curID=G->bgnID; curID<=G->endID; curID++) {
    gkRead *read = gkpStore->gkStore_getRead(curID);

    G->basesLen += read->gkRead_sequenceLength() + 1;
  }

  for (uint64 c=0; c<Clen; c++) {
    switch (C[c].type) {
      case DELETE:
      case A_INSERT:
      case C_INSERT:
      case G_INSERT:
      case T_INSERT:
        G->adjustsLen++;
        break;
    }
  }

  fprintf(stderr, "Correcting "F_U64" bases with "F_U64" indel adjustments.\n", G->basesLen, G->adjustsLen);

  fprintf(stderr, "--Allocate "F_U64" + "F_U64" + "F_U64" MB for bases, adjusts and reads.\n",
          (sizeof(char) * G->basesLen) >> 20,
          (sizeof(Adjust_t) * G->adjustsLen) >> 20,
          (sizeof(Frag_Info_t) * (G->endID - G->bgnID + 1)) >> 20);

  G->bases        = new char          [G->basesLen];
  G->adjusts      = new Adjust_t      [G->adjustsLen];
  G->reads        = new Frag_Info_t   [G->endID - G->bgnID + 1];
  G->readsLen     = 0;

  G->basesLen   = 0;
  G->adjustsLen = 0;

  uint64   changes[12] = {0};

  //  Load reads and apply corrections for each one.

  gkReadData *readData = new gkReadData;

  for (uint32 curID=G->bgnID; curID<=G->endID; curID++) {
    gkRead *read       = gkpStore->gkStore_getRead(curID);

    gkpStore->gkStore_loadReadData(read, readData);

    uint32  readLength = read->gkRead_sequenceLength();
    char   *readBases  = readData->gkReadData_getSequence();

    //  Save pointers to the bases and adjustments.

    G->reads[G->readsLen].bases       = G->bases   + G->basesLen;
    G->reads[G->readsLen].basesLen    = 0;
    G->reads[G->readsLen].adjusts     = G->adjusts + G->adjustsLen;
    G->reads[G->readsLen].adjustsLen  = 0;

    //  Find the correct corrections.

    while ((Cpos < Clen) && (C[Cpos].readID < curID))
      Cpos++;

    //  We should be at the IDENT message.

    if (C[Cpos].type != IDENT) {
      fprintf(stderr, "ERROR: didn't find IDENT at Cpos="F_U64" for read "F_U32"\n", Cpos, curID);
      fprintf(stderr, "       C[Cpos] = keep_left=%u keep_right=%u type=%u pos=%u readID=%u\n",
              C[Cpos].keep_left,
              C[Cpos].keep_right,
              C[Cpos].type,
              C[Cpos].pos,
              C[Cpos].readID);
    }
    assert(C[Cpos].type == IDENT);

    G->reads[G->readsLen].keep_left  = C[Cpos].keep_left;
    G->reads[G->readsLen].keep_right = C[Cpos].keep_right;

    //Cpos++;

    //  Now do the corrections.

    correctRead(curID,
                G->reads[G->readsLen].bases,
                G->reads[G->readsLen].basesLen,
                G->reads[G->readsLen].adjusts,
                G->reads[G->readsLen].adjustsLen,
                readData->gkReadData_getSequence(),
                read->gkRead_sequenceLength(),
                C,
                Cpos,
                Clen,
                changes);

    //  Update the lengths in the globals.

    G->basesLen   += G->reads[G->readsLen].basesLen   + 1;
    G->adjustsLen += G->reads[G->readsLen].adjustsLen;
    G->readsLen   += 1;
  }

  delete readData;
  delete Cfile;

  fprintf(stderr, "Corrected "F_U64" bases with "F_U64" substitutions, "F_U64" deletions and "F_U64" insertions.\n",
          G->basesLen,
          changes[A_SUBST] + changes[C_SUBST] + changes[G_SUBST] + changes[T_SUBST],
          changes[DELETE],
          changes[A_INSERT] + changes[C_INSERT] + changes[G_INSERT] + changes[T_INSERT]);
}
