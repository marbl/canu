
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


//  A filter to convert upper- to lower-case, and to map all undefined bases to 'a'.
//
//    for (uint32 i = 0; i < 256; i++)
//      filter[i] = 'a';
//
//    filter['A'] = filter['a'] = 'a';
//    filter['C'] = filter['c'] = 'c';
//    filter['G'] = filter['g'] = 'g';
//    filter['T'] = filter['t'] = 't';
//
char const
filter[256] = { 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'c', 'a', 'a', 'a', 'g', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 't', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'c', 'a', 'a', 'a', 'g', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 't', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
                'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a' };



void
correctRead(uint32 curID,
            char *fseq, uint32 &fseqLen, Adjust_t *fadj, uint32 &fadjLen,
            const char *oseq, uint32  oseqLen,
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
    //fprintf(stderr, "SKIP Cpos=%u Clen=%u for read %u, want read %u\n", Cpos, Clen, C[Cpos].readID, curID);
    Cpos++;
  }

  //  Skip any IDENT message.

  assert(C[Cpos].type == IDENT);

  //G.reads[G.readsLen].keep_left  = C[Cpos].keep_left;
  //G.reads[G.readsLen].keep_right = C[Cpos].keep_right;

  Cpos++;
  assert(Cpos <= Clen);

  //fprintf(stderr, "Start at Cpos=%d position=%d type=%d id=%d\n", Cpos, C[Cpos].pos, C[Cpos].type, C[Cpos].readID);

  int32   adjVal = 0;

  for (uint32 i = 0; i < oseqLen; ) {

    //  No more corrections OR no more corrections for this read OR no correction at position -- just copy base
    if (Cpos == Clen || C[Cpos].readID != curID || i < C[Cpos].pos) {
      //fprintf(stderr, "Introducing IDENT '%c' read=%u i=%u \n", filter[oseq[i]], C[Cpos].readID, i);
      fseq[fseqLen++] = filter[oseq[i++]];
      continue;
    }

    assert(Cpos < Clen);
    assert(i == C[Cpos].pos);

    if (changes)
      changes[C[Cpos].type]++;

    switch (C[Cpos].type) {
      case DELETE:  //  Delete base
        //fprintf(stderr, "DELETE %u pos %u adjust %d\n", fadjLen, i+1, adjVal-1);
        //fprintf(stderr, "Introducing DELETION read=%u i=%u \n", C[Cpos].readID, i);
        fadj[fadjLen].adjpos = i + 1;
        fadj[fadjLen].adjust = --adjVal;
        fadjLen++;
        i++;
        break;

      case A_SUBST:
      case C_SUBST:
      case G_SUBST:
      case T_SUBST:
        //fprintf(stderr, "Introducing SUBST '%c' -> '%c' read=%u i=%u \n", filter[oseq[i]], VoteChar(C[Cpos].type), C[Cpos].readID, i);
        fseq[fseqLen++] = VoteChar(C[Cpos].type);
        i++;
        break;

      case A_INSERT:
      case C_INSERT:
      case G_INSERT:
      case T_INSERT:
        //fprintf(stderr, "Introducing INSERTION '%c' read=%u i=%u \n", VoteChar(C[Cpos].type), C[Cpos].readID, i);
        fseq[fseqLen++] = VoteChar(C[Cpos].type);

        fadj[fadjLen].adjpos = i + 1;
        fadj[fadjLen].adjust = ++adjVal;
        fadjLen++;
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

//  Load reads from seqStore, and apply corrections.

void
Correct_Frags(coParameters *G,
              sqStore      *seqStore,
              FILE *correctedReads) {

  //  Open the corrections, as an array.

  memoryMappedFile     *Cfile = new memoryMappedFile(G->correctionsName);
  Correction_Output_t  *C     = (Correction_Output_t *)Cfile->get();
  uint64                Cpos  = 0;
  uint64                Clen  = Cfile->length() / sizeof(Correction_Output_t);

  uint64     firstRecord   = 0;
  uint64     currentRecord = 0;

  fprintf(stderr, "Reading " F_U64 " corrections from '%s'.\n", Clen, G->correctionsName);

  //  Count the number of bases, so we can do two gigantic allocations for bases and adjustments.
  //  Adjustments are always less than the number of corrections; we could also count exactly.

  G->basesLen   = 0;

  for (uint32 curID=G->bgnID; curID<=G->endID; curID++)
    G->basesLen += seqStore->sqStore_getReadLength(curID) + 1;

  uint64 del_cnt = 0;
  uint64 ins_cnt = 0;
  for (uint64 c = 0; c < Clen; c++) {
    switch (C[c].type) {
      case DELETE:
        del_cnt++;
        break;
      case A_INSERT:
      case C_INSERT:
      case G_INSERT:
      case T_INSERT:
        ins_cnt++;
        break;
      default: {}
    }
  }
  G->basesLen += ins_cnt; // allow extra space for insertions in case reads get longer
  G->adjustsLen = ins_cnt + del_cnt;

  fprintf(stderr, "Correcting " F_U64 " bases with " F_U64 " indel adjustments.\n", G->basesLen, G->adjustsLen);

  fprintf(stderr, "--Allocate " F_U64 " + " F_U64 " + " F_U64 " MB for bases, adjusts and reads.\n",
          (sizeof(char)        * (uint64)(G->basesLen))             / 1048576,   //  MacOS GCC 4.9.4 can't decide if these three
          (sizeof(Adjust_t)    * (uint64)(G->adjustsLen))           / 1048576,   //  values are %u, %lu or %llu.  We force cast
          (sizeof(Frag_Info_t) * (uint64)(G->endID - G->bgnID + 1)) / 1048576);  //  them to be uint64.

  G->bases        = new char          [G->basesLen];
  G->adjusts      = new Adjust_t      [G->adjustsLen];
  G->reads        = new Frag_Info_t   [G->endID - G->bgnID + 1];
  G->readsLen     = 0;

  G->basesLen   = 0;
  G->adjustsLen = 0;

  uint64   changes[12] = {0};

  //  Load reads and apply corrections for each one.

  for (uint32 curID=G->bgnID; curID<=G->endID; curID++) {

    auto &read = G->reads[G->readsLen];
    //  Save pointers to the bases and adjustments.
    read.bases       = G->bases   + G->basesLen;
    read.basesLen    = 0;
    read.adjusts     = G->adjusts + G->adjustsLen;
    read.adjustsLen  = 0;

    //  Find the correct corrections.

    while ((Cpos < Clen) && (C[Cpos].readID < curID))
      Cpos++;

    //  We should be at the IDENT message.

    if (C[Cpos].type != IDENT) {
      fprintf(stderr, "ERROR: didn't find IDENT at Cpos=" F_U64 " for read " F_U32 "\n", Cpos, curID);
      fprintf(stderr, "       C[Cpos] = keep_left=%u keep_right=%u type=%u pos=%u readID=%u\n",
              C[Cpos].keep_left,
              C[Cpos].keep_right,
              C[Cpos].type,
              C[Cpos].pos,
              C[Cpos].readID);
    }
    assert(C[Cpos].type == IDENT);

    read.keep_left  = C[Cpos].keep_left;
    read.keep_right = C[Cpos].keep_right;

    //Cpos++;

    //  Now actually load the read and do the corrections.

    if (seqStore->sqStore_getReadLength(curID) > 0) {
      sqRead stored;
      seqStore->sqStore_getRead(curID, &stored);

      correctRead(curID,
                  read.bases,
                  read.basesLen,
                  read.adjusts,
                  read.adjustsLen,
                  stored.sqRead_sequence(),
                  stored.sqRead_length(),
                  C,
                  Cpos,
                  Clen,
                  changes);

      if (correctedReads != NULL) {
        AS_UTL_writeFastA(correctedReads, read.bases, read.basesLen, 60, ">%d\n", curID);
      }
    }

    //  Update the lengths in the globals.

    G->basesLen   += G->reads[G->readsLen].basesLen   + 1;
    G->adjustsLen += G->reads[G->readsLen].adjustsLen;
    G->readsLen   += 1;
  }

  delete Cfile;

  fprintf(stderr, "Corrected " F_U64 " bases with " F_U64 " substitutions, " F_U64 " deletions and " F_U64 " insertions.\n",
          G->basesLen,
          changes[A_SUBST] + changes[C_SUBST] + changes[G_SUBST] + changes[T_SUBST],
          changes[DELETE],
          changes[A_INSERT] + changes[C_INSERT] + changes[G_INSERT] + changes[T_INSERT]);
}
