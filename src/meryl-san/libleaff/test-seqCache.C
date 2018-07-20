
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

#include "util.h"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

#include "test-correctSequence.H"


uint32
testSeqVsCorrect(seqInCore *S, uint32 testID) {
  uint32 err = 0;

  if (S == 0L) {
    fprintf(stderr, "testID:"uint32FMT" - empty sequence\n", testID);
    return(1);
  }

  uint32 sid = S->getIID();

  if (strcmp(S->header(),   correctSequence[sid].header) != 0) {
    fprintf(stderr, "testID:"uint32FMT" - header differs '%s' vs '%s'\n", testID, S->header(), correctSequence[sid].header);
    err++;
  }
  if (S->headerLength() != correctSequence[sid].headerLength) {
    fprintf(stderr, "testID:"uint32FMT" - header length differs "uint32FMT" vs "uint32FMT"\n", testID, S->headerLength(), correctSequence[sid].headerLength);
    err++;
  }
  if (strcmp(S->sequence(), correctSequence[sid].sequence) != 0) {
    fprintf(stderr, "testID:"uint32FMT" - sequence differs\n", testID);
    err++;
  }
  if (strlen(S->sequence()) != correctSequence[sid].sequenceLength) {
    fprintf(stderr, "testID:"uint32FMT" - sequence length differs strlen "uint32FMT" vs "uint32FMT"\n", testID, (uint32)strlen(S->sequence()), correctSequence[sid].sequenceLength);
    err++;
  }
  if (S->sequenceLength() != correctSequence[sid].sequenceLength) {
    fprintf(stderr, "testID:"uint32FMT" - sequence length differs "uint32FMT" vs "uint32FMT"\n", testID, S->sequenceLength(), correctSequence[sid].sequenceLength);
    err++;
  }

  return(err);
}


uint32
testSeqCacheIDLookups(seqCache *SC) {
  uint32      err    = 0;
  uint32      numSeq = SC->getNumberOfSequences();
  double      start  = getTime();

  //  1 - getSequenceIID()
  fprintf(stderr, "1 - getSequenceIID()\n");
  for (uint32 sid=0; sid<numSeq; sid++) {
    if (sid != SC->getSequenceIID(correctSequence[sid].header)) {
      fprintf(stderr, "2 - failed to find name '%s'\n", correctSequence[sid].header);
      err++;
    }
  }

  fprintf(stderr, "Test took %f seconds.\n", getTime() - start);

  return(err);
}


uint32
testSeqCache(seqCache *SC) {
  uint32      err    = 0;
  uint32      numSeq = SC->getNumberOfSequences();
  seqInCore  *S      = 0L;
  double      start  = getTime();

  //  0 - getSequenceLength()
  fprintf(stderr, "0 - getSequenceLength()\n");
  for (uint32 sid=0; sid<numSeq; sid++)
    if (SC->getSequenceLength(sid) != correctSequence[sid].sequenceLength) {
      fprintf(stderr, "1 - length differs.\n");
      err++;
    }

  //  2 - stream with getSequenceInCore()
  fprintf(stderr, "2 - stream with getSequenceInCore()\n");
  S = SC->getSequenceInCore();
  while (S != 0L) {
    err += testSeqVsCorrect(S, 2);
    delete S;
    S = SC->getSequenceInCore();
  }

  //  3 - iterate with getSequenceInCore(sid++)
  fprintf(stderr, "3 - iterate with getSequenceInCore(sid++)\n");
  for (uint32 sid=0; sid<numSeq; sid++) {
    S = SC->getSequenceInCore(sid);
    err += testSeqVsCorrect(S, 3);
    delete S;
  }

  //  4 - random with getSequenceInCore(sid)
  fprintf(stderr, "4 - random with getSequenceInCore(sid)\n");
  for (uint32 cnt=0; cnt<4*numSeq; cnt++) {
    uint32 sid = mtRandom32(mtctx) % numSeq;
    S = SC->getSequenceInCore(sid);
    err += testSeqVsCorrect(S, 4);
    delete S;
  }

  fprintf(stderr, "Test took %f seconds.\n", getTime() - start);

  return(err);
}


int
main(int argc, char **argv) {
  uint32     minLen = 100;
  uint32     maxLen = 2000;
  uint32     numSeq = 100000;
  seqCache  *SC     = 0L;
  uint32     err    = 0;

  generateCorrectSequence(minLen, maxLen, numSeq);

  fprintf(stderr, "seqCache(file, 0, true) (ID lookups)\n");
  SC = new seqCache("test-correctSequence.fasta", 0, true);
  //err += testSeqCacheIDLookups(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 0, true)\n");
  SC = new seqCache("test-correctSequence.fasta", 0, true);
  err += testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 1, true)\n");
  SC = new seqCache("test-correctSequence.fasta", 1, true);
  err += testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 2, true)\n");
  SC = new seqCache("test-correctSequence.fasta", 2, true);
  err += testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 4, true)\n");
  SC = new seqCache("test-correctSequence.fasta", 4, true);
  err += testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 8, true)\n");
  SC = new seqCache("test-correctSequence.fasta", 8, true);
  err += testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 32, true)\n");
  SC = new seqCache("test-correctSequence.fasta", 32, true);
  err += testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 200, true)\n");
  SC = new seqCache("test-correctSequence.fasta", 200, true);
  err += testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 1000000, true)\n");
  SC = new seqCache("test-correctSequence.fasta", 1000000, true);
  err += testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 0, true) -- loadAllSequence\n");
  SC = new seqCache("test-correctSequence.fasta", 0, true);
  SC->loadAllSequences();
  err += testSeqCache(SC);
  delete SC;

  removeCorrectSequence(numSeq);

  if (err == 0)
    fprintf(stderr, "Success!\n");

  exit(err > 0);
}
