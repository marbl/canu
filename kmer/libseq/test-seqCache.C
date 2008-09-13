#include "util.h"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"



struct correctSequence_t {
  char      header[256];
  u32bit    headerLength;
  char     *sequence;
  u32bit    sequenceLength;
};

correctSequence_t  *correctSequence = 0L;
mt_s               *mtctx           = 0L;


void
generateCorrectSequence(u32bit minLen, u32bit maxLen, u32bit numSeq) {
  char      bases[4] = {'A', 'C', 'G', 'T'};

  u32bit    n = numSeq;
  u32bit    s = minLen;
  u32bit    l = maxLen;

  correctSequence = new correctSequence_t [n];
  mtctx           = mtInit(58145);

  FILE *F = fopen("test-seqCache.fasta", "w");

  for (u32bit i=0; i<n; i++) {
    u32bit j = s + (mtRandom32(mtctx) % (l-s));
    u32bit p = 0;

    sprintf(correctSequence[i].header, "sequence%d", i);

    correctSequence[i].headerLength   = strlen(correctSequence[i].header);
    correctSequence[i].sequence       = new char [j+1];
    correctSequence[i].sequenceLength = j;

    while (p < j)
      correctSequence[i].sequence[p++] = bases[mtRandom32(mtctx) & 0x3];            

    correctSequence[i].sequence[p] = 0;

    fprintf(F, ">%s\n", correctSequence[i].header);
    fprintf(F, "%s\n",  correctSequence[i].sequence);
  }

  fclose(F);
}



u32bit
testSeqVsCorrect(seqInCore *S, u32bit testID) {
  u32bit err = 0;

  if (S == 0L) {
    fprintf(stderr, "testID:"u32bitFMT" - empty sequence\n", testID);
    //return(1);
  }

  u32bit sid = S->getIID();

  if (strcmp(S->header(),   correctSequence[sid].header) != 0) {
    fprintf(stderr, "testID:"u32bitFMT" - header differs '%s' vs '%s'\n", testID, S->header(), correctSequence[sid].header);
    err++;
  }
  if (S->headerLength() != correctSequence[sid].headerLength) {
    fprintf(stderr, "testID:"u32bitFMT" - header length differs "u32bitFMT" vs "u32bitFMT"\n", testID, S->headerLength(), correctSequence[sid].headerLength);
    err++;
  }
  if (strcmp(S->sequence(), correctSequence[sid].sequence) != 0) {
    fprintf(stderr, "testID:"u32bitFMT" - sequence differs\n", testID);
    err++;
  }
  if (S->sequenceLength() != correctSequence[sid].sequenceLength) {
    fprintf(stderr, "testID:"u32bitFMT" - sequence length differs "u32bitFMT" vs "u32bitFMT"\n", testID, S->sequenceLength(), correctSequence[sid].sequenceLength);
    err++;
  }

  return(err);
}


void
testSeqCache(seqCache *SC) {
  u32bit      err    = 0;
  u32bit      numSeq = SC->getNumberOfSequences();
  seqInCore  *S      = 0L;
  double      start  = getTime();

  //  0 - getSequenceLength()
  for (u32bit sid=0; sid<numSeq; sid++)
    if (SC->getSequenceLength(sid) != correctSequence[sid].sequenceLength) {
      fprintf(stderr, "1 - length differs.\n");
      err++;
    }

  //  1 - getSequenceIID()
  for (u32bit sid=0; sid<numSeq; sid++) {
    if (sid != SC->getSequenceIID(correctSequence[sid].header)) {
      fprintf(stderr, "2 - failed to find name '%s'\n", correctSequence[sid].header);
      err++;
    }
  }

  //  2 - stream with getSequenceInCore()
  S = SC->getSequenceInCore();
  while (S != 0L) {
    err += testSeqVsCorrect(S, 2);
    delete S;
    S = SC->getSequenceInCore();
  }

  //  3 - iterate with getSequenceInCore(sid++)
  for (u32bit sid=0; sid<numSeq; sid++) {
    S = SC->getSequenceInCore(sid);
    err += testSeqVsCorrect(S, 3);
    delete S;
  }

  //  4 - random with getSequenceInCore(sid)
  for (u32bit cnt=0; cnt<4*numSeq; cnt++) {
    u32bit sid = mtRandom32(mtctx) % numSeq;
    S = SC->getSequenceInCore(sid);
    err += testSeqVsCorrect(S, 4);
    delete S;
  }

  fprintf(stderr, "Test took %f seconds.\n", getTime() - start);
}




int
main(int argc, char **argv) {
  u32bit     minLen = 100;
  u32bit     maxLen = 2000;
  u32bit     numSeq = 1000000;
  seqCache  *SC     = 0L;


  generateCorrectSequence(minLen, maxLen, numSeq);


  fprintf(stderr, "seqCache(file, 0, true)\n");
  SC = new seqCache("test-seqCache.fasta", 0, true);
  testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 1, true)\n");
  SC = new seqCache("test-seqCache.fasta", 1, true);
  testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 2, true)\n");
  SC = new seqCache("test-seqCache.fasta", 2, true);
  testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 4, true)\n");
  SC = new seqCache("test-seqCache.fasta", 4, true);
  testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 8, true)\n");
  SC = new seqCache("test-seqCache.fasta", 8, true);
  testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 32, true)\n");
  SC = new seqCache("test-seqCache.fasta", 32, true);
  testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 200, true)\n");
  SC = new seqCache("test-seqCache.fasta", 200, true);
  testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 1000000, true)\n");
  SC = new seqCache("test-seqCache.fasta", 1000000, true);
  testSeqCache(SC);
  delete SC;

  fprintf(stderr, "seqCache(file, 0, true) -- loadAllSequence\n");
  SC = new seqCache("test-seqCache.fasta", 0, true);
  SC->loadAllSequences();
  testSeqCache(SC);
  delete SC;




  for (u32bit i=0; i<numSeq; i++)
    delete [] correctSequence[i].sequence;

  delete [] correctSequence;
}
