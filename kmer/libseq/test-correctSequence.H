#ifndef TEST_CORRECTSEQUENCE_H
#define TEST_CORRECTSEQUENCE_H


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

  u32bit    seed = (u32bit)(getTime() * 1000);

  fprintf(stderr, "generateCorrectSequence()-- Using seed "u32bitFMT"\n", seed);

  correctSequence = new correctSequence_t [n];
  mtctx           = mtInit(seed);

  FILE *F = fopen("test-correctSequence.fasta", "w");

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


void
removeCorrectSequence(u32bit numSeq) {

  unlink("test-correctSequence.fasta");
  unlink("test-correctSequence.fastaidx");

  for (u32bit i=0; i<numSeq; i++)
    delete [] correctSequence[i].sequence;

  delete [] correctSequence;
}


#endif  //  TEST_CORRECTSEQUENCE_H
