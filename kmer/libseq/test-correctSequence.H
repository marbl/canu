#ifndef TEST_CORRECTSEQUENCE_H
#define TEST_CORRECTSEQUENCE_H

//#define WITH_WHITESPACE

struct correctSequence_t {
  char      header[256];
  u32bit    headerLength;
  char     *sequence;
  u32bit    sequenceLength;
};


correctSequence_t  *correctSequence = 0L;
mt_s               *mtctx           = 0L;

char     *chainSeq;
u32bit   *chainSeqPos;
u32bit   *chainSeqIID;
u64bit   *chainStrPos;


void
generateCorrectSequence(u32bit minLen, u32bit maxLen, u32bit numSeq) {
  char      bases[4] = {'A', 'C', 'G', 'T'};

  u32bit    n = numSeq;
  u32bit    s = minLen;
  u32bit    l = maxLen;

  u32bit    seed = (u32bit)(getTime() * 1000);

  fprintf(stderr, "generateCorrectSequence()-- Using seed "u32bitFMT"\n", seed);
  fprintf(stderr, "generateCorrectSequence()-- Generating "u32bitFMT" sequences of length "u32bitFMT" to "u32bitFMT"\n", numSeq, minLen, maxLen);

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

    //  Spend lots of pain adding whitespace to the file.

#ifdef WITH_WHITESPACE
    for (u32bit r=mtRandom32(mtctx) % 4; r--; )
      fprintf(F, "\n");

    fprintf(F, ">%s\n", correctSequence[i].header);

    for (u32bit r=mtRandom32(mtctx) % 4; r--; )
      fprintf(F, "\n");

    for (u32bit p=0; p<correctSequence[i].sequenceLength; p++) {
      fprintf(F, "%c", correctSequence[i].sequence[p]);
      if ((mtRandom32(mtctx) % 100) == 0)
        for (u32bit r=mtRandom32(mtctx) % 4; r--; )
          fprintf(F, "\n");
    }
#else
    fprintf(F, ">%s\n", correctSequence[i].header);
    fprintf(F, "%s\n",  correctSequence[i].sequence);
#endif

  }

  for (u32bit r=mtRandom32(mtctx) % 4; r--; )
    fprintf(F, "\n");

  fclose(F);
}


void
generateChainedAnswer(u32bit numSeq, char sep, u32bit sepLen) {
  u32bit maxLen = 0;

  for (u32bit i=0; i<numSeq; i++)
    maxLen += correctSequence[i].sequenceLength + sepLen;

  maxLen -= sepLen;

  delete [] chainSeq;
  delete [] chainSeqPos;
  delete [] chainSeqIID;
  delete [] chainStrPos;

  chainSeq    = new char   [maxLen + 1];
  chainSeqPos = new u32bit [maxLen];
  chainSeqIID = new u32bit [maxLen];
  chainStrPos = new u64bit [maxLen];

  u32bit p      = 0;
  u64bit strpos = 0;

  for (u32bit sid=0; sid<numSeq; sid++) {
    for (u32bit ppp=0; ppp<correctSequence[sid].sequenceLength; ppp++, p++, strpos++) {
      chainSeq[p]    = correctSequence[sid].sequence[ppp];
      chainSeqPos[p] = ppp;
      chainSeqIID[p] = sid;
      chainStrPos[p] = strpos;
    }
    if (sid+1 < numSeq) {
      for (u32bit ppp=0; ppp<sepLen; ppp++, p++) {
        chainSeq[p]    = sep;
        chainSeqPos[p] = ~u32bitZERO;
        chainSeqIID[p] = ~u32bitZERO;
        chainStrPos[p] = ~u64bitZERO;
      }
    }
  }

  assert(p == maxLen);

  //  Used to get the length of the string.
  chainSeq[maxLen] = 0;
}


void
removeCorrectSequence(u32bit numSeq) {

  unlink("test-correctSequence.fasta");
  unlink("test-correctSequence.fastaidx");

  for (u32bit i=0; i<numSeq; i++)
    delete [] correctSequence[i].sequence;

  delete [] correctSequence;

  delete [] chainSeq;
  delete [] chainSeqPos;
  delete [] chainSeqIID;
  delete [] chainStrPos;
}


#endif  //  TEST_CORRECTSEQUENCE_H
