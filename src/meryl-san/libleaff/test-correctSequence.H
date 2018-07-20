
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

#ifndef TEST_CORRECTSEQUENCE_H
#define TEST_CORRECTSEQUENCE_H

//#define WITH_WHITESPACE

struct correctSequence_t {
  char      header[256];
  uint32    headerLength;
  char     *sequence;
  uint32    sequenceLength;
};


correctSequence_t  *correctSequence = 0L;
mt_s               *mtctx           = 0L;

char     *chainSeq;
uint32   *chainSeqPos;
uint32   *chainSeqIID;
uint64   *chainStrPos;


void
generateCorrectSequence(uint32 minLen, uint32 maxLen, uint32 numSeq) {
  char      bases[4] = {'A', 'C', 'G', 'T'};

  uint32    n = numSeq;
  uint32    s = minLen;
  uint32    l = maxLen;

  uint32    seed = (uint32)(getTime() * 1000);

  fprintf(stderr, "generateCorrectSequence()-- Using seed "uint32FMT"\n", seed);
  fprintf(stderr, "generateCorrectSequence()-- Generating "uint32FMT" sequences of length "uint32FMT" to "uint32FMT"\n", numSeq, minLen, maxLen);

  correctSequence = new correctSequence_t [n];
  mtctx           = mtInit(seed);

  FILE *F = fopen("test-correctSequence.fasta", "w");

  for (uint32 i=0; i<n; i++) {
    uint32 j = s + (mtRandom32(mtctx) % (l-s));
    uint32 p = 0;

    sprintf(correctSequence[i].header, "sequence%d", i);

    correctSequence[i].headerLength   = strlen(correctSequence[i].header);
    correctSequence[i].sequence       = new char [j+1];
    correctSequence[i].sequenceLength = j;

    while (p < j)
      correctSequence[i].sequence[p++] = bases[mtRandom32(mtctx) & 0x3];

    correctSequence[i].sequence[p] = 0;

    //  Spend lots of pain adding whitespace to the file.

#ifdef WITH_WHITESPACE
    for (uint32 r=mtRandom32(mtctx) % 4; r--; )
      fprintf(F, "\n");

    fprintf(F, ">%s\n", correctSequence[i].header);

    for (uint32 r=mtRandom32(mtctx) % 4; r--; )
      fprintf(F, "\n");

    for (uint32 p=0; p<correctSequence[i].sequenceLength; p++) {
      fprintf(F, "%c", correctSequence[i].sequence[p]);
      if ((mtRandom32(mtctx) % 100) == 0)
        for (uint32 r=mtRandom32(mtctx) % 4; r--; )
          fprintf(F, "\n");
    }
#else
    fprintf(F, ">%s\n", correctSequence[i].header);
    fprintf(F, "%s\n",  correctSequence[i].sequence);
#endif

  }

  for (uint32 r=mtRandom32(mtctx) % 4; r--; )
    fprintf(F, "\n");

  fclose(F);
}


void
generateChainedAnswer(uint32 numSeq, char sep, uint32 sepLen) {
  uint32 maxLen = 0;

  for (uint32 i=0; i<numSeq; i++)
    maxLen += correctSequence[i].sequenceLength + sepLen;

  maxLen -= sepLen;

  delete [] chainSeq;
  delete [] chainSeqPos;
  delete [] chainSeqIID;
  delete [] chainStrPos;

  chainSeq    = new char   [maxLen + 1];
  chainSeqPos = new uint32 [maxLen];
  chainSeqIID = new uint32 [maxLen];
  chainStrPos = new uint64 [maxLen];

  uint32 p      = 0;
  uint64 strpos = 0;

  for (uint32 sid=0; sid<numSeq; sid++) {
    for (uint32 ppp=0; ppp<correctSequence[sid].sequenceLength; ppp++, p++, strpos++) {
      chainSeq[p]    = correctSequence[sid].sequence[ppp];
      chainSeqPos[p] = ppp;
      chainSeqIID[p] = sid;
      chainStrPos[p] = strpos;
    }
    if (sid+1 < numSeq) {
      for (uint32 ppp=0; ppp<sepLen; ppp++, p++) {
        chainSeq[p]    = sep;
        chainSeqPos[p] = ~uint32ZERO;
        chainSeqIID[p] = ~uint32ZERO;
        chainStrPos[p] = ~uint64ZERO;
      }
    }
  }

  assert(p == maxLen);

  //  Used to get the length of the string.
  chainSeq[maxLen] = 0;
}


void
removeCorrectSequence(uint32 numSeq) {

  unlink("test-correctSequence.fasta");
  unlink("test-correctSequence.fastaidx");

  for (uint32 i=0; i<numSeq; i++)
    delete [] correctSequence[i].sequence;

  delete [] correctSequence;

  delete [] chainSeq;
  delete [] chainSeqPos;
  delete [] chainSeqIID;
  delete [] chainStrPos;
}


#endif  //  TEST_CORRECTSEQUENCE_H
