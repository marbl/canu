
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
 *  This file is derived from:
 *
 *    kmer/leaff/stats.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2009-FEB-07 to 2014-APR-11
 *      are Copyright 2009,2012-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-08 to 2015-MAR-21
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "seqCache.H"

#include <algorithm>

using namespace std;


void
stats(char *filename, uint64 refLen) {
  seqCache    *F = new seqCache(filename);

  bool                  V[256];
  for (uint32 i=0; i<256; i++)
    V[i] = false;
  V['n'] = true;
  V['N'] = true;

  uint32  numSeq = F->getNumberOfSequences();

  uint64    Ss = 0;  //  actual length of span
  uint64    Rs = 0;  //  reference length of span
  uint32   *Ls = new uint32 [numSeq];

  uint64    Sb = 0;
  uint64    Rb = 0;
  uint32   *Lb = new uint32 [numSeq];

  for (uint32 i=0; i<numSeq; i++)
    Ls[i] = Lb[i] = 0;

  for (uint32 s=0; s<numSeq; s++) {
    seqInCore  *S      = F->getSequenceInCore(s);
    uint32      len    = S->sequenceLength();
    uint32      span   = len;
    uint32      base   = len;

    for (uint32 pos=1; pos<len; pos++) {
      if (V[S->sequence()[pos]])
        base--;
    }

    Ss += span;
    Sb += base;

    Ls[S->getIID()] = span;
    Lb[S->getIID()] = base;

    delete S;
  }

  if (refLen > 0) {
    Rs = refLen;
    Rb = refLen;
  } else {
    Rs = Ss;
    Rb = Sb;
  }

  //qsort(Ls, numSeq, sizeof(uint32), uint32_compare);
  //qsort(Lb, numSeq, sizeof(uint32), uint32_compare);

  sort(Ls, Ls + numSeq);
  sort(Lb, Lb + numSeq);

  reverse(Ls, Ls + numSeq);
  reverse(Lb, Lb + numSeq);

  uint32  n50s[11] = {0};
  uint32  l50s[11] = {0};

  uint32  n50b[11] = {0};
  uint32  l50b[11] = {0};

  uint32  sizes[11] = {0};
  uint32  sizeb[11] = {0};

  for (uint32 i=0; i<11; i++) {
    sizes[i] = i * Rs / 10;
    sizeb[i] = i * Rb / 10;
    //fprintf(stderr, "SIZE %2d  s=%d b=%d\n", i, sizes[i], sizeb[i]);
  }

  for (uint32 i=0, sum=0, n=1; (i < numSeq) && (n < 11); i++) {
    if ((sum <  sizes[n]) && (sizes[n] <= sum + Ls[i])) {
      n50s[n]  = Ls[i];
      l50s[n]  = i;
      n++;
    }

    sum += Ls[i];
  }


  for (uint32 i=0, sum=0, n=1; (i < numSeq) && (n < 11); i++) {
    if ((sum <  sizeb[n]) && (sizeb[n] <= sum + Lb[i])) {
      n50b[n]  = Ls[i];
      l50b[n]  = i;
      n++;
    }

    sum += Lb[i];
  }

  //for (uint32 i=0, sum=0; sum < Rb/2; i++) {
  //}

  fprintf(stdout, "%s\n", F->getSourceName());
  fprintf(stdout, "\n");
  fprintf(stdout, "numSeqs  "F_U32"\n", numSeq);
  fprintf(stdout, "\n");
  fprintf(stdout, "SPAN (smallest "F_U32" largest "F_U32")\n", Ls[numSeq-1], Ls[0]);
  for (uint32 i=1; i<10; i++)
    fprintf(stdout, "n"F_U32"    %10"F_U32P" at index "F_U32"\n", 10 * i, n50s[i], l50s[i]);
  fprintf(stdout, "totLen %10"F_U64P"\n", Ss);
  fprintf(stdout, "refLen %10"F_U64P"\n", Rs);
  fprintf(stdout, "\n");
  fprintf(stdout, "BASES (smallest "F_U32" largest "F_U32")\n", Lb[numSeq-1], Lb[0]);
  for (uint32 i=1; i<10; i++)
    fprintf(stdout, "n"F_U32"    %10"F_U32P" at index "F_U32"\n", 10 * i, n50b[i], l50b[i]);
  fprintf(stdout, "totLen %10"F_U64P"\n", Sb);
  fprintf(stdout, "refLen %10"F_U64P"\n", Rb);

  delete [] Ls;
  delete [] Lb;
}
