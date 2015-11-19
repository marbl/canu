
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
 *    kmer/leaff/dups.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2009-FEB-07 to 2014-APR-11
 *      are Copyright 2009,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "seqCache.H"

#include "md5.H"

md5_s *
computeMD5ForEachSequence(seqCache *F) {
  uint32   numSeqs = F->getNumberOfSequences();
  md5_s   *result  = new md5_s [numSeqs];

  for (uint32 idx=0; idx < numSeqs; idx++) {
    seqInCore *s1 = F->getSequenceInCore(idx);
    md5_string(result+idx, s1->sequence(), s1->sequenceLength());
    result[idx].i = s1->getIID();
    delete s1;
  }

  return(result);
}


void
mapDuplicates_Print(char *filea, seqInCore *sa,
                    char *fileb, seqInCore *sb) {

  if (strcmp(sa->sequence(), sb->sequence()) == 0)
    fprintf(stdout, F_U32" <-> "F_U32"\n", sa->getIID(), sb->getIID());
  else
    fprintf(stderr, "COLLISION DETECTED BETWEEN %s:"F_U32" AND %s:"F_U32"!\nPLEASE REPORT THIS TO bri@walenz.org!\n",
            filea, sa->getIID(), fileb, sb->getIID());
}



void
findDuplicates(char *filename) {
  seqInCore  *s1 = 0L;
  seqInCore  *s2 = 0L;
  seqCache   *A = new seqCache(filename);

  uint32 numSeqs = A->getNumberOfSequences();

  fprintf(stderr, "Computing MD5's for each sequence in '%s'.\n", filename);
  md5_s *result = computeMD5ForEachSequence(A);

  fprintf(stderr, "Sorting MD5's.\n");
  qsort(result, numSeqs, sizeof(md5_s), md5_compare);

  fprintf(stderr, "Verifying identity, and output\n");
  for (uint32 idx=1; idx<numSeqs; idx++) {
    if (md5_compare(result+idx-1, result+idx) == 0) {
      if (result[idx-1].i == result[idx].i) {
        fprintf(stderr, "Internal error: found two copies of the same sequence iid ("F_U32")!\n", result[idx].i);
        exit(1);
      }

      s1 = A->getSequenceInCore(result[idx-1].i);
      s2 = A->getSequenceInCore(result[idx].i);

      if (strcmp(s1->sequence(), s2->sequence()) == 0) {
        fprintf(stdout, F_U32":%s\n"F_U32":%s\n\n",
                result[idx-1].i, s1->header(),
                result[idx  ].i, s2->header());
      } else {
        fprintf(stderr, "COLLISION DETECTED BETWEEN IID "F_U32" AND "F_U32"!\nPLEASE REPORT THIS TO bri@walenz.org!\n",
                result[idx-1].i, result[idx].i);
      }

      delete s1;
      delete s2;
    }
  }

  delete [] result;
  delete    A;
}



void
mapDuplicates(char *filea, char *fileb) {
  fprintf(stderr, "Computing MD5's for each sequence in '%s'.\n", filea);
  seqCache  *A = new seqCache(filea);
  md5_s     *resultA = computeMD5ForEachSequence(A);

  fprintf(stderr, "Computing MD5's for each sequence in '%s'.\n", fileb);
  seqCache  *B = new seqCache(fileb);
  md5_s     *resultB = computeMD5ForEachSequence(B);

  uint32  numSeqsA = A->getNumberOfSequences();
  uint32  numSeqsB = B->getNumberOfSequences();
  uint32  idxA = 0;
  uint32  idxB = 0;

  fprintf(stderr, "Sorting MD5's.\n");
  qsort(resultA, numSeqsA, sizeof(md5_s), md5_compare);
  qsort(resultB, numSeqsB, sizeof(md5_s), md5_compare);

  fprintf(stderr, "Finding duplicates.\n");
  while ((idxA<numSeqsA) && (idxB<numSeqsB)) {
    int res = md5_compare(resultA+idxA, resultB+idxB);

    if (res == 0) {
      seqInCore *sa = A->getSequenceInCore(resultA[idxA].i);
      seqInCore *sb = B->getSequenceInCore(resultB[idxB].i);

      mapDuplicates_Print(filea, sa, fileb, sb);

      //  While the B sequence matches the current A sequence, output a match
      //
      uint32 idxBb = idxB+1;
      int resb = md5_compare(resultA+idxA, resultB+idxBb);
      while (resb == 0) {
        seqInCore *sbb = B->getSequenceInCore(resultB[idxBb].i);

        mapDuplicates_Print(filea, sa, fileb, sbb);

        delete sbb;

        idxBb++;
        resb = md5_compare(resultA+idxA, resultB+idxBb);
      }

      //  And likewise for A
      //
      uint32 idxAa = idxA+1;
      int resa = md5_compare(resultA+idxAa, resultB+idxB);
      while (resa == 0) {
        seqInCore *saa = A->getSequenceInCore(resultA[idxAa].i);

        mapDuplicates_Print(filea, saa, fileb, sb);

        delete saa;

        idxAa++;
        resa = md5_compare(resultA+idxAa, resultB+idxB);
      }

      delete sa;
      delete sb;

      idxA++;
      idxB++;
    } else {
      if (res < 0)
        idxA++;
      else
        idxB++;
    }
  }

  delete A;
  delete B;
}

