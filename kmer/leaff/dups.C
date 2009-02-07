#include "bio++.H"
#include "seqCache.H"


md5_s *
computeMD5ForEachSequence(seqCache *F) {
  u32bit   numSeqs = F->getNumberOfSequences();
  md5_s   *result  = new md5_s [numSeqs];

  for (u32bit idx=0; idx < numSeqs; idx++) {
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
    fprintf(stdout, u32bitFMT" <-> "u32bitFMT"\n", sa->getIID(), sb->getIID());
  else
    fprintf(stderr, "COLLISION DETECTED BETWEEN %s:"u32bitFMT" AND %s:"u32bitFMT"!\nPLEASE REPORT THIS TO bri@walenz.org!\n",
            filea, sa->getIID(), fileb, sb->getIID());
}



void
findDuplicates(char *filename) {
  seqInCore  *s1 = 0L;
  seqInCore  *s2 = 0L;
  seqCache   *A = new seqCache(filename);

  u32bit numSeqs = A->getNumberOfSequences();

  fprintf(stderr, "Computing MD5's for each sequence in '%s'.\n", filename);
  md5_s *result = computeMD5ForEachSequence(A);

  fprintf(stderr, "Sorting MD5's.\n");
  qsort(result, numSeqs, sizeof(md5_s), md5_compare);

  fprintf(stderr, "Verifying identity, and output\n");
  for (u32bit idx=1; idx<numSeqs; idx++) {
    if (md5_compare(result+idx-1, result+idx) == 0) {
      if (result[idx-1].i == result[idx].i) {
        fprintf(stderr, "Internal error: found two copies of the same sequence iid ("u32bitFMT")!\n", result[idx].i);
        exit(1);
      }

      s1 = A->getSequenceInCore(result[idx-1].i);
      s2 = A->getSequenceInCore(result[idx].i);

      if (strcmp(s1->sequence(), s2->sequence()) == 0) {
        fprintf(stdout, u32bitFMT":%s\n"u32bitFMT":%s\n\n",
                result[idx-1].i, s1->header(),
                result[idx  ].i, s2->header());
      } else {
        fprintf(stderr, "COLLISION DETECTED BETWEEN IID "u32bitFMT" AND "u32bitFMT"!\nPLEASE REPORT THIS TO bri@walenz.org!\n",
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

  u32bit  numSeqsA = A->getNumberOfSequences();
  u32bit  numSeqsB = B->getNumberOfSequences();
  u32bit  idxA = 0;
  u32bit  idxB = 0;

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
      u32bit idxBb = idxB+1;
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
      u32bit idxAa = idxA+1;
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

