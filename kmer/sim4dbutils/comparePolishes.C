#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "bri++.H"
#include "sim4polish.h"
#include "sim4polishBuilder.H"

#define MAX_POLISHES   10000000
#define MAX_ESTS         100000  //  Should be set to exactly the number of ESTs
#define REGION_TOLERANCE     15
#define EXON_TOLERANCE       15

//  This will compare two runs of ESTmapper.  The est set and the
//  genomic must be the same for both runs.  It was used to evaluate
//  parameter choices during ESTmapper development.

int *
readRepeats(char *path) {
  int        *repeatList = new int [MAX_ESTS];
  char        repeatFile[2084];

  for (int i=0; i<MAX_ESTS; i++)
    repeatList[i] = 0;

  sprintf(repeatFile, "%s/2-filter/repeats", path);

  fprintf(stderr, "Reading repeats from %s\n", repeatFile);

  errno = 0;

  FILE *F = fopen(repeatFile, "r");
  if (errno) {
    fprintf(stderr, "Can't open repeats file '%s': %s\n", repeatFile, strerror(errno));
    fprintf(stderr, "Not labelling things as repeats.\n");
    return(repeatList);
  }

  int anInt;

  while (!feof(F)) {
    fscanf(F, " %d ", &anInt);
    if (!feof(F))
      repeatList[anInt] = 1;
  }

  fclose(F);

  return(repeatList);
}





//  Returns 1 if A is approximately the same polish as B
//          0 otherwise
//
int
isApproximatelySamePolish(sim4polish *A,
                          sim4polish *B) {

  //  If not from the same EST/GEN pair, or mapped to different
  //  strands, they aren't compatible.
  //
  if ((A->estID != B->estID) ||
      (A->genID != B->genID) ||
      (A->matchOrientation != B->matchOrientation))
    return(0);

  if (s4p_IsSameRegion(A, B, REGION_TOLERANCE) &&
      s4p_IsSameExonModel(A, B, EXON_TOLERANCE)) {
    return(1);
  } else {
    return(0);
  }
}



void
compareESTmapperDirectories(int argc, char **argv) {
  u32bit      minI  = atoi(argv[1]);
  u32bit      minC  = atoi(argv[2]);
  char       *Fname = argv[3];
  char       *Apath = argv[4];
  char       *Bpath = argv[5];

  //int        *Arepeat = readRepeats(Apath);
  int        *Brepeat = readRepeats(Bpath);

  char   filename[1025];

  sprintf(filename, "%s/%s", Apath, Fname);
  sim4polishList     A(filename);

  sprintf(filename, "%s/%s", Bpath, Fname);
  sim4polishList     B(filename);

  //  Construct lists of the ESTs found in either A or B.
  //
  u32bit            *Afound = new u32bit [MAX_ESTS];
  u32bit            *Bfound = new u32bit [MAX_ESTS];

  for (int i=0; i<MAX_ESTS; i++)
    Afound[i] = Bfound[i] = 0;

  for (u32bit i=0; i<A.length(); i++)
    if ((A[i]->percentIdentity  >= minI) ||
        (A[i]->querySeqIdentity >= minC))
      Afound[A[i]->estID]++;

  for (u32bit i=0; i<B.length(); i++)
    if ((B[i]->percentIdentity  >= minI) ||
        (B[i]->querySeqIdentity >= minC))
      Bfound[B[i]->estID]++;


  u32bit  numMatched  = 0;
  u32bit  numRepeat   = 0;
  u32bit  numProbable = 0;
  u32bit  numExtra    = 0;
  u32bit  numBelow    = 0;

  u32bit  minB = 0;

  for (u32bit i=0; i<A.length(); i++) {
    int  found = 0;

    //  If this polish is below our quality level, we don't care about it.
    //
    if ((A[i]->percentIdentity  < minI) ||
        (A[i]->querySeqIdentity < minC)) {
      numBelow++;
      continue;
    }

    u32bit numFound   = 0;
    u32bit numMissing = 0;
    u32bit numExtra   = 0;
    u32bit bestMatch  = 0;

    u32bit j = minB;
    while ((j < B.length()) &&
           (A[i]->estID >= B[j]->estID)) {
      if (A[i]->estID > B[j]->estID)
        minB = j;

      if (isApproximatelySamePolish(A[i], B[j]))
        found++;

      //  Save the polish in B with the highest number of found exons.
      //
      u32bit  f, m, e;

      s4p_compareExons_Overlap(A[i], B[j], 1e-10, &f, &m, &e);

      //  Check that we've not already found a best??

      if (numFound < f) {
        if (numFound != 0)
          fprintf(stderr, "WARNING:  Found multiple best matches for A=%d\n", i);

        numFound   = f;
        numMissing = m;
        numExtra   = e;
        bestMatch  = j;
      }

      j++;
    }

    if (found) {

      //  Yea!  We found it!  Make sure that we only found it once.

      numMatched++;

      if (found > 2) {
        fprintf(stderr, "Found %d matches for A=%d?\n", found, i);
        s4p_printPolish(stderr, A[i], S4P_PRINTPOLISH_FULL);
      }
    } else {

      //  See if it's a repeat.
      //
      if (Brepeat[i]) {
        numRepeat++;

#if 0
        fprintf(stdout, "----------Found ZERO matches because it's a repeat.\n", i);
        s4p_printPolish(stdout, A[i], S4P_PRINTPOLISH_FULL);
#endif
      } else if (numFound > 0) {
        numProbable++;

        fprintf(stdout, "----------Found PROBABLE match with %d -- numExonsFound=%d numExonsMissing=%d numExonsExtra=%d\n",
                bestMatch, numFound, numMissing, numExtra);
        s4p_printPolish(stdout, A[i], S4P_PRINTPOLISH_NORMALIZED);
        s4p_printPolish(stdout, B[bestMatch], S4P_PRINTPOLISH_NORMALIZED);
      } else {
        numExtra++;

        if ((A[i]->percentIdentity  >= minI) &&
            (A[i]->querySeqIdentity >= 80)) {
          fprintf(stdout, "----------Found ZERO matches HIGH-ID -- numExonsFound=%d numExonsMissing=%d numExonsExtra=%d\n",
                  numFound, numMissing, numExtra);
          s4p_printPolish(stdout, A[i], S4P_PRINTPOLISH_FULL);
        } else {
          fprintf(stdout, "----------Found ZERO matches LOW-ID -- numExonsFound=%d numExonsMissing=%d numExonsExtra=%d\n",
                  numFound, numMissing, numExtra);
          s4p_printPolish(stdout, A[i], S4P_PRINTPOLISH_FULL);
        }
      }
    }
  }

  //
  //  We should do the same statistics for only good?
  //

  fprintf(stdout, "--------------------------------------------------------------------------------\n");
  fprintf(stdout, "numBelow        = %d\n", numBelow);
  fprintf(stdout, "numMatched      = %d\n", numMatched);
  fprintf(stdout, "numRepeat       = %d\n", numRepeat);
  fprintf(stdout, "numProbable     = %d\n", numProbable);
  fprintf(stdout, "numExtra (in A) = %d\n", numExtra);

  //  Count the number of times an EST was found by one but not the
  //  other.  Don't care if both found it, or both missed it.
  //
  int   onlyAfound = 0;
  int   onlyBfound = 0;

  for (int i=0; i<MAX_ESTS; i++) {
    if (Afound[i] && !Bfound[i])
      onlyAfound++;
    if (Bfound[i] && !Afound[i])
      onlyBfound++;
  }

  fprintf(stdout, "only A found    = %d (cdna)\n", onlyAfound);
  fprintf(stdout, "only B found    = %d (cdna)\n", onlyBfound);
}



void
comparePolishFiles(int argc, char **argv) {
  u32bit      minI  = atoi(argv[1]);
  u32bit      minC  = atoi(argv[2]);
  char       *Apath = argv[3];
  char       *Bpath = argv[4];

  fprintf(stderr, "reading A from %s\n", Apath);
  sim4polishList     A(Apath);
  A.sortBycDNAIID();

  fprintf(stderr, "reading B from %s\n", Bpath);
  sim4polishList     B(Bpath);
  B.sortBycDNAIID();

  u32bit  thisA = 0;
  u32bit  thisB = 0;
  u32bit  lastB = 0;

  u32bit  Bmissed  = 0;
  u32bit  Boverlap = 0;
  u32bit  Abetter  = 0;
  u32bit  Bbetter  = 0;
  u32bit  Nbetter  = 0;
  u32bit  Equal    = 0;

  while (thisA < A.length()) {
    u32bit  largestF = 0;
    u32bit  largestA = 0;
    u32bit  largestB = 0;
    u32bit  bestB    = 0;

#if 1
    //  Remember the first B that has the correct estID
    //
    while ((A[thisA]->estID > B[lastB]->estID) && (lastB < B.length()))
      lastB++;
#endif

    thisB = lastB;

    bool    overlapped = false;

    //  Scan forward in B comparing matches with the same IID.
    //
    while ((A[thisA]->estID == B[thisB]->estID) && (thisB < B.length())) {
      //while (thisB < B.length()) {
      if ((A[thisA]->estID            == B[thisB]->estID) &&
          (A[thisA]->genID            == B[thisB]->genID) &&
          (A[thisA]->matchOrientation == B[thisB]->matchOrientation)) {

        overlapped = s4p_IsRegionOverlap(A[thisA], B[thisB]);

        u32bit  f, a, b;

        //s4p_compareExons_Ends(A[thisA], B[thisB], 15, &f, &a, &b);
        s4p_compareExons_Overlap(A[thisA], B[thisB], 0.75, &f, &a, &b);

        //  If we haven't found a good match yet, and we
        //  are overlapping, save thisB
        //
        if ((largestF == 0) && (bestB == 0) && (overlapped)) {
          bestB = thisB;
        }

        //  We want to maximize F while minimizing A+B
        //
        //  This hopefully will get us around the problem of having
        //  duplicate matches that are subsets:
        //    match 1 has exons 1 2 3 4 5 6
        //    match 2 has exons     3 4 5
        //
        //  When finding the correct pair for match 2, we need to
        //  ignore match 1, because it has 3 missing exons, and
        //  pick match 2 with no missing exons.
        //
        if ((largestF < f) ||
            ((largestF <= f) && (largestA + largestB >= a + b))) {
          largestF = f;
          largestA = a;
          largestB = b;
          bestB = thisB;
        }
      }

      thisB++;
    }

    if (largestF == 0) {
      if (overlapped) {
        //  Thing in A overlapped something in B
        fprintf(stdout, "A not found in B, but overlapped\n");
        s4p_printPolish(stdout, A[thisA], S4P_PRINTPOLISH_FULL | S4P_PRINTPOLISH_NORMALIZED);
        s4p_printPolish(stdout, B[bestB], S4P_PRINTPOLISH_FULL | S4P_PRINTPOLISH_NORMALIZED);
        Boverlap++;
      } else {
        //  Thing in A was not found in B
        fprintf(stdout, "A not found in B\n");
        s4p_printPolish(stdout, A[thisA], S4P_PRINTPOLISH_FULL | S4P_PRINTPOLISH_NORMALIZED);
        Bmissed++;
      }
    } else if ((largestF == A[thisA]->numExons) && (largestF == B[bestB]->numExons)) {
      //  We matched all exons
      Equal++;
    } else if ((largestA == 0) && (largestB == 0)) {
      //  Shouldn't happen; didn't match all exons, and didn't miss any.
      fprintf(stderr, "\nBoth A and B are zero and f != num exons  f=%d  a=%d  b=%d??\n", largestF, A[thisA]->numExons, B[bestB]->numExons);
    } else if ((largestA > 0) && (largestB > 0)) {
      //  They both have exons not matched in the other one.
      Nbetter++;
    } else if (largestA > 0) {
      //  A has extra exons, so A is better.
      Abetter++;
      //fprintf(stdout, "Abetter f=%f a=%d b=%d\n", largestF, largestA, largestB);
      //s4p_printPolish(stdout, A[thisA], S4P_PRINTPOLISH_FULL | S4P_PRINTPOLISH_NORMALIZED);
      //s4p_printPolish(stdout, B[bestB], S4P_PRINTPOLISH_FULL | S4P_PRINTPOLISH_NORMALIZED);
    } else if (largestB > 0) {
      //  B has extra exons, so B is better.
      Bbetter++;
      //fprintf(stdout, "Bbetter f=%f a=%d b=%d\n", largestF, largestA, largestB);
      //s4p_printPolish(stdout, A[thisA], S4P_PRINTPOLISH_FULL | S4P_PRINTPOLISH_NORMALIZED);
      //s4p_printPolish(stdout, B[bestB], S4P_PRINTPOLISH_FULL | S4P_PRINTPOLISH_NORMALIZED);
    } else {
      fprintf(stderr, "\nUnmatched case at a=%d b=%d!\n", thisA, bestB);
    }

#if 0
    s4p_printPolish(stderr, A[thisA], S4P_PRINTPOLISH_FULL | S4P_PRINTPOLISH_NORMALIZED);
    s4p_printPolish(stderr, A[bestB], S4P_PRINTPOLISH_FULL | S4P_PRINTPOLISH_NORMALIZED);
#endif

    thisA++;

#if 1
    fprintf(stderr, "A=%6d  B=%6d  Equal %6d  Bmissed %6d  Boverlap %6d  Abetter %6d  Bbetter %6d  Nbetter %6d\r",
            thisA, lastB, Equal, Bmissed, Boverlap, Abetter, Bbetter, Nbetter);
    fflush(stderr);
#endif
  }

  fprintf(stderr, "A=%6d  B=%6d  Equal %6d  Bmissed %6d  Boverlap %6d  Abetter %6d  Bbetter %6d  Nbetter %6d\n",
          thisA, lastB, Equal, Bmissed, Boverlap, Abetter, Bbetter, Nbetter);
}








int
main(int argc, char **argv) {

  if ((argc != 6) && (argc != 5)) {
    fprintf(stderr, "usage: %s <minid> <mincov> <polishes-filename> <path-to-A> <path-to-B>\n", argv[0]);
    fprintf(stderr, "       %s <minid> <mincov> <polishes-file-1> <polishes-file-2>\n", argv[0]);
    exit(1);
  }

  if (argc == 6)
    compareESTmapperDirectories(argc, argv);

  if (argc == 5)
    comparePolishFiles(argc, argv);
}


