#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "libbri.H"
#include "sim4polish.h"
#include "sim4polishBuilder.H"

#define MAX_POLISHES   10000000
#define MAX_ESTS        5000000  //  Should be set to exactly the number of ESTs
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
    repeatList[anInt] = 1;
  }

  fclose(F);

  return(repeatList);
}





//  Returns 1 if A is approximately the same polish as B
//          0 otherwise
//
int
comparePolish(sim4polish *A,
              sim4polish *B) {

  //  If not from the same EST/GEN pair, or mapped to different
  //  strands, they aren't compatible.
  //
  if ((A->estID != B->estID) ||
      (A->genID != B->genID) ||
      (A->matchOrientation != B->matchOrientation))
    return(0);

  if (s4p_IsSameRegion(A, B, REGION_TOLERANCE) &&
      s4p_IsSameExonModel(A, B, EXON_TOLERANCE))
    return(1);
  else
    return(0);
}




int
main(int argc, char **argv) {

  if (argc != 6) {
    fprintf(stderr, "usage: %s <polishes-filename> <minid> <mincov> <path-to-A> <path-to-B>\n", argv[0]);
    exit(1);
  }

  char       *Fname = argv[1];
  int         minI  = atoi(argv[2]);
  int         minC  = atoi(argv[3]);
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
  int               *Afound = new int [MAX_ESTS];
  int               *Bfound = new int [MAX_ESTS];

  for (int i=0; i<MAX_ESTS; i++)
    Afound[i] = Bfound[i] = 0;

  for (int i=0; i<A.length(); i++)
    if ((A[i]->percentIdentity  >= minI) ||
        (A[i]->querySeqIdentity >= minC))
      Afound[A[i]->estID]++;

  for (int i=0; i<B.length(); i++)
    if ((B[i]->percentIdentity  >= minI) ||
        (B[i]->querySeqIdentity >= minC))
      Bfound[B[i]->estID]++;


  int  numMatched  = 0;
  int  numRepeat   = 0;
  int  numProbable = 0;
  int  numExtra    = 0;
  int  numBelow    = 0;

  int  minB = 0;

  for (int i=0; i<A.length(); i++) {
    int  found = 0;

    //  If this polish is below our quality level, we don't care about it.
    //
    if ((A[i]->percentIdentity  < minI) ||
        (A[i]->querySeqIdentity < minC)) {
      numBelow++;
      continue;
    }

    int numFound   = 0;
    int numMissing = 0;
    int numExtra   = 0;
    int bestMatch  = 0;

    int j = minB;
    while ((j < B.length()) &&
           (A[i]->estID >= B[j]->estID)) {
      if (A[i]->estID > B[j]->estID)
        minB = j;

      if (comparePolish(A[i], B[j]))
        found++;

      //  Save the polish in B with the highest number of found exons.
      //
      int  f, m, e;

      s4p_compareExons_Overlap(A[i], B[j], &f, &m, &e);

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
        s4p_printPolish(stderr, A[i]);
      }
    } else {

      //  See if it's a repeat.
      //
      if (Brepeat[i]) {
        numRepeat++;

#if 0
        fprintf(stdout, "----------Found ZERO matches because it's a repeat.\n", i);
        s4p_printPolish(stdout, A[i]);
#endif
      } else if (numFound > 0) {
        numProbable++;

        fprintf(stdout, "----------Found PROBABLE match with %d -- numExonsFound=%d numExonsMissing=%d numExonsExtra=%d\n",
                bestMatch, numFound, numMissing, numExtra);
        s4p_printPolishNormalized(stdout, A[i]);
        s4p_printPolishNormalized(stdout, B[bestMatch]);
      } else {
        numExtra++;

        if ((A[i]->percentIdentity  >= minI) &&
            (A[i]->querySeqIdentity >= 80)) {
          fprintf(stdout, "----------Found ZERO matches HIGH-ID -- numExonsFound=%d numExonsMissing=%d numExonsExtra=%d\n",
                  numFound, numMissing, numExtra);
          s4p_printPolish(stdout, A[i]);
        } else {
          fprintf(stdout, "----------Found ZERO matches LOW-ID -- numExonsFound=%d numExonsMissing=%d numExonsExtra=%d\n",
                  numFound, numMissing, numExtra);
          s4p_printPolish(stdout, A[i]);
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


