#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "libbri.H"
#include "sim4reader.h"

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



struct sortedPolishSet {
  sim4polish **polishes;
  int          num;
};



sortedPolishSet *
readPolishes(char *path, char *name) {
  sortedPolishSet  *p = new sortedPolishSet;
  char              polishFile[2084];

  sprintf(polishFile, "%s/%s", path, name);

  fprintf(stderr, "Reading polished from %s\n", polishFile);

  errno = 0;

  FILE *F = fopen(polishFile, "r");
  if (errno) {
    fprintf(stderr, "Can't open polishes '%s': %s\n", polishFile, strerror(errno));
    exit(1);
  }

  p->polishes = new sim4polish * [MAX_POLISHES];
  p->num      = 0;

  while (!feof(F)) {
    p->polishes[p->num] = readPolish(F);
    p->num++;

    if (p->num > MAX_POLISHES) {
      fprintf(stderr, "ERROR:  MAX_POLISHES too small!\n");
      exit(1);
    }
  }
  fclose(F);

  p->num--;

  fprintf(stderr, "Read %d polishes -- sorting.\n", p->num);

  qsort(p->polishes, p->num, sizeof(sim4polish *), estIDcompare);

  return(p);
}






//  Returns 1 if A is approximately the same polish as B
//          0 otherwise
//
int
comparePolish(sim4polish   *A,
              sim4polish   *B) {

  //  If not from the same EST/GEN pair, or mapped to different
  //  strands, they aren't compatible.
  //
  if ((A->estID != B->estID) ||
      (A->genID != B->genID) ||
      (A->matchOrientation != B->matchOrientation))
    return(0);

  //
  //  We don't really care about the search region, but we do want
  //  to ensure that the polishes are on (roughly) the same region.
  //
  int Alo, Ahi;
  int Blo, Bhi;
  int Dlo, Dhi;

  Alo = A->genLo + A->exons[0].genFrom;
  Ahi = A->genLo + A->exons[A->numExons-1].genTo;

  Blo = B->genLo + B->exons[0].genFrom;
  Bhi = B->genLo + B->exons[B->numExons-1].genTo;

  Dlo = Blo - Alo;
  Dhi = Bhi - Ahi;

  if ((Dlo < -REGION_TOLERANCE) || (Dlo > REGION_TOLERANCE) ||
      (Dhi < -REGION_TOLERANCE) || (Dhi > REGION_TOLERANCE)) {
    return(0);
  }

  //
  //  We are the same EST/GEN pair, and are on about the same region.
  //  Check the exons.
  //

  if (A->numExons != B->numExons)
    return(0);

  for (int i=0; i<A->numExons; i++) {
    int Alo, Ahi;
    int Blo, Bhi;
    int Dlo, Dhi;

    Alo = A->genLo + A->exons[i].genFrom;
    Ahi = A->genLo + A->exons[i].genTo;

    Blo = B->genLo + B->exons[i].genFrom;
    Bhi = B->genLo + B->exons[i].genTo;

    Dlo = Blo - Alo;
    Dhi = Bhi - Ahi;

    if ((Dlo < -EXON_TOLERANCE) || (Dlo > EXON_TOLERANCE) ||
        (Dhi < -EXON_TOLERANCE) || (Dhi > EXON_TOLERANCE)) {
      return(0);
    }
  }

  return(1);
}




//  Returns the number of exons that are the same, missing or extra
//  (relative to A).
//
void
compareExons(sim4polish   *A,
             sim4polish   *B,
             int          &numSame,
             int          &numMissing,
             int          &numExtra) {

  numSame    = 0;
  numMissing = 0;
  numExtra   = 0;

  //  If not from the same EST/GEN pair, or mapped to different
  //  strands, call the sort function given.
  //
  if ((A->estID != B->estID) ||
      (A->genID != B->genID) ||
      (A->matchOrientation != B->matchOrientation))
    return;

  int  *foundA = new int [A->numExons];
  int  *foundB = new int [B->numExons];

  for (int i=0; i<A->numExons; i++)
    foundA[i] = 0;

  for (int i=0; i<B->numExons; i++)
    foundB[i] = 0;

  for (int i=0; i<A->numExons; i++) {
    for (int j=0; j<B->numExons; j++) {

#if 0
      //  If they have similar end points, declare a match
      //
      int Dlo = (B->genLo + B->exons[j].genFrom) - (A->genLo + A->exons[i].genFrom);
      int Dhi = (B->genLo + B->exons[j].genTo)   - (A->genLo + A->exons[i].genTo);

      if ((Dlo > -EXON_TOLERANCE) && (Dlo < EXON_TOLERANCE) &&
          (Dhi > -EXON_TOLERANCE) && (Dhi < EXON_TOLERANCE)) {
        foundA[i]++;
        foundB[j]++;
        numSame++;
      }
#else
      //  If they overlap, declare a match
      //
      int al = A->genLo + A->exons[i].genFrom;
      int ah = A->genLo + A->exons[i].genTo;
      int bl = B->genLo + B->exons[j].genFrom;
      int bh = B->genLo + B->exons[j].genTo;

      if (((al <= bl) && (bl <= ah) && (ah <= bh)) ||
          ((bl <= al) && (al <= bh) && (bh <= ah)) ||
          ((al <= bl) && (bh <= ah)) ||
          ((bl <= al) && (ah <= bh))) {
        foundA[i]++;
        foundB[j]++;
        numSame++;
      }
#endif

    }
  }

  for (int i=0; i<A->numExons; i++) {
    if (foundA[i] == 0)
      numExtra++;
    if (foundA[i] > 1)
      fprintf(stderr, "WARNING: Found exon %d %d times in A!\n", i, foundA[i]);
  }

  for (int i=0; i<B->numExons; i++) {
    if (foundB[i] == 0)
      numMissing++;
    if (foundB[i] > 1)
      fprintf(stderr, "WARNING: Found exon %d %d times in B!\n", i, foundB[i]);
  }
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

  sortedPolishSet   *A = readPolishes(Apath, Fname);
  sortedPolishSet   *B = readPolishes(Bpath, Fname);

  //  Construct lists of the ESTs found in either A or B.
  //
  int               *Afound = new int [MAX_ESTS];
  int               *Bfound = new int [MAX_ESTS];

  for (int i=0; i<MAX_ESTS; i++)
    Afound[i] = Bfound[i] = 0;

  for (int i=0; i<A->num; i++)
    if ((A->polishes[i]->percentIdentity  >= minI) ||
        (A->polishes[i]->querySeqIdentity >= minC))
      Afound[A->polishes[i]->estID]++;

  for (int i=0; i<B->num; i++)
    if ((B->polishes[i]->percentIdentity  >= minI) ||
        (B->polishes[i]->querySeqIdentity >= minC))
      Bfound[B->polishes[i]->estID]++;


  int  numMatched  = 0;
  int  numRepeat   = 0;
  int  numProbable = 0;
  int  numExtra    = 0;
  int  numBelow    = 0;

  int  minB = 0;

  for (int i=0; i<A->num; i++) {
    int  found = 0;

    //  If this polish is below our quality level, we don't care about it.
    //
    if ((A->polishes[i]->percentIdentity  < minI) ||
        (A->polishes[i]->querySeqIdentity < minC)) {
      numBelow++;
      continue;
    }

    int numFound   = 0;
    int numMissing = 0;
    int numExtra   = 0;
    int bestMatch  = 0;

    int j = minB;
    while ((j < B->num) &&
           (A->polishes[i]->estID >= B->polishes[j]->estID)) {
      if (A->polishes[i]->estID > B->polishes[j]->estID)
        minB = j;

      if (comparePolish(A->polishes[i], B->polishes[j]))
        found++;

      //  Save the polish in B with the highest number of found exons.
      //
      int  f, m, e;

      compareExons(A->polishes[i], B->polishes[j], f, m, e);

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
        printPolish(stderr, A->polishes[i]);
      }
    } else {

      //  See if it's a repeat.
      //
      if (Brepeat[i]) {
        numRepeat++;

#if 0
        fprintf(stdout, "----------Found ZERO matches because it's a repeat.\n", i);
        printPolish(stdout, A->polishes[i]);
#endif
      } else if (numFound > 0) {
        numProbable++;

        fprintf(stdout, "----------Found PROBABLE match with %d -- numExonsFound=%d numExonsMissing=%d numExonsExtra=%d\n",
                bestMatch, numFound, numMissing, numExtra);
        printPolishNormalized(stdout, A->polishes[i]);
        printPolishNormalized(stdout, B->polishes[bestMatch]);
      } else {
        numExtra++;

        if ((A->polishes[i]->percentIdentity  >= minI) &&
            (A->polishes[i]->querySeqIdentity >= 80)) {
          fprintf(stdout, "----------Found ZERO matches HIGH-ID -- numExonsFound=%d numExonsMissing=%d numExonsExtra=%d\n",
                  numFound, numMissing, numExtra);
          printPolish(stdout, A->polishes[i]);
        } else {
          fprintf(stdout, "----------Found ZERO matches LOW-ID -- numExonsFound=%d numExonsMissing=%d numExonsExtra=%d\n",
                  numFound, numMissing, numExtra);
          printPolish(stdout, A->polishes[i]);
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


