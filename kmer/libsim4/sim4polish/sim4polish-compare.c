#include "sim4polish.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

//
//  Routines for comparing sim4polish structures.
//
//  Many of these routines assume that the iid's are consistent for
//  the pair of polishes.  In particular, that they are mapped to the
//  same set of genomic sequences.
//


int
s4p_estIDcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->estID < B->estID) return(-1);
  if (A->estID > B->estID) return(1);
  if (A->genID < B->genID) return(-1);
  if (A->genID > B->genID) return(1);
  if (A->genLo < B->genLo) return(-1);
  if (A->genLo > B->genLo) return(1);
  if (A->exons[0].genFrom < B->exons[0].genFrom) return(-1);
  if (A->exons[0].genFrom > B->exons[0].genFrom) return(1);

  return(0);
}



int
s4p_genIDcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->genID < B->genID) return(-1);
  if (A->genID > B->genID) return(1);
  if (A->estID < B->estID) return(-1);
  if (A->estID > B->estID) return(1);
  if (A->genLo < B->genLo) return(-1);
  if (A->genLo > B->genLo) return(1);
  if (A->exons[0].genFrom < B->exons[0].genFrom) return(-1);
  if (A->exons[0].genFrom > B->exons[0].genFrom) return(1);

  return(0);
}





//  Return false if not from the same EST/GEN pair, or mapped to
//  different strands, true otherwise.
//
int
s4p_compatable(sim4polish *A, sim4polish *B) {
  if ((A->estID != B->estID) ||
      (A->genID != B->genID) ||
      (A->matchOrientation != B->matchOrientation))
    return(0);
  else
    return(1);
}




//  Returns true if the two polishes are on about the same genomic
//  region
//
int
s4p_IsSameRegion(sim4polish *A, sim4polish *B, u32bit tolerance) {
  int Alo=0, Ahi=0;
  int Blo=0, Bhi=0;
  int Dlo=0, Dhi=0;

  if (A->numExons > 0) {
    Alo = A->genLo + A->exons[0].genFrom;
    Ahi = A->genLo + A->exons[A->numExons-1].genTo;
  }

  if (B->numExons > 0) {
    Blo = B->genLo + B->exons[0].genFrom;
    Bhi = B->genLo + B->exons[B->numExons-1].genTo;
  }

  Dlo = Blo - Alo;
  Dhi = Bhi - Ahi;

  if ((Dlo < -tolerance) || (Dlo > tolerance) ||
      (Dhi < -tolerance) || (Dhi > tolerance))
    return(0);

  return(1);
}



//  Returns true if the two polishes have the same number of exons,
//  and each exon is mapped to about the same genomic region.
//
int
s4p_IsSameExonModel(sim4polish *A, sim4polish *B, u32bit tolerance) {
  int i;
  int Alo=0, Ahi=0;
  int Blo=0, Bhi=0;
  int Dlo=0, Dhi=0;

  if (A->numExons != B->numExons)
    return(0);

  for (i=0; i<A->numExons; i++) {
    Alo = A->genLo + A->exons[i].genFrom;
    Ahi = A->genLo + A->exons[i].genTo;

    Blo = B->genLo + B->exons[i].genFrom;
    Bhi = B->genLo + B->exons[i].genTo;

    Dlo = Blo - Alo;
    Dhi = Bhi - Ahi;

    if ((Dlo < -tolerance) || (Dlo > tolerance) ||
        (Dhi < -tolerance) || (Dhi > tolerance)) {
      return(0);
    }
  }
  return(1);
}




void
s4p_compareExons_Overlap(sim4polish *A,
                         sim4polish *B,
                         int        *numSame,
                         int        *numMissing,
                         int        *numExtra) {
  int       i, j;
  //int       Dlo=0, Dhi=0;
  int       al=0, ah=0, bl=0, bh=0;
  int      *foundA = 0L;
  int      *foundB = 0L;

  if (numSame)     *numSame    = 0;
  if (numMissing)  *numMissing = 0;
  if (numExtra)    *numExtra   = 0;

  errno = 0;

  foundA = (int *)malloc(sizeof(int) * (A->numExons + B->numExons));
  foundB = foundA + A->numExons;

  if (errno) {
    fprintf(stderr, "s4p_compareExons()-- Can't allocate "u32bitFMT" + "u32bitFMT" words for counting exons.\n%s\n", A->numExons, B->numExons, strerror(errno));
    exit(1);
  }

  for (i=0; i<A->numExons; i++)
    foundA[i] = 0;

  for (i=0; i<B->numExons; i++)
    foundB[i] = 0;

  //  If they overlap, declare a match
  //
  for (i=0; i<A->numExons; i++) {
    for (j=0; j<B->numExons; j++) {
      al = A->genLo + A->exons[i].genFrom;
      ah = A->genLo + A->exons[i].genTo;
      bl = B->genLo + B->exons[j].genFrom;
      bh = B->genLo + B->exons[j].genTo;

      if (((al <= bl) && (bl <= ah) && (ah <= bh)) ||
          ((bl <= al) && (al <= bh) && (bh <= ah)) ||
          ((al <= bl) && (bh <= ah)) ||
          ((bl <= al) && (ah <= bh))) {
        foundA[i]++;
        foundB[j]++;

        if (numSame)
          (*numSame)++;
      }
    }
  }

  for (i=0; i<A->numExons; i++) {
    if (numExtra && (foundA[i] == 0))
      (*numExtra)++;
#if 0
    if (foundA[i] > 1)
      fprintf(stderr, "WARNING: Found exon %d %d times in A!\n", i, foundA[i]);
#endif
  }

  for (i=0; i<B->numExons; i++) {
    if (numMissing && (foundB[i] == 0))
      (*numMissing)++;
#if 0
    if (foundB[i] > 1)
      fprintf(stderr, "WARNING: Found exon %d %d times in B!\n", i, foundB[i]);
#endif
  }

  free(foundA);
}





void
s4p_compareExons_Ends(sim4polish *A,
                      sim4polish *B,
                      int         tolerance,
                      int        *numSame,
                      int        *numMissing,
                      int        *numExtra) {
  int       i, j;
  int       Dlo=0, Dhi=0;
  //int       al=0, ah=0, bl=0, bh=0;
  int      *foundA = 0L;
  int      *foundB = 0L;

  if (numSame)     *numSame    = 0;
  if (numMissing)  *numMissing = 0;
  if (numExtra)    *numExtra   = 0;

  foundA = (int *)malloc(sizeof(int) * (A->numExons + B->numExons));
  foundB = foundA + A->numExons;

  if (errno) {
    fprintf(stderr, "s4p_compareExons()-- Can't allocate "u32bitFMT" + "u32bitFMT" words for counting exons.\n%s\n", A->numExons, B->numExons, strerror(errno));
    exit(1);
  }

  for (i=0; i<A->numExons; i++)
    foundA[i] = 0;

  for (i=0; i<B->numExons; i++)
    foundB[i] = 0;

  //  If they have similar end points, declare a match
  //
  for (i=0; i<A->numExons; i++) {
    for (j=0; j<B->numExons; j++) {
      Dlo = (int)(B->genLo + B->exons[j].genFrom) - (int)(A->genLo + A->exons[i].genFrom);
      Dhi = (int)(B->genLo + B->exons[j].genTo)   - (int)(A->genLo + A->exons[i].genTo);

      if ((Dlo > -tolerance) && (Dlo < tolerance) &&
          (Dhi > -tolerance) && (Dhi < tolerance)) {
        foundA[i]++;
        foundB[j]++;

        if (numSame)
          (*numSame)++;
      }
    }
  }

  for (i=0; i<A->numExons; i++) {
    if (numExtra && (foundA[i] == 0))
      (*numExtra)++;
#if 0
    if (foundA[i] > 1)
      fprintf(stderr, "WARNING: Found exon %d %d times in A!\n", i, foundA[i]);
#endif
  }

  for (i=0; i<B->numExons; i++) {
    if (numMissing && (foundB[i] == 0))
      (*numMissing)++;
#if 0
    if (foundB[i] > 1)
      fprintf(stderr, "WARNING: Found exon %d %d times in B!\n", i, foundB[i]);
#endif
  }

  free(foundA);
}




