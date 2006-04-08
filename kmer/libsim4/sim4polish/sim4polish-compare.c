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
  if (A->genLo + A->exons[0].genFrom < B->genLo + B->exons[0].genFrom) return(-1);
  if (A->genLo + A->exons[0].genFrom > B->genLo + B->exons[0].genFrom) return(1);

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
  if (A->genLo + A->exons[0].genFrom < B->genLo + B->exons[0].genFrom) return(-1);
  if (A->genLo + A->exons[0].genFrom > B->genLo + B->exons[0].genFrom) return(1);
  if (A->estID < B->estID) return(-1);
  if (A->estID > B->estID) return(1);

  return(0);
}



int
s4p_estDEFcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);
  int         e = 0;

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->estDefLine == 0L)  return(1);
  if (B->estDefLine == 0L)  return(-1);
  e = strcmp(A->estDefLine, B->estDefLine);
  if (e < 0) return(-1);
  if (e > 0) return(1);

  if (A->genDefLine == 0L)  return(1);
  if (B->genDefLine == 0L)  return(-1);
  e = strcmp(A->genDefLine, B->genDefLine);
  if (e < 0) return(-1);
  if (e > 0) return(1);

  if (A->genLo < B->genLo) return(-1);
  if (A->genLo > B->genLo) return(1);
  if (A->exons[0].genFrom < B->exons[0].genFrom) return(-1);
  if (A->exons[0].genFrom > B->exons[0].genFrom) return(1);

  return(0);
}



int
s4p_genDEFcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);
  int         e = 0;

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->genDefLine == 0L)  return(1);
  if (B->genDefLine == 0L)  return(-1);
  e = strcmp(A->genDefLine, B->genDefLine);
  if (e < 0) return(-1);
  if (e > 0) return(1);

  if (A->estDefLine == 0L)  return(1);
  if (B->estDefLine == 0L)  return(-1);
  e = strcmp(A->estDefLine, B->estDefLine);
  if (e < 0) return(-1);
  if (e > 0) return(1);

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
s4p_IsSameRegion(sim4polish *A, sim4polish *B, int tolerance) {
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



//  Returns true if the two polishes overlap genomic regions
//
int
s4p_IsRegionOverlap(sim4polish *A, sim4polish *B) {
  int Alo=0, Ahi=0;
  int Blo=0, Bhi=0;

  if (A->genID != B->genID)
    return(0);

  if (A->numExons > 0) {
    Alo = A->genLo + A->exons[0].genFrom;
    Ahi = A->genLo + A->exons[A->numExons-1].genTo;
  }

  if (B->numExons > 0) {
    Blo = B->genLo + B->exons[0].genFrom;
    Bhi = B->genLo + B->exons[B->numExons-1].genTo;
  }

  if (((Alo <= Blo) && (Blo <= Ahi)) ||
      ((Blo <= Alo) && (Alo <= Bhi)))
    return(1);
  return(0);
}



//  Returns true if the two polishes have the same number of exons,
//  and each exon is mapped to about the same genomic region.
//
int
s4p_IsSameExonModel(sim4polish *A, sim4polish *B, int tolerance) {
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
        (Dhi < -tolerance) || (Dhi > tolerance))
      return(0);
  }
  return(1);
}



void
s4p_compareExons_Overlap(sim4polish *A,
                         sim4polish *B,
                         double      overlapThreshold,
                         u32bit     *numSame,
                         u32bit     *numAMissed,
                         u32bit     *numBMissed) {
  u32bit    i, j;
  u32bit    al=0, ah=0, bl=0, bh=0;
  u32bit   *foundA = 0L;
  u32bit   *foundB = 0L;
  double    overlap = 0;

  if (numSame)     *numSame    = 0;
  if (numAMissed)  *numAMissed = 0;
  if (numBMissed)  *numBMissed = 0;

  errno = 0;

  foundA = (u32bit *)malloc(sizeof(u32bit) * (A->numExons + B->numExons));
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

      overlap = 0;

      //  Compute the percent overlapping as:
      //
      //     ----------
      //            ----------
      //            ^^^ = 3
      //     ^^^^^^^^^^^^^^^^^ = 17
      //
      //  overlap = 3/17
      //

      if ((al <= bl) && (bl <= ah)) {
        //  B starts somewhere in A
        //
        if (ah < bh) {
          //  B ends outside A
          //
          //  aaaaaaaaaaa
          //     bbbbbbbbbbbbb
          overlap = (double)(ah-bl) / (double)(bh-al);
        } else {
          //  B ends inside A
          //
          //  aaaaaaaaaaa
          //     bbbbb
          overlap = (double)(bh-bl) / (double)(ah-al);
        }
      }
      if ((bl <= al) && (al <= bh)) {
        //  B ends somewhere in A
        //
        if (bh < ah) {
          //  B starts outside A
          //
          //       aaaaaaaaaaa
          //  bbbbbbbbbbbbb
          overlap = (double)(bh-al) / (double)(ah-bl);
        } else {
          //  B starts inside A
          //
          //       aaaa
          //  bbbbbbbbbbbbb
          overlap = (double)(ah-al) / (double)(bh-bl);
        }
      }

      if (overlap >= overlapThreshold) {
        foundA[i]++;
        foundB[j]++;

        if (numSame)
          (*numSame)++;
      }
    }

  }

  for (i=0; i<A->numExons; i++) {
    //if (foundA[i] > 1)  fprintf(stderr, "WARNING: Found exon %d %d times in A!\n", i, foundA[i]);
    if (numAMissed && (foundA[i] == 0))
      (*numAMissed)++;
  }

  for (i=0; i<B->numExons; i++) {
    //if (foundB[i] > 1)  fprintf(stderr, "WARNING: Found exon %d %d times in B!\n", i, foundB[i]);
    if (numBMissed && (foundB[i] == 0))
      (*numBMissed)++;
  }

  free(foundA);
}





void
s4p_compareExons_Ends(sim4polish *A,
                      sim4polish *B,
                      s32bit     tolerance,
                      u32bit     *numSame,
                      u32bit     *numAMissed,
                      u32bit     *numBMissed) {
  u32bit    i, j;
  u32bit    Dlo=0, Dhi=0;
  u32bit   *foundA = 0L;
  u32bit   *foundB = 0L;

  if (numSame)     *numSame    = 0;
  if (numAMissed)  *numAMissed = 0;
  if (numBMissed)  *numBMissed = 0;

  foundA = (u32bit *)malloc(sizeof(u32bit) * (A->numExons + B->numExons));
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
    //if (foundA[i] > 1)  fprintf(stderr, "WARNING: Found exon %d %d times in A!\n", i, foundA[i]);
    if (numAMissed && (foundA[i] == 0))
      (*numAMissed)++;
  }

  for (i=0; i<B->numExons; i++) {
    //if (foundB[i] > 1)  fprintf(stderr, "WARNING: Found exon %d %d times in B!\n", i, foundB[i]);
    if (numBMissed && (foundB[i] == 0))
      (*numBMissed)++;
  }

  free(foundA);
}
