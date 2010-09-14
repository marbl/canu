#include "sim4polish.H"
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

  if (A->_estID < B->_estID) return(-1);
  if (A->_estID > B->_estID) return(1);
  if (A->_genID < B->_genID) return(-1);
  if (A->_genID > B->_genID) return(1);
  if (A->_exons[0]._genFrom < B->_exons[0]._genFrom) return(-1);
  if (A->_exons[0]._genFrom > B->_exons[0]._genFrom) return(1);

  return(0);
}



int
s4p_genIDcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->_genID < B->_genID) return(-1);
  if (A->_genID > B->_genID) return(1);
  if (A->_exons[0]._genFrom < B->_exons[0]._genFrom) return(-1);
  if (A->_exons[0]._genFrom > B->_exons[0]._genFrom) return(1);
  if (A->_estID < B->_estID) return(-1);
  if (A->_estID > B->_estID) return(1);

  return(0);
}



int
s4p_estDEFcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);
  int         e = 0;

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->_estDefLine == 0L)  return(1);
  if (B->_estDefLine == 0L)  return(-1);
  e = strcmp(A->_estDefLine, B->_estDefLine);
  if (e < 0) return(-1);
  if (e > 0) return(1);

  if (A->_genDefLine == 0L)  return(1);
  if (B->_genDefLine == 0L)  return(-1);
  e = strcmp(A->_genDefLine, B->_genDefLine);
  if (e < 0) return(-1);
  if (e > 0) return(1);

  if (A->_exons[0]._genFrom < B->_exons[0]._genFrom) return(-1);
  if (A->_exons[0]._genFrom > B->_exons[0]._genFrom) return(1);

  return(0);
}



int
s4p_genDEFcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);
  int         e = 0;

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->_genDefLine == 0L)  return(1);
  if (B->_genDefLine == 0L)  return(-1);
  e = strcmp(A->_genDefLine, B->_genDefLine);
  if (e < 0) return(-1);
  if (e > 0) return(1);

  if (A->_estDefLine == 0L)  return(1);
  if (B->_estDefLine == 0L)  return(-1);
  e = strcmp(A->_estDefLine, B->_estDefLine);
  if (e < 0) return(-1);
  if (e > 0) return(1);

  if (A->_exons[0]._genFrom < B->_exons[0]._genFrom) return(-1);
  if (A->_exons[0]._genFrom > B->_exons[0]._genFrom) return(1);

  return(0);
}







//  Return false if not from the same EST/GEN pair, or mapped to
//  different strands, true otherwise.
//
bool
s4p_compatable(sim4polish *A, sim4polish *B) {
  if ((A->_estID != B->_estID) ||
      (A->_genID != B->_genID) ||
      (A->_matchOrientation != B->_matchOrientation))
    return(false);
  else
    return(true);
}




//  Returns true if the two polishes are on about the same genomic
//  region
//
bool
s4p_IsSameRegion(sim4polish *A, sim4polish *B, int tolerance) {
  s32bit Alo=0, Ahi=0;
  s32bit Blo=0, Bhi=0;
  s32bit Dlo=0, Dhi=0;

  if (A->_numExons > 0) {
    Alo = (s32bit)A->_exons[0]._genFrom;
    Ahi = (s32bit)A->_exons[A->_numExons-1]._genTo;
  }

  if (B->_numExons > 0) {
    Blo = (s32bit)B->_exons[0]._genFrom;
    Bhi = (s32bit)B->_exons[B->_numExons-1]._genTo;
  }

  Dlo = Blo - Alo;
  Dhi = Bhi - Ahi;

  if ((Dlo < -tolerance) || (Dlo > tolerance) ||
      (Dhi < -tolerance) || (Dhi > tolerance))
    return(false);
  else
    return(true);
}



//  Returns true if the two polishes overlap genomic regions
//
bool
s4p_IsRegionOverlap(sim4polish *A, sim4polish *B) {
  s32bit Alo=0, Ahi=0;
  s32bit Blo=0, Bhi=0;

  if (A->_genID != B->_genID)
    return(false);

  if (A->_numExons > 0) {
    Alo = (s32bit)A->_exons[0]._genFrom;
    Ahi = (s32bit)A->_exons[A->_numExons-1]._genTo;
  }

  if (B->_numExons > 0) {
    Blo = (s32bit)B->_exons[0]._genFrom;
    Bhi = (s32bit)B->_exons[B->_numExons-1]._genTo;
  }

  if (((Alo <= Blo) && (Blo <= Ahi)) ||
      ((Blo <= Alo) && (Alo <= Bhi)))
    return(true);
  else
    return(false);
}



//  Returns true if the two polishes have the same number of exons,
//  and each exon is mapped to about the same genomic region.
//
bool
s4p_IsSameExonModel(sim4polish *A, sim4polish *B, int tolerance) {
  s32bit Alo=0, Ahi=0;
  s32bit Blo=0, Bhi=0;
  s32bit Dlo=0, Dhi=0;

  if (A->_numExons != B->_numExons)
    return(0);

  for (u32bit i=0; i<A->_numExons; i++) {
    Alo = (s32bit)A->_exons[i]._genFrom;
    Ahi = (s32bit)A->_exons[i]._genTo;

    Blo = (s32bit)B->_exons[i]._genFrom;
    Bhi = (s32bit)B->_exons[i]._genTo;

    Dlo = Blo - Alo;
    Dhi = Bhi - Ahi;

    if ((Dlo < -tolerance) || (Dlo > tolerance) ||
        (Dhi < -tolerance) || (Dhi > tolerance))
      return(false);
  }

  return(true);
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

  foundA = new u32bit [A->_numExons + B->_numExons];
  foundB = foundA + A->_numExons;

  if (errno) {
    fprintf(stderr, "s4p_compareExons()-- Can't allocate "u32bitFMT" + "u32bitFMT" words for counting exons.\n%s\n", A->_numExons, B->_numExons, strerror(errno));
    exit(1);
  }

  for (i=0; i<A->_numExons; i++)
    foundA[i] = 0;

  for (i=0; i<B->_numExons; i++)
    foundB[i] = 0;

  //  If they overlap, declare a match
  //
  for (i=0; i<A->_numExons; i++) {
    for (j=0; j<B->_numExons; j++) {
      al = A->_exons[i]._genFrom;
      ah = A->_exons[i]._genTo;
      bl = B->_exons[j]._genFrom;
      bh = B->_exons[j]._genTo;

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

  for (i=0; i<A->_numExons; i++) {
    //if (foundA[i] > 1)  fprintf(stderr, "WARNING: Found exon %d %d times in A!\n", i, foundA[i]);
    if (numAMissed && (foundA[i] == 0))
      (*numAMissed)++;
  }

  for (i=0; i<B->_numExons; i++) {
    //if (foundB[i] > 1)  fprintf(stderr, "WARNING: Found exon %d %d times in B!\n", i, foundB[i]);
    if (numBMissed && (foundB[i] == 0))
      (*numBMissed)++;
  }

  delete [] foundA;
}





void
s4p_compareExons_Ends(sim4polish *A,
                      sim4polish *B,
                      s32bit     tolerance,
                      u32bit     *numSame,
                      u32bit     *numAMissed,
                      u32bit     *numBMissed) {
  u32bit    i, j;
  s32bit    Dlo=0, Dhi=0;
  u32bit   *foundA = 0L;
  u32bit   *foundB = 0L;

  if (numSame)     *numSame    = 0;
  if (numAMissed)  *numAMissed = 0;
  if (numBMissed)  *numBMissed = 0;

  foundA = new u32bit [A->_numExons + B->_numExons];
  foundB = foundA + A->_numExons;

  if (errno) {
    fprintf(stderr, "s4p_compareExons()-- Can't allocate "u32bitFMT" + "u32bitFMT" words for counting exons.\n%s\n", A->_numExons, B->_numExons, strerror(errno));
    exit(1);
  }

  for (i=0; i<A->_numExons; i++)
    foundA[i] = 0;

  for (i=0; i<B->_numExons; i++)
    foundB[i] = 0;

  //  If they have similar end points, declare a match
  //
  for (i=0; i<A->_numExons; i++) {
    for (j=0; j<B->_numExons; j++) {
      Dlo = (s32bit)(B->_exons[j]._genFrom) - (s32bit)(A->_exons[i]._genFrom);
      Dhi = (s32bit)(B->_exons[j]._genTo)   - (s32bit)(A->_exons[i]._genTo);

      if ((Dlo > -tolerance) && (Dlo < tolerance) &&
          (Dhi > -tolerance) && (Dhi < tolerance)) {
        foundA[i]++;
        foundB[j]++;

        if (numSame)
          (*numSame)++;
      }
    }
  }

  for (i=0; i<A->_numExons; i++) {
    //if (foundA[i] > 1)  fprintf(stderr, "WARNING: Found exon %d %d times in A!\n", i, foundA[i]);
    if (numAMissed && (foundA[i] == 0))
      (*numAMissed)++;
  }

  for (i=0; i<B->_numExons; i++) {
    //if (foundB[i] > 1)  fprintf(stderr, "WARNING: Found exon %d %d times in B!\n", i, foundB[i]);
    if (numBMissed && (foundB[i] == 0))
      (*numBMissed)++;
  }

  delete [] foundA;
}
