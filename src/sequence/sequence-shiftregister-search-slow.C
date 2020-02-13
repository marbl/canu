
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
 *    src/sequence/sequence.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sequence/sequence.H"
#include "sequence/sequence-shiftregister-gf4.H"

#include "bits.H"


char    srprint[65];
char    svprint[65];


//  Print the array backwards.  [len]....[2][1][0].
char *
printSV(uint32 len, uint32 *sv) {
  uint32  pp = len;

  for (uint32 ii=0; ii<len; ii++)
    svprint[--pp] = sv[ii] + '0';

  svprint[len] = 0;

  return(svprint);
}



//  Print the array backwards.  [len]....[2][1][0].
char *
printMult(uint32 len, uint32 *sv, uint32 out) {
  uint32  pp = len;

  for (uint32 kk=0; kk<len; kk++)
    svprint[--pp] = gf4mult[sv[kk]][out] + '0';

  svprint[len] = 0;

  return(svprint);
}



bool
lessThanEqual(uint32 len, uint32 *A, uint32 *B) {
  uint32 pp = len-1;

  for (uint32 kk=0; kk<len; kk++, pp--) {
    if (A[pp] < B[pp])   //  Yup, A is less than B.
      return(true);
    if (B[pp] < A[pp])   //  Nope, A is more than B.
      return(false);
  }

  return(true);
}



//  Increment the array, [0] low order.
void
incrementSV(uint32 len, uint32 *sv) {
  uint32 pp = 0;

  sv[pp]++;

  for (uint32 ll=len; --ll > 0; ) {
    if (sv[pp] > 3) {
      sv[pp] = 0;
      sv[++pp]++;
    } else {
      break;
    }
  }
}





void
searchShiftRegisterSlow(shiftRegisterParameters &srPar) {

  //  Note.
  //
  //    [0] is the output value
  //    [i] <- [i+1]

  uint32  SR[65]    = { 0 };
  uint32  sr[65]    = { 0 };
  uint64  cyclen    =   0;
  uint64  cycmax    =   0;
  uint32  sv[65]    = { 0 };
  uint32  svmax[65] = { 0 };

  //  Make the tap vector a monic polynomial, else we're guaranteed to 
  //  never find a maximal size cycle.

  cycmax             = ((uint64)1) << (2 * srPar.order);
  SR[0]              = 1;
  sr[0]              = 1;
  sv[srPar.order-1]  = 1;

  for (uint32 ii=0; ii<srPar.order; ii++)
    svmax[ii] = 3;

  //  If we're given intiial values, use those.

  if (srPar.sr[0] != 0) { 
    for (uint32 ii=0; ii<srPar.order; ii++)
      SR[ii] = sr[ii] = srPar.sr[ii] - '0';
  }

  if (srPar.svmin[0] != 0) {
    for (uint32 ii=0; ii<srPar.order; ii++)
      sv[ii] = srPar.svmin[ii] - '0';
  }

  if (srPar.svmax[0] != 0) {
    for (uint32 ii=0; ii<srPar.order; ii++)
      svmax[ii] = srPar.svmax[ii] - '0';
  }

  //  Check that sv and sr are plausible.

  ;

  //  Log.

  fprintf(stderr, "Finding cycles for length %u bases (slow method).\n", srPar.order);
  fprintf(stderr, "  sr     %s\n", printSV(srPar.order, sr));
  fprintf(stderr, "  cyclen %lu\n", cyclen);
  fprintf(stderr, "  cycmax %lu\n", cycmax);
  fprintf(stderr, "  sv     %s\n", printSV(srPar.order, sv));
  fprintf(stderr, "  svmax  %s\n", printSV(srPar.order, svmax));
  fprintf(stderr, "\n");

  //

  bitArray  *detect = new bitArray(cycmax);
  uint64     kmer   = 0;

  while (lessThanEqual(srPar.order, sv, svmax) == true) {
    //fprintf(stderr, "SV %s\r", printSV(srPar.order, sv));

    //  Reset the shift register, rebuild the kmer.

    kmer = 0;

    for (uint32 ii=0; ii<srPar.order; ii++) {
      sr[ii] =  SR[ii];

      kmer <<= 2;
      kmer  |= sr[ii];
    }

    //  Loop until we hit a cycle.

    detect->clear();

    for (cyclen=1; (detect->flipBit(kmer) == false); cyclen++) {
      uint32  out = sr[0];

#ifdef DEBUG
      fprintf(stderr, "cycle %8lu out %c sr %s\n", cyclen, out + '0', printSV(srPar.order, sr + 1));
      fprintf(stderr, "                    add %s\n", printMult(srPar.order, sv, out));
#endif

      //  Shift, add in the taps, and rebuild the kmer

      kmer = 0;

      for (uint32 kk=0; kk<srPar.order; kk++) {
        sr[kk] = sr[kk+1] ^ gf4mult[sv[kk]][out];

        kmer <<= 2;
        kmer  |= sr[kk];
      }

#ifdef DEBUG
      fprintf(stderr, "                  final %s\n", printSV(srPar.order, sr));
#endif
    }

    //  report the cycle.
    if (cyclen + 1 >= srPar.report * cycmax)
      fprintf(stdout, "%14lu/%14lu %7.3f%% for vector %s\n",
              cyclen,
              cycmax,
              100.0 * cyclen / cycmax,
              printSV(srPar.order, sv));

    incrementSV(srPar.order, sv);
  }

  fprintf(stderr, "SV %s\n", printSV(srPar.order, sv));

  delete detect;
}
