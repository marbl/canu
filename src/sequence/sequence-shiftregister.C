
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


uint32  gf4mult[4][4] = { { 0, 0, 0, 0 },
                          { 0, 1, 2, 3 },
                          { 0, 2, 3, 1 },
                          { 0, 3, 1, 2 } };

uint32  gf4add[4][4] =  { { 0, 1, 2, 3 },
                          { 1, 0, 3, 2 },
                          { 2, 3, 0, 1 },
                          { 3, 2, 1, 0 } };

char  srprint[17];
char  svprint[17];

char
printBase(uint32 sr) {

  sr = (sr << 1) + 'A';

  if (sr == 'E')
    sr = 'T';

  return(sr);
}

char *
printSR(uint32 len, uint32 *sr) {

  for (uint32 ii=0; ii<len; ii++) {
    srprint[ii] = (sr[ii] << 1) + 'A';

    if (sr[ii] == 'E')
      sr[ii] = 'T';
  }

  srprint[len] = 0;

  return(srprint);
}

char *
printSV(uint32 len, uint32 *sv) {
  for (uint32 ii=0; ii<len; ii++)
    svprint[ii] = sv[ii] + '0';

  svprint[len] = 0;

  return(svprint);
}





void
incrementSV(uint32 len, uint32 maxv, uint32 *sv) {

  sv[--len]++;

  while (len > 0) {
    if (sv[len] > maxv) {
      sv[len] = 0;
      sv[len-1]++;
    } else {
      break;
    }

    len--;
  }
}

#if 0
    sv[7]++;

    if (sv[7] > maxv)  { sv[7] = 0;  sv[6]++; }
    if (sv[6] > maxv)  { sv[6] = 0;  sv[5]++; }
    if (sv[5] > maxv)  { sv[5] = 0;  sv[4]++; }
    if (sv[4] > maxv)  { sv[4] = 0;  sv[3]++; }
    if (sv[3] > maxv)  { sv[3] = 0;  sv[2]++; }
    if (sv[2] > maxv)  { sv[2] = 0;  sv[1]++; }
    if (sv[1] > maxv)  { sv[1] = 0;  sv[0]++; }
#endif




void
doShiftRegister(shiftRegisterParameters &srPar) {

  srPar.initialize();

  fprintf(stderr, "VERSION 2\n");

  //  Allocate space for the loop detector, set local variables.
  uint32  *detect = new uint32 [65536];

  uint64   kmer  = 0;
  uint32  len    = srPar.len;
  uint32  SR[17] = {0};
  uint32  sr[17] = {0};
  uint32  sv[17] = {0};

  uint32  shift  = 0;
  uint64  kmask  = 0;
  uint32  mask   = 0;
  uint32  maxv   = 0;

#if 0
  shift  = 1;
  mask   = 0x0001;
  maxv   = 1;
#else
  shift  = 2;
  mask   = 0x0003;
  maxv   = 3;
#endif

  for (uint32 ii=0; ii<len; ii++) {
    SR[ii] = srPar.sr[ii] - '0';
    sr[ii] = srPar.sr[ii] - '0';
    sv[ii] = srPar.sv[ii] - '0';

    kmask <<= shift;
    kmask  |= maxv;
  }

#if 1
  while (sv[0] <= maxv) {
#endif

    for (uint32 ii=0; ii<len; ii++) {    //  RESET the string, rebuild the kmer
      sr[ii] =  SR[ii];

      kmer <<= shift;
      kmer  |= sr[ii];
      kmer  &= kmask;
    }

    memset(detect, 0, sizeof(uint32) * 65536);   //  RESET detect

    fprintf(stdout, ">%s_%s\n", printSR(len, sr), printSV(len, sv));


    //  Loop until we hit a cycle.

    uint32  cycleLen = 1;

    while (detect[kmer] == 0) {
      assert(kmer < 65536);

#if 0
      fprintf(stderr, "%05u %c%c%c%c%c%c%c%c 0x%04lx\n",
              cycleLen, printSR(len, sr), kmer);
#endif

      //  Mark we've been here.
      detect[kmer] = cycleLen++;

      //  Save the output letter
      uint32  out = sr[0];

      fprintf(stdout, "%c", printBase(out));

      //  left circular shift, adding to the taps
      //  the non-tap cells get just the cell before
      //  the tap cells get the sum of the cell before it and the output
      for (uint32 kk=0; kk<len; kk++) {   //  REQUIRES sr[len] = 0
        sr[kk]  = sr[kk+1];
        sr[kk] += gf4mult[sv[kk]][out];
        sr[kk] &= mask;
      }

      //  Rebuild the kmer
      for (uint32 kk=0; kk<len; kk++) {
        kmer <<= shift;
        kmer  |= sr[kk];
        kmer  &= kmask;
      }
    }

    //  report the cycle.
    if (cycleLen >= kmask - 2)
      fprintf(stderr, "%05u %s 0x%04lx %5u - cycle for vector %s\n",
              cycleLen, printSR(len, sr), kmer, detect[kmer], printSV(len, sv));

    fprintf(stdout, "\n");

#if 1
    incrementSV(len, maxv, sv);
  }
#endif

  delete [] detect;
}


