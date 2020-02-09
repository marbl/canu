
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
#include "bits.H"

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





void
searchShiftRegister_9(shiftRegisterParameters &srPar) {
  //  Allocate space for the loop detector, set local variables.

  uint64  kmer   = 0;
  uint32  len    = 9;
  uint32  SR[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 };
  uint32  sr[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  uint32  sv[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  uint32  shift  = 2;
  uint64  kmask  = 0x0003ffff;
  uint32  mask   = 0x00000003;
  uint32  maxv   = 3;

  bitArray  *detect = new bitArray(kmask + 1);

  while (sv[0] <= maxv) {

    //  Reset the string, rebuild the kmer.

    kmer = 0;

    for (uint32 ii=0; ii<9; ii++) {
      sr[ii] = SR[ii];

      kmer <<= shift;
      kmer  |= sr[ii];
    }

    //  Reset the cycle detector.

    detect->clear();

    //  Loop until we hit a cycle.

    uint32  cycleLen = 1;

    while (detect->flipBit(kmer) == false) {
      cycleLen++;

      //  Save the output.

      uint32  out = sr[0];

      //  Do a left shift, adding to each tap (if sv[kk] is zero, the multiply is zero).
      //  Note that this requires sr[len] == 0.
      //
      //  At the same time, rebuild the kmer.

      kmer = 0;

      for (uint32 kk=0; kk<9; kk++) {
        sr[kk]  = sr[kk+1];
        sr[kk] += gf4mult[sv[kk]][out];
        sr[kk] &= mask;

        kmer <<= shift;
        kmer  |= sr[kk];
      }
    }

    //  Report the cycle if it's sufficiently large.

    if (cycleLen >= kmask + 1)
      fprintf(stderr, "%05u %s 0x%04lx - cycle for vector %s\n",
              cycleLen, printSR(9, sr), kmer, printSV(9, sv));

    //  And move onto the next SV.

    incrementSV(9, maxv, sv);
  }

  delete detect;
}









void
searchShiftRegister_9bits(shiftRegisterParameters &srPar) {
  //  Allocate space for the loop detector, set local variables.

  uint64  sr     = 0x0000000000000001llu;  //  This is 3-bit packed
  uint64  sv     = 0x0000000000000000llu;  //  This is 2-bit packed
  uint64  svmax  = 0x000000000003ffffllu;  //  This is 2-bit packed
  uint64  kmer   = 0x0000000000000001llu;  //  This is 2-bit packed

  bitArray  *detect = new bitArray(svmax + 1);

  //  We can precompute the result of gf4mult[sv[kk]][out] used in the addition.
  //  There are four possible 'out' values, and sv is fixed for each cycle.
  //
  //  The addition itself needs an extra bit for overflow, which
  //  does make our detect array 50% larger.

  uint64  gf4widemult[4];
  uint64  cycleLen;

  while (sv <= svmax) {
    kmer = 0x0000000000000001llu;        //  Reset the string, rebuild the kmer.
    sr   = 0x0000000000000001llu;

    detect->clear();                     //  Reset the cycle detector.

    for (uint32 out=0; out<4; out++) {   //  Build the wide multiplication table.
      uint64 mult = 0;

      mult  = gf4mult[(sv >> 16) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 14) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 12) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 10) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  8) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  6) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  4) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  2) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  0) & 0x03][out];

      gf4widemult[out] = mult;
    }

    //  Loop until we hit a cycle.
    //
    //  Note that cycleLen will be one less than the max possible, because we
    //  cannot ever get the all zero state.  This is nice for the
    //  implementation, because we can then just compare against svmax, the
    //  last possible state vector, to decide if we cycled through all
    //  possible 'svmax + 1' states (including the all zero state we can't
    //  ever get).

    for (cycleLen=0; (detect->flipBit(kmer) == false); cycleLen++) {
      uint64  out = sr & 0x03;                  //  Save the output value

      sr >>= 3;                                 //  Shift right.
      sr   += gf4widemult[out];                 //  Add in the taps.
      sr   &= 0333333333333333333333llu;        //  Mask out the extra carries (octal = 011011...)

      kmer  = (sr & 0x0000000000000003) >> 0;   //  Compress the 3-bit wide sr into a 2-bit wide kmer
      kmer |= (sr & 0x0000000000000018) >> 1;
      kmer |= (sr & 0x00000000000000c0) >> 2;
      kmer |= (sr & 0x0000000000000600) >> 3;
      kmer |= (sr & 0x0000000000003000) >> 4;
      kmer |= (sr & 0x0000000000018000) >> 5;
      kmer |= (sr & 0x00000000000c0000) >> 6;
      kmer |= (sr & 0x0000000000600000) >> 7;
      kmer |= (sr & 0x0000000003000000) >> 8;
    }

    //  Report the cycle if it's sufficiently large.

    if (cycleLen >= svmax)
      fprintf(stderr, "%-8lu 0x%016lx - cycle for vector 0x%016lx\n",
              cycleLen, sr, sv);

    //  And move the next SV.

    sv++;
  }

  delete detect;
}




void
searchShiftRegister_10bits(shiftRegisterParameters &srPar) {
  //  Allocate space for the loop detector, set local variables.

  uint64  sr     = 0x0000000000000001llu;  //  This is 3-bit packed
  uint64  sv     = 0x0000000000000000llu;  //  This is 2-bit packed
  uint64  svmax  = 0x00000000000fffffllu;  //  This is 2-bit packed
  uint64  kmer   = 0x0000000000000001llu;  //  This is 2-bit packed

  bitArray  *detect = new bitArray(svmax + 1);

  //  We can precompute the result of gf4mult[sv[kk]][out] used in the addition.
  //  There are four possible 'out' values, and sv is fixed for each cycle.
  //
  //  The addition itself needs an extra bit for overflow, which
  //  does make our detect array 50% larger.

  uint64  gf4widemult[4];
  uint64  cycleLen;

  while (sv <= svmax) {
    kmer = 0x0000000000000001llu;        //  Reset the string, rebuild the kmer.
    sr   = 0x0000000000000001llu;

    detect->clear();                     //  Reset the cycle detector.

    for (uint32 out=0; out<4; out++) {   //  Build the wide multiplication table.
      uint64 mult = 0;

      mult  = gf4mult[(sv >> 18) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 16) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 14) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 12) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 10) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  8) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  6) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  4) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  2) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  0) & 0x03][out];

      gf4widemult[out] = mult;
    }

    //  Loop until we hit a cycle.
    //
    //  Note that cycleLen will be one less than the max possible, because we
    //  cannot ever get the all zero state.  This is nice for the
    //  implementation, because we can then just compare against svmax, the
    //  last possible state vector, to decide if we cycled through all
    //  possible 'svmax + 1' states (including the all zero state we can't
    //  ever get).

    for (cycleLen=0; (detect->flipBit(kmer) == false); cycleLen++) {
      uint64  out = sr & 0x03;                  //  Save the output value

      sr >>= 3;                                 //  Shift right.
      sr   += gf4widemult[out];                 //  Add in the taps.
      sr   &= 0333333333333333333333llu;        //  Mask out the extra carries (octal = 011011...)

      kmer  = (sr & 0x0000000000000003) >> 0;   //  Compress the 3-bit wide sr into a 2-bit wide kmer
      kmer |= (sr & 0x0000000000000018) >> 1;
      kmer |= (sr & 0x00000000000000c0) >> 2;
      kmer |= (sr & 0x0000000000000600) >> 3;
      kmer |= (sr & 0x0000000000003000) >> 4;
      kmer |= (sr & 0x0000000000018000) >> 5;
      kmer |= (sr & 0x00000000000c0000) >> 6;
      kmer |= (sr & 0x0000000000600000) >> 7;
      kmer |= (sr & 0x0000000003000000) >> 8;
      kmer |= (sr & 0x0000000018000000) >> 9;
    }

    //  Report the cycle if it's sufficiently large.

    if (cycleLen >= svmax)
      fprintf(stderr, "%-8lu 0x%016lx - cycle for vector 0x%016lx\n",
              cycleLen, sr, sv);

    //  And move the next SV.

    sv++;
  }

  delete detect;
}





void
searchShiftRegister_11bits(shiftRegisterParameters &srPar) {
  //  Allocate space for the loop detector, set local variables.

  uint64  sr     = 0x0000000000000001llu;  //  This is 3-bit packed
  uint64  sv     = 0x0000000000000000llu;  //  This is 2-bit packed
  uint64  svmax  = 0x00000000003fffffllu;  //  This is 2-bit packed
  uint64  kmer   = 0x0000000000000001llu;  //  This is 2-bit packed

  //  Make the polynomial monic, else we're guaranteed to
  //  be not a maximal size cycle.
  sv  = 1 << 20;
  sv &= svmax;

  assert(sv != 0);

  bitArray  *detect = new bitArray(svmax + 1);

  //  We can precompute the result of gf4mult[sv[kk]][out] used in the addition.
  //  There are four possible 'out' values, and sv is fixed for each cycle.
  //
  //  The addition itself needs an extra bit for overflow, which
  //  does make our detect array 50% larger.

  uint64  gf4widemult[4];
  uint64  cycleLen;

  while (sv <= svmax) {
    kmer = 0x0000000000000001llu;        //  Reset the string, rebuild the kmer.
    sr   = 0x0000000000000001llu;

    detect->clear();                     //  Reset the cycle detector.

    for (uint32 out=0; out<4; out++) {   //  Build the wide multiplication table.
      uint64 mult = 0;

      mult  = gf4mult[(sv >> 20) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 18) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 16) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 14) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 12) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >> 10) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  8) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  6) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  4) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  2) & 0x03][out];  mult <<= 3;
      mult |= gf4mult[(sv >>  0) & 0x03][out];

      gf4widemult[out] = mult;

      //fprintf(stderr, "widemult[%016lx][%02x] %022lo\n", sv, out, mult);
    }

    //  Loop until we hit a cycle.
    //
    //  Note that cycleLen will be one less than the max possible, because we
    //  cannot ever get the all zero state.  This is nice for the
    //  implementation, because we can then just compare against svmax, the
    //  last possible state vector, to decide if we cycled through all
    //  possible 'svmax + 1' states (including the all zero state we can't
    //  ever get).

    for (cycleLen=0; (detect->flipBit(kmer) == false); cycleLen++) {
      uint64  out = sr & 0x03;                  //  Save the output value

      sr >>= 3;                                 //  Shift right.

      //fprintf(stderr, "cycle %8lu out %02lx sr %022lo\n", cycleLen, out, sr);

      sr   ^= gf4widemult[out];                 //  Add in the taps.
      sr   &= 0333333333333333333333llu;        //  Mask out the extra carries (octal = 011011...)

      //fprintf(stderr, "                     add %022lo\n", gf4widemult[out]);
      //fprintf(stderr, "                   final %022lo\n", sr);

      kmer  = (sr & 0x0000000000000003) >> 0;   //  Compress the 3-bit wide sr into a 2-bit wide kmer
      kmer |= (sr & 0x0000000000000018) >> 1;
      kmer |= (sr & 0x00000000000000c0) >> 2;
      kmer |= (sr & 0x0000000000000600) >> 3;
      kmer |= (sr & 0x0000000000003000) >> 4;
      kmer |= (sr & 0x0000000000018000) >> 5;
      kmer |= (sr & 0x00000000000c0000) >> 6;
      kmer |= (sr & 0x0000000000600000) >> 7;
      kmer |= (sr & 0x0000000003000000) >> 8;
      kmer |= (sr & 0x0000000018000000) >> 9;
      kmer |= (sr & 0x00000000c0000000) >> 10;
    }

    //  Report the cycle if it's sufficiently large.

    if (cycleLen >= 0.99 * svmax)
      fprintf(stderr, "%8lu/%8lu %7.3f%% for vector 0%016lx\n",
              cycleLen,
              svmax,
              100.0 * cycleLen / svmax,
              sv);

    //  And move the next SV.

    sv++;
  }

  delete detect;
}






void
searchShiftRegister_10(shiftRegisterParameters &srPar) {
  //  Allocate space for the loop detector, set local variables.

  uint64  kmer   = 0;
  uint32  len    = 10;
  uint32  SR[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 };
  uint32  sr[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  uint32  sv[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  uint32  shift  = 2;
  uint64  kmask  = 0x000fffff;
  uint32  mask   = 0x00000003;
  uint32  maxv   = 3;

  bitArray  *detect = new bitArray(kmask + 1);

  while (sv[0] <= maxv) {

    //  Reset the string, rebuild the kmer.

    kmer = 0;

    for (uint32 ii=0; ii<10; ii++) {
      sr[ii] = SR[ii];

      kmer <<= shift;
      kmer  |= sr[ii];
    }

    //  Reset the cycle detector.

    detect->clear();

    //  Loop until we hit a cycle.

    uint32  cycleLen = 1;

    while (detect->flipBit(kmer) == false) {
      cycleLen++;

      //  Save the output.

      uint32  out = sr[0];

      //  Do a left shift, adding to each tap (if sv[kk] is zero, the multiply is zero).
      //  Note that this requires sr[len] == 0.
      //
      //  At the same time, rebuild the kmer.

      kmer = 0;

      for (uint32 kk=0; kk<10; kk++) {
        sr[kk]  = sr[kk+1];
        sr[kk] += gf4mult[sv[kk]][out];
        sr[kk] &= mask;

        kmer <<= shift;
        kmer  |= sr[kk];
      }
    }

    //  Report the cycle if it's sufficiently large.

    if (cycleLen >= kmask + 1)
      fprintf(stderr, "%05u %s 0x%04lx - cycle for vector %s\n",
              cycleLen, printSR(10, sr), kmer, printSV(10, sv));

    //  And move onto the next SV.

    incrementSV(10, maxv, sv);
  }

  delete detect;
}






void
doShiftRegister(shiftRegisterParameters &srPar) {

  srPar.initialize();

  fprintf(stderr, "VERSION 2\n");

#if 0
  if (srPar.len == 9) {
    searchShiftRegister_9bits(srPar);
    return;
  }

  if (srPar.len == 10) {
    searchShiftRegister_10bits(srPar);
    return;
  }

  if (srPar.len == 11) {
    searchShiftRegister_11bits(srPar);
    return;
  }
#endif

  //  Allocate space for the loop detector, set local variables.

  uint64   kmer  = 0;
  uint32  len    = srPar.len;
  uint32  SR[17] = {0};
  uint32  sr[17] = {0};
  uint32  sv[17] = {0};

  uint32  shift  = 2;
  uint64  kmask  = 0;
  uint32  mask   = 0x0003;
  uint32  maxv   = 3;

  uint64  svmax  = 1;

  for (uint32 ii=0; ii<len; ii++) {
    SR[ii] = srPar.sr[ii] - '0';
    sr[ii] = srPar.sr[ii] - '0';
    sv[ii] = srPar.sv[ii] - '0';

    kmask <<= shift;
    kmask  |= maxv;

    svmax *= 4;
  }

  sv[0] = 1;   //  Else we're trivially not a maximal cycle.

  bitArray  *detect = new bitArray(kmask + 1);

#define SEQOUT
#ifdef SEQOUT
  fprintf(stderr, "Allocating space for a sequence of length %lu\n", kmask);
  char      *seq    = new char    [kmask];
#endif

#if 1
  while (sv[0] <= maxv) {
#endif

    for (uint32 ii=0; ii<len; ii++) {    //  RESET the string, rebuild the kmer
      sr[ii] =  SR[ii];

      kmer <<= shift;
      kmer  |= sr[ii];
      kmer  &= kmask;
    }

    //fprintf(stderr, "SR %s\n", printSV(len, sr));   //  prints as 0123
    fprintf(stderr, "SV %s\n", printSV(len, sv));

    detect->clear();

    //  Loop until we hit a cycle.

    uint64  cycleLen = 1;

    while (detect->flipBit(kmer) == false) {
      cycleLen++;

      uint32  out = sr[0];

#ifdef SEQOUT
      seq[cycleLen-2] = printBase(out);
      seq[cycleLen-1] = 0;
#endif

      //  Rotate and shift in a new zero.  Requres sr[len] = 0;
      for (uint32 kk=0; kk<len; kk++)
        sr[kk]  = sr[kk+1];

      //fprintf(stdout, "%c", printBase(out));

      //fprintf(stderr, "cycle %8lu out %c sr %s\n", cycleLen, out + '0', printSV(len, sr));
      //fprintf(stderr, "                    add ");

      //  Add in the taps.
      for (uint32 kk=0; kk<len; kk++) {
        //fprintf(stderr, "%c", gf4mult[sv[kk]][out] + '0');

        sr[kk] ^= gf4mult[sv[kk]][out];
        sr[kk] &= mask;
      }

      //fprintf(stderr, "\n");
      //fprintf(stderr, "                  final %s\n", printSV(len, sr));

      //  Rebuild the kmer
      for (uint32 kk=0; kk<len; kk++) {
        kmer <<= shift;
        kmer  |= sr[kk];
        kmer  &= kmask;
      }
    }

    //  report the cycle.
    fprintf(stderr, "%8lu/%8lu %7.3f%% for vector %s\n",
            cycleLen,
            svmax,
            100.0 * cycleLen / svmax,
            printSV(len, sv));

#ifdef SEQOUT
    if (cycleLen >= 0.99 * svmax) {
      fprintf(stdout, ">%s\n", printSV(len, sv));
      fprintf(stdout, "%s\n", seq);
      fflush(stdout);
    }
#endif

    //fprintf(stdout, "\n");

#if 1
    incrementSV(len, maxv, sv);
  }
#endif

  delete detect;
}


