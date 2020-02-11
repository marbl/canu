
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


void
searchShiftRegisterFast(shiftRegisterParameters &srPar) {
 //  Allocate space for the loop detector, set local variables.

  uint64  sr     = srPar.getEncodedSR();      //  The shift register state 
  uint64  cyclen = srPar.getCycleLen();       //
  uint64  cycmax = srPar.getCycleMax();       //  The length of the maximum cycle
  uint64  sv     = srPar.getEncodedSVmin();   //  The tap vector
  uint64  svmax  = srPar.getEncodedSVmax();   //  The first vector we don't want to examine
  uint64  svmask = srPar.getEncodedSVmask();
  uint64  gf4widemult[4];

  //  Log.

  fprintf(stderr, "Finding cycles for length %u bases.\n", srPar.len);
  fprintf(stderr, "  sr     %022lo\n", expandTo3(sr));
  fprintf(stderr, "  cyclen %lu\n", cyclen);
  fprintf(stderr, "  cycmax %lu\n", cycmax);
  fprintf(stderr, "  sv     %022lo\n", expandTo3(sv));
  fprintf(stderr, "  svmax  %022lo\n", expandTo3(svmax));
  fprintf(stderr, "  svmask %022lo\n", expandTo3(svmask));
  fprintf(stderr, "\n");

  //  We can precompute the result of gf4mult[sv[kk]][out] used in the addition.
  //  There are four possible 'out' values, and sv is fixed for each cycle.
  //
  //  The addition itself needs an extra bit for overflow, which
  //  does make our detect array 50% larger.

  bitArray  *detect = new bitArray(cycmax);

  while (sv <= svmax) {
    sr   = 0x0000000000000001llu;

    for (uint64 out=0; out<4; out++) {   //  Build the wide multiplication table.
      uint64 mult = 0;

      //  Including all these for small sizes (e.g., size 9) doesn't
      //  seem to result in any slow down.

      mult |= gf4mult[(sv >> 42) & 0x03][out];  mult <<= 2;   // 22
      mult |= gf4mult[(sv >> 40) & 0x03][out];  mult <<= 2;   // 21
      mult |= gf4mult[(sv >> 38) & 0x03][out];  mult <<= 2;   // 20
      mult |= gf4mult[(sv >> 36) & 0x03][out];  mult <<= 2;   // 19
      mult |= gf4mult[(sv >> 34) & 0x03][out];  mult <<= 2;   // 18
      mult |= gf4mult[(sv >> 32) & 0x03][out];  mult <<= 2;   // 17
      mult |= gf4mult[(sv >> 30) & 0x03][out];  mult <<= 2;   // 16
      mult |= gf4mult[(sv >> 28) & 0x03][out];  mult <<= 2;   // 15
      mult |= gf4mult[(sv >> 26) & 0x03][out];  mult <<= 2;   // 14
      mult |= gf4mult[(sv >> 24) & 0x03][out];  mult <<= 2;   // 13
      mult |= gf4mult[(sv >> 22) & 0x03][out];  mult <<= 2;   // 12
      mult |= gf4mult[(sv >> 20) & 0x03][out];  mult <<= 2;   // 11
      mult |= gf4mult[(sv >> 18) & 0x03][out];  mult <<= 2;   // 10
      mult |= gf4mult[(sv >> 16) & 0x03][out];  mult <<= 2;   //  9
      mult |= gf4mult[(sv >> 14) & 0x03][out];  mult <<= 2;   //  8
      mult |= gf4mult[(sv >> 12) & 0x03][out];  mult <<= 2;   //  7
      mult |= gf4mult[(sv >> 10) & 0x03][out];  mult <<= 2;   //  6
      mult |= gf4mult[(sv >>  8) & 0x03][out];  mult <<= 2;   //  5
      mult |= gf4mult[(sv >>  6) & 0x03][out];  mult <<= 2;   //  4
      mult |= gf4mult[(sv >>  4) & 0x03][out];  mult <<= 2;   //  3
      mult |= gf4mult[(sv >>  2) & 0x03][out];  mult <<= 2;   //  2
      mult |= gf4mult[(sv >>  0) & 0x03][out];                //  1

      mult &= svmask;

#ifdef DEBUG
      fprintf(stderr, "widemult[%022lo][%02x] %022lo\n", expandTo3(sv), out, expandTo3(mult));
#endif

      gf4widemult[out] = mult;
    }

    //  Loop until we hit a cycle.
    //
    //  Note that cyclen will be one less than the max possible, because we
    //  cannot ever get the all zero state.  This is nice for the
    //  implementation, because we can then just compare against svmax, the
    //  last possible state vector, to decide if we cycled through all
    //  possible 'svmax + 1' states (including the all zero state we can't
    //  ever get).

    detect->clear();                     //  Reset the cycle detector.

    for (cyclen=0; (detect->flipBit(sr) == false); cyclen++) {
      uint64  out = sr & 0x03;           //  Save the output value
      uint64  mul = gf4widemult[out];    //  Compute the multiplier

#ifdef DEBUG
      fprintf(stderr, "cycle %8lu out %02lx sr %022lo\n", cyclen, out, expandTo3(sr >> 2));
      fprintf(stderr, "                     add %022lo\n", expandTo3(mul));
      fprintf(stderr, "                   final %022lo\n", expandTo3((sr >> 2) ^ mul));
#endif

      sr = (sr >> 2) ^ mul;
    }

    //  Report the cycle if it's sufficiently large.

#if 1
    if (cyclen + 1 >= 1.00 * cycmax)
      fprintf(stdout, "%12lu/%12lu %7.3f%% for vector %022lo\n",
              cyclen,
              cycmax,
              100.0 * cyclen / cycmax,
              expandTo3(sv));
#endif

    //  And move the next SV.

    sv++;
  }

  delete detect;
}


