
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
emitShiftRegisterFast(shiftRegisterParameters &srPar) {
 //  Allocate space for the loop detector, set local variables.

  uint64  sr     = srPar.getEncodedSR();      //  The shift register state 
  uint64  cyclen = srPar.getCycleLen();       //
  uint64  cycmax = srPar.getCycleMax();       //  The length of the maximum cycle
  uint64  sv     = srPar.getEncodedSVmin();   //  The tap vector
  uint64  svmax  = srPar.getEncodedSVmax();   //  The first vector we don't want to examine
  uint64  svmask = srPar.getEncodedSVmask();
  uint64  gf4widemult[4];

  //  If no length supplied, set it so we'll emit every possible kmer (except
  //  the all-zero kmer that we manually insert).

  if (srPar.length == 0)
    srPar.length = cycmax + srPar.order - 1 - 1;

  //  Log.

  fprintf(stderr, "Emitting %lu bases, plus an extra %u.\n", srPar.length - srPar.order - 1 - 1, srPar.order - 1);
  fprintf(stderr, "  sr     %0*lo\n", srPar.order, expandTo3(sr));
  fprintf(stderr, "  cyclen %lu\n",   cyclen);
  fprintf(stderr, "  cycmax %lu\n",   cycmax);
  fprintf(stderr, "  sv     %0*lo\n", srPar.order, expandTo3(sv));
  fprintf(stderr, "  svmax  %0*lo\n", srPar.order, expandTo3(svmax));
  fprintf(stderr, "  svmask %0*lo\n", srPar.order, expandTo3(svmask));
  fprintf(stderr, "\n");

  //  Write a header.

  fprintf(stdout, ">bases\n");

  //  Precompute the result of gf4mult[sv[kk]][out] used in the addition.
  //  There are four possible 'out' values, and sv is fixed for each cycle.

  for (uint64 out=0; out<4; out++) {
    uint64 mult = 0;

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

    gf4widemult[out] = mult & svmask;
  }

  //  Loop until we have emitted the requested number of bases.

  uint32  bufferMax = 100;
  uint32  bufferLen = 0;
  char   *buffer    = new char [bufferMax];

  bool    emitExtra = true;
  uint32  emitCount = 1;

  for (cyclen=0; cyclen < srPar.length; cyclen++) {
    uint64  out = sr & 0x03;           //  Save the output value
    uint64  mul = gf4widemult[out];    //  Compute the multiplier

    //  Emit the base.

    buffer[bufferLen++] = srPar.numberToBase(out);

    if (bufferLen == bufferMax) {
      writeToFile(buffer, "bases", sizeof(char), bufferLen, stdout);

      bufferLen = 0;
    }

    //  If this is the first time we've seen order-1 zero's in a row, output
    //  an extra one to make the all-zero kmer.

    if (out == 0)
      emitCount++;
    else
      emitCount = 1;

    if ((emitExtra == true) && (emitCount == srPar.order)) {
      emitExtra = false;
      buffer[bufferLen++] = srPar.numberToBase(0);
    }

    //  Log.

#ifdef DEBUG
    fprintf(stderr, "cycle %8lu out %02lx sr %0*lo\n", cyclen, out, srPar.order, expandTo3(sr >> 2));
    fprintf(stderr, "                     add %0*lo\n", srPar.order, expandTo3(mul));
    fprintf(stderr, "                   final %0*lo\n", srPar.order, expandTo3((sr >> 2) ^ mul));
#endif

    //  And do the work.

    sr = (sr >> 2) ^ mul;
  }

  writeToFile(buffer, "bases", sizeof(char), bufferLen, stdout);
  fprintf(stdout, "\n");

  delete [] buffer;
}


