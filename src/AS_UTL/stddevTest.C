
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2015-MAR-10
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "stddev.H"

//  g++ -Wall -o stddevTest -I. -I.. stddevTest.C

int
main(int argc, char **argv) {

  stdDev<uint32>   sdu;
  stdDev<int32>    sdi;
  stdDev<double>   sdd;

  sdu.update((uint32)2);
  sdu.update((uint32)4);
  sdu.update((uint32)4);
  sdu.update((uint32)4);
  sdu.update((uint32)5);
  sdu.update((uint32)5);
  sdu.update((uint32)7);
  sdu.update((uint32)9);

  sdi.update((uint32)2);
  sdi.update((uint32)4);
  sdi.update((uint32)4);
  sdi.update((uint32)4);
  sdi.update((uint32)5);
  sdi.update((uint32)5);
  sdi.update((uint32)7);
  sdi.update((uint32)9);

  sdd.update((uint32)2);
  sdd.update((uint32)4);
  sdd.update((uint32)4);
  sdd.update((uint32)4);
  sdd.update((uint32)5);
  sdd.update((uint32)5);
  sdd.update((uint32)7);
  sdd.update((uint32)9);


  fprintf(stderr, "uint32  size %u mean %f variance %f stddev %f\n",
          sdu.size(), sdu.mean(), sdu.variance(), sdu.stddev());
  fprintf(stderr, "int32   size %u mean %f variance %f stddev %f\n",
          sdi.size(), sdi.mean(), sdi.variance(), sdi.stddev());
  fprintf(stderr, "double  size %u mean %f variance %f stddev %f\n",
          sdd.size(), sdd.mean(), sdd.variance(), sdd.stddev());

  assert(sdu.variance() == 32.0 / 7.0);
  assert(sdi.variance() == 32.0 / 7.0);
  assert(sdd.variance() == 32.0 / 7.0);

  exit(0);
}
