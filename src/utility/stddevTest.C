
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
 *    Brian P. Walenz beginning on 2016-MAR-10
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "stddev.H"

//  g++ -Wall -o stddevTest -I. -I.. stddevTest.C

void
testInsert(void) {
  stdDev<uint32>   sdu;
  stdDev<int32>    sdi;
  stdDev<double>   sdd;

  sdu.insert((uint32)2);
  sdu.insert((uint32)4);
  sdu.insert((uint32)4);
  sdu.insert((uint32)4);
  sdu.insert((uint32)5);
  sdu.insert((uint32)5);
  sdu.insert((uint32)7);
  sdu.insert((uint32)9);

  sdi.insert((uint32)2);
  sdi.insert((uint32)4);
  sdi.insert((uint32)4);
  sdi.insert((uint32)4);
  sdi.insert((uint32)5);
  sdi.insert((uint32)5);
  sdi.insert((uint32)7);
  sdi.insert((uint32)9);

  sdd.insert((uint32)2);
  sdd.insert((uint32)4);
  sdd.insert((uint32)4);
  sdd.insert((uint32)4);
  sdd.insert((uint32)5);
  sdd.insert((uint32)5);
  sdd.insert((uint32)7);
  sdd.insert((uint32)9);

  fprintf(stderr, "Expect mean=5, variance=%f, stddev=%f\n", 32.0 / 7.0, sqrt(32.0 / 7.0));

  fprintf(stderr, "  uint32  size %u mean %f variance %f stddev %f\n",
          sdu.size(), sdu.mean(), sdu.variance(), sdu.stddev());
  fprintf(stderr, "  int32   size %u mean %f variance %f stddev %f\n",
          sdi.size(), sdi.mean(), sdi.variance(), sdi.stddev());
  fprintf(stderr, "  double  size %u mean %f variance %f stddev %f\n",
          sdd.size(), sdd.mean(), sdd.variance(), sdd.stddev());

  assert(sdu.variance() == 32.0 / 7.0);
  assert(sdi.variance() == 32.0 / 7.0);
  assert(sdd.variance() == 32.0 / 7.0);

  fprintf(stderr, "\n\n");
}



void
testRemove(void) {
  double   values[10] = { 1, 2, 3, 4, 9, 8, 7, 6, 20, 30 };

  stdDev<double>  sd;

  fprintf(stderr, "Expect final to be zero, and insert() == remove().\n");

  for (int ii=0; ii<10; ii++) {
    sd.insert(values[ii]);
    fprintf(stderr, "insert[%2d] mean %8.4f stddev %8.4f\n", ii+1, sd.mean(), sd.stddev());
  }

  assert(sd.mean() == 9.0);

  fprintf(stderr, "\n");

  for (int ii=9; ii>=0; ii--) {
    sd.remove(values[ii]);
    fprintf(stderr, "remove[%2d] mean %8.4f stddev %8.4f\n", ii, sd.mean(), sd.stddev());
  }

  assert(sd.mean()   == 0.0);
  assert(sd.stddev() == 0.0);
}



int
main(int argc, char **argv) {

  testInsert();
  testRemove();

  exit(0);
}
