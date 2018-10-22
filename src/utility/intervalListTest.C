
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
 *    Brian P. Walenz beginning on 2018-JUL-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "intervalList.H"

//  g++ -o intervalListTest -I.. -I. intervalListTest.C

int
main(int argc, char **argv) {

  if (0) {
    intervalList<int32>  t1;

    t1.add(0, 10);
    t1.add(11,7);
    t1.add(20, 8);

    fprintf(stderr, "BEFORE:\n");
    for (uint32 ii=0; ii<t1.numberOfIntervals(); ii++)
      fprintf(stderr, "%2d %3d-%3d\n", ii, t1.lo(ii), t1.hi(ii));

    t1.merge(-1);

    fprintf(stderr, "AFTER:\n");
    for (uint32 ii=0; ii<t1.numberOfIntervals(); ii++)
      fprintf(stderr, "%2d %3d-%3d\n", ii, t1.lo(ii), t1.hi(ii));
  }

  if (1) {
    intervalList<uint32>  il;

    il.add(1, -1);

    intervalList<uint32>  de(il);

    il.merge();

    for (uint32 ii=0; ii<il.numberOfIntervals(); ii++)
      fprintf(stderr, "il %2u %4u-%4u\n", ii, il.lo(ii), il.hi(ii));

    for (uint32 ii=0; ii<de.numberOfIntervals(); ii++)
      fprintf(stderr, "de %2u %4u-%4u %4d\n", ii, de.lo(ii), de.hi(ii), de.depth(ii));
  }

  exit(0);
}
