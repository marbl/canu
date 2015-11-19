
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
 *    Brian P. Walenz on 2004-MAY-06
 *      are Copyright 2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2006-OCT-23 to 2014-APR-11
 *      are Copyright 2006,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-OCT-14
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>

#include "util++.H"

mt_s  *mt = 0L;

void
test(void) {
  int           e = 0;
  intervalList<uint32>  I;

  I.add(71,  3);
  I.add( 5,  3);
  I.add(32,  5);
  I.add(73,  3);
  I.add(55, 10);
  I.add( 5,  3);
  I.add(10,  5);
  I.add(20, 10);
  I.add(30, 10);
  I.add(50, 10);
  I.add(70,  3);
  I.add(72,  3);
  I.add( 5,  3);
  I.add(15,  5);

#if 0
  for (uint32 i=0; i<I.numberOfIntervals(); i++)
    fprintf(stderr, "%2d] %2d %2d\n", i, I.lo(i), I.hi(i));
#endif

  I.sort();
  I.merge();

  if (I.sumOfLengths() != 54)
    fprintf(stderr, "Lengths don't add up.\n"), e++;

  if (I.numberOfIntervals() != 4)
    fprintf(stderr, "Wrong number of intervals.\n"), e++;

  if ((I.lo(0) != 5) || (I.hi(0) != 8))
    fprintf(stderr, "Interval 0 is wrong.\n"), e++;
  if ((I.lo(1) != 10) || (I.hi(1) != 40))
    fprintf(stderr, "Interval 1 is wrong.\n"), e++;
  if ((I.lo(2) != 50) || (I.hi(2) != 65))
    fprintf(stderr, "Interval 2 is wrong.\n"), e++;
  if ((I.lo(3) != 70) || (I.hi(3) != 76))
    fprintf(stderr, "Interval 3 is wrong.\n"), e++;

  if (e)
    exit(e);
}



void
testIntersect(uint32 type) {
  uint32   numTests = 1000000;
  uint32  *beg     = new uint32 [numTests];
  uint32  *len     = new uint32 [numTests];
  uint32  *end     = new uint32 [numTests];
  uint32  *abegh   = new uint32 [numTests];
  uint32  *aendh   = new uint32 [numTests];
  uint32  *bbegh   = new uint32 [numTests];
  uint32  *bendh   = new uint32 [numTests];
  uint32   errors  = 0;
  uint32   passed  = 0;

  intervalList<uint32>  A;
  intervalList<uint32>  B;

  //
  //  Build two interval lists
  //
  //  type == 0 --> all pairwise
  //  type == 1 --> A sequence is solid
  //  type == 2 --> B sequence is solid
  //

  if (type == 1)
    A.add(1, 1500000000);
  if (type == 2)
    B.add(1, 1500000000);

  for (uint32 i=0; i<numTests; i++) {

    //  Compute the result we want to get
    //
    len[i] = mtRandom32(mt) % 200;
    if (len[i] < 100) {
      beg[i] = end[i] = mtRandom32(mt) % 100 + 100;
    } else {
      beg[i] = mtRandom32(mt) % 100 + 100;
      end[i] = beg[i] + len[i];
    }

    //  Reset if the type is 1 or 2.
    //
    if ((type == 1) || (type == 2)) {
      len[i] = mtRandom32(mt) % 100 + 100;
      beg[i] = mtRandom32(mt) % 100 + 100;
      end[i] = beg[i] + len[i];
    }

    //  Extend it to an interval -- we can extend exactly one end, or
    //  two opposite ends.
    //
    abegh[i] = 0;
    aendh[i] = 0;
    bbegh[i] = 0;
    bendh[i] = 0;

    if (type == 0) {
      switch (mtRandom32(mt) % 8) {
        case 0:
          abegh[i] = mtRandom32(mt) % 50;
          break;
        case 1:
          aendh[i] = mtRandom32(mt) % 50;
          break;
        case 2:
          bbegh[i] = mtRandom32(mt) % 50;
          break;
        case 3:
          bendh[i] = mtRandom32(mt) % 50;
          break;
        case 4:
          abegh[i] = mtRandom32(mt) % 50;
          aendh[i] = mtRandom32(mt) % 50;
          break;
        case 5:
          bbegh[i] = mtRandom32(mt) % 50;
          bendh[i] = mtRandom32(mt) % 50;
          break;
        case 6:
          abegh[i] = mtRandom32(mt) % 50;
          bendh[i] = mtRandom32(mt) % 50;
          break;
        case 7:
          aendh[i] = mtRandom32(mt) % 50;
          bbegh[i] = mtRandom32(mt) % 50;
          break;
      }
    }

    //  Add it to the lists -- if type == 1 or 2, these should then
    //  get merged into the one big thing.
    //
    A.add(1000 * i + beg[i] - abegh[i], abegh[i] + end[i] - beg[i] + aendh[i]);
    B.add(1000 * i + beg[i] - bbegh[i], bbegh[i] + end[i] - beg[i] + bendh[i]);
  }

  intervalList<uint32> I;
  I.intersect(A, B);

  //
  //  Check the result.
  //

  for (uint32 i=0, j=0; i<numTests; i++) {
    uint32  b = I.lo(j) - 1000 * i;
    uint32  e = I.hi(j) - 1000 * i;

    if (len[i] < 100) {
      //
      //  Expect no result here.  We ca only test that the stuff that
      //  should be intersecting is correct, and if all that is, then
      //  I guess the non-intersection stuff is correct too.
      //
    } else {
      if ((b != beg[i]) || (e != end[i])) {
        fprintf(stderr, "FAILED[%4d]: "uint32FMT"-"uint32FMT" X "uint32FMT"-"uint32FMT" -> "uint32FMT","uint32FMT" ("uint32FMT","uint32FMT") (should have been "uint32FMT","uint32FMT")\n",
                i,
                beg[i] - abegh[i], beg[i] - abegh[i] + abegh[i] + end[i] - beg[i] + aendh[i],
                beg[i] - bbegh[i], beg[i] - bbegh[i] + bbegh[i] + end[i] - beg[i] + bendh[i],
                b, e, (uint32)I.lo(j), (uint32)I.hi(j),
                beg[i], end[i]);
        errors++;
      } else {
        passed++;
      }
      j++;
    }
  }

  fprintf(stderr, "intersection test had "uint32FMT" successes and "uint32FMT" errors.\n", passed, errors);
}




void
testMerge(void) {
  intervalList<uint32,double>  IL;
  intervalList<uint32,double>  ID;

  //  Test 1:  one long sequence containing lots of little non-overlapping sequences
  //  Test 2:  three long overlapping sequences, containing lots of non-overlapping sequences
  //  Test 3:  dense random
  //  Test 4:  special cases

  fprintf(stderr, "Merge test 1\n");
  IL.clear();
  IL.add(0, 100000);
  for (uint32 i=0; i<999; i++)
    IL.add(100 + 100 * i, 50);
  IL.merge();
  for (uint32 i=0; i<IL.numberOfIntervals(); i++)
    fprintf(stderr, "IL["uint32FMTW(3)"] %6u %6u\n", i, IL.lo(i), IL.hi(i));

  IL.clear();
  for (uint32 i=0; i<999; i++)
    IL.add(100 + 100 * i, 50);
  IL.add(0, 100000);
  IL.merge();
  for (uint32 i=0; i<IL.numberOfIntervals(); i++)
    fprintf(stderr, "IL["uint32FMTW(3)"] %6u %6u\n", i, IL.lo(i), IL.hi(i));

  fprintf(stderr, "Merge test 2\n");
  IL.clear();
  IL.add(0, 25000);
  IL.add(25000, 25000);
  IL.add(50000, 50000);
  for (uint32 i=0; i<999; i++)
    IL.add(100 + 100 * i, 50);
  IL.merge();
  for (uint32 i=0; i<IL.numberOfIntervals(); i++)
    fprintf(stderr, "IL["uint32FMTW(3)"] %6u %6u\n", i, IL.lo(i), IL.hi(i));

  fprintf(stderr, "Merge test 3\n");
  IL.clear();
  uint32 lo = 200;
  uint32 hi = 0;
  for (uint32 i=0; i<999; i++) {
    uint32 beg = mtRandom32(mt) % 100;
    uint32 end = mtRandom32(mt) % 100 + 100;
    if (beg < lo)  lo = beg;
    if (end > hi)  hi = end;
    IL.add(beg, end - beg);
  }
  IL.merge();
  if ((IL.lo(0) != lo) || (IL.hi(0) != hi))
    fprintf(stderr, "ERROR!\n");
  for (uint32 i=0; i<IL.numberOfIntervals(); i++)
    fprintf(stderr, "IL["uint32FMTW(3)"] %6u %6u\n", i, IL.lo(i), IL.hi(i));

  fprintf(stderr, "Merge test 4a\n");
  IL.clear();
  IL.add(0, 25000);
  IL.add(25000, 25000);
  IL.add(50000, 50000);
  IL.merge();
  for (uint32 i=0; i<IL.numberOfIntervals(); i++)
    fprintf(stderr, "IL["uint32FMTW(3)"] %6u %6u\n", i, IL.lo(i), IL.hi(i));

  fprintf(stderr, "Merge test 4b\n");
  IL.clear();
  IL.add(    0, 25000, 1);
  IL.add(25000, 25000, 2);
  IL.add(50000, 50000, 4);
  IL.add(20000,  5000, 8);
  IL.add(45000,  5000, 16);
  IL.add(95000,  5000, 32);
  ID.depth(IL);
  IL.merge();
  for (uint32 i=0; i<IL.numberOfIntervals(); i++)
    fprintf(stderr, "IL["uint32FMTW(3)"] %6u %6u count %6u value %f\n", i, IL.lo(i), IL.hi(i), IL.count(i), IL.value(i));
  for (uint32 i=0; i<ID.numberOfIntervals(); i++)
    fprintf(stderr, "ID["uint32FMTW(3)"] %6u %6u depth %6u value %f\n", i, ID.lo(i), ID.hi(i), ID.count(i), ID.value(i));

  fprintf(stderr, "Merge test 5\n");
  IL.clear();
  ID.clear();
  IL.add(    0, 25000, 1);
  IL.add(25000, 25000, 2);
  IL.add(50000, 50000, 4);
  IL.add(20000, 20000, 8);
  IL.add(30000, 40000, 16);
  IL.add(50000, 10000, 32);
  //IL.add(10000, 90000, 32);
  ID.depth(IL);
  for (uint32 i=0; i<IL.numberOfIntervals(); i++)
    fprintf(stderr, "IL["uint32FMTW(3)"] %6u %6u count %6u value %f\n", i, IL.lo(i), IL.hi(i), IL.count(i), IL.value(i));
  for (uint32 i=0; i<ID.numberOfIntervals(); i++)
    fprintf(stderr, "ID["uint32FMTW(3)"] %6u %6u depth %6u value %f\n", i, ID.lo(i), ID.hi(i), ID.count(i), ID.value(i));

  fprintf(stderr, "Merge test 6 (same as 5, but val = default)\n");
  IL.clear();
  ID.clear();
  IL.add(    0, 25000);
  IL.add(25000, 25000);
  IL.add(50000, 50000);
  IL.add(20000, 20000);
  IL.add(30000, 40000);
  IL.add(50000, 10000);
  //IL.add(10000, 90000);
  ID.depth(IL);
  for (uint32 i=0; i<IL.numberOfIntervals(); i++)
    fprintf(stderr, "IL["uint32FMTW(3)"] %6u %6u count %6u value %f\n", i, IL.lo(i), IL.hi(i), IL.count(i), IL.value(i));
  for (uint32 i=0; i<ID.numberOfIntervals(); i++)
    fprintf(stderr, "ID["uint32FMTW(3)"] %6u %6u depth %6u value %f\n", i, ID.lo(i), ID.hi(i), ID.count(i), ID.value(i));


}



int
main(int argc, char **argv) {

  mt = mtInit(time(NULL));

  test();

  testIntersect(0);
  testIntersect(1);
  testIntersect(2);

  testMerge();

  exit(0);
}
