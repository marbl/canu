
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
#include "mt19937ar.H"

//  g++ -o intervalListTest -I.. -I. intervalListTest.C

int
main(int argc, char **argv) {

  AS_configure(argc, argv);

  //  Boring basic test.
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
    mtRandom  mt(strtouint32(argv[1]));

    //  About 6.5 minutes per million, so this should be about an hour.
    uint32  iterMax = 935000;

    for (uint32 iter=0; iter<iterMax; iter++) {
      uint32  numIntervals =     mt.mtRandom32() % 5000;
      uint32  maxLen       = 1 + mt.mtRandom32() % 1000;
      uint32  maxBgn       = 1 + mt.mtRandom32() % 50000;
      uint32 *depth        = new uint32 [maxBgn + maxLen];

      memset(depth, 0, sizeof(uint32) * (maxBgn + maxLen));

      if (iter % 1000 == 0)
        fprintf(stderr, "%10u/%10u: %3u intervals, each up to %4u long, coords up to %4u\n",
                iter, iterMax,
                numIntervals, maxLen, maxBgn);

      intervalList<uint32>  il;

      //  Add intervals to the list.
      //  Sum depths explicitly.
      for (uint32 ii=0; ii<numIntervals; ii++) {
        uint32  bgn = mt.mtRandom32() % maxBgn;   //  bgn between 0 and maxBgn
        uint32  len = mt.mtRandom32() % maxLen;   //  len between 0 and maxLen

        il.add(bgn, len);

        for (uint32 xx=bgn; xx<bgn+len; xx++)
          depth[xx]++;

        //fprintf(stderr, "IL %u - %u\n", bgn, bgn+len);
      }

      //  Convert intervals to depths.
      intervalDepth<uint32>  de(il);

      //  Over all the depth regions, subtract the computed depth from
      //  the explicit depth.
      for (uint32 xx=0; xx<de.numberOfIntervals(); xx++) {
        uint32  bgn = de.lo(xx);
        uint32  end = de.hi(xx);
        uint32  dpt = de.depth(xx);

        //fprintf(stderr, "ID %u - %u depth %u\n", bgn, end, dpt);

        for (uint32 cc=bgn; cc<end; cc++) {
          //if (cc < 30)
          //  fprintf(stderr, "depth[%u] = %u -> %u\n", cc, depth[cc], depth[cc] - dpt);
          depth[cc] -= dpt;
        }
      }

      //  Every explicit depth should now be zero, even the ones
      //  not covered.
      for (uint32 cc=0; cc<maxBgn + maxLen; cc++) {
        if (depth[cc] != 0)
          fprintf(stderr, "ERROR: depth[%u] = %u in iter %u\n", cc, depth[cc], iter);
        assert(depth[cc] == 0);
      }

      delete [] depth;
    }
  }

  exit(0);
}
