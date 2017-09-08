
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
 *    Sergey Koren beginning on 2016-FEB-24
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Brian P. Walenz beginning on 2016-OCT-24
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "ovStore.H"
#include "splitToWords.H"

#include <vector>

using namespace std;

int
main(int argc, char **argv) {
  char           *outName  = NULL;
  char           *gkpName  = NULL;
  bool		  partialOverlaps = false;
  uint32          minOverlapLength = 0;
  uint32          tolerance = 0;

  vector<char *>  files;

  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      outName = argv[++arg];

    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-tolerance") == 0) {
      tolerance = atoi(argv[++arg]);;

    } else if (strcmp(argv[arg], "-partial") == 0) {
      partialOverlaps = true;

    } else if (strcmp(argv[arg], "-len") == 0) {
      minOverlapLength = atoi(argv[++arg]);

    } else if (AS_UTL_fileExists(argv[arg])) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((err) || (gkpName == NULL) || (outName == NULL) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] file.mhap[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  Converts mhap native output to ovb\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o out.ovb     output file\n");
    fprintf(stderr, "\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR:  no gkpStore (-G) supplied\n");
    if (files.size() == 0)
      fprintf(stderr, "ERROR:  no overlap files supplied\n");

    exit(1);
  }

  char        *ovStr = new char [1024*1024];

  gkStore    *gkpStore = gkStore::gkStore_open(gkpName);
  ovOverlap   ov(gkpStore);
  ovFile      *of = new ovFile(NULL, outName, ovFileFullWrite);

  for (uint32 ff=0; ff<files.size(); ff++) {
    compressedFileReader  *in = new compressedFileReader(files[ff]);

    //  $1        $2     $3     $4     $5     $6         $7      $8    $9     $10      $11          $12        $13
    //  0         1      2      3      4      5          6       7     8      9        10           11         12
    //  0f1bd7b6  8189   1310   8014   +      b74d9367   14205   7340  14051  277      6711         255        cm:i:32
    //  0f1bd7b6  8189   1152   7272   -      a3026aca   7731    1642  7547   157      6120         255        cm:i:24
    //  aiid      alen   bgn    end    bori   biid       blen    bgn   end    #match   minimizers   alnlen     cm:i:errori

    while (fgets(ovStr, 1024*1024, in->file()) != NULL) {
      splitToWords  W(ovStr);

      ov.a_iid = W(0);
      ov.b_iid = W(5);

      if (ov.a_iid == ov.b_iid)
        continue;

      ov.dat.ovl.ahg5 = W(2);
      ov.dat.ovl.ahg3 = W(1) - W(3);

      if (W[4][0] == '+') {
        ov.dat.ovl.bhg5 = W(7);
        ov.dat.ovl.bhg3 = W(6) - W(8);
        ov.flipped(false);
      } else {
        ov.dat.ovl.bhg3 = W(7);
        ov.dat.ovl.bhg5 = W(6) - W(8);
        ov.flipped(true);
      }

      ov.erate(1-((double)W(9)/W(10)));

      //  Check the overlap - the hangs must be less than the read length.

      uint32  alen = gkpStore->gkStore_getRead(ov.a_iid)->gkRead_sequenceLength();
      uint32  blen = gkpStore->gkStore_getRead(ov.b_iid)->gkRead_sequenceLength();

      if ((alen < ov.dat.ovl.ahg5 + ov.dat.ovl.ahg3) ||
          (blen < ov.dat.ovl.bhg5 + ov.dat.ovl.bhg3)) {
        fprintf(stderr, "INVALID OVERLAP %8u (len %6d) %8u (len %6d) hangs %6lu %6lu - %6lu %6lu flip %lu\n",
                ov.a_iid, alen,
                ov.b_iid, blen,
                ov.dat.ovl.ahg5, ov.dat.ovl.ahg3,
                ov.dat.ovl.bhg5, ov.dat.ovl.bhg3,
                ov.dat.ovl.flipped);
        exit(1);
      }

      if (!ov.overlapIsDovetail() && partialOverlaps == false) {
         if (alen <= blen && ov.dat.ovl.ahg5 >= 0 && ov.dat.ovl.ahg3 >= 0 && ov.dat.ovl.bhg5 >= ov.dat.ovl.ahg5 && ov.dat.ovl.bhg3 >= ov.dat.ovl.ahg3 && ((ov.dat.ovl.ahg5 + ov.dat.ovl.ahg3)) < tolerance) {
              ov.dat.ovl.bhg5 = max(0, ov.dat.ovl.bhg5 - ov.dat.ovl.ahg5); ov.dat.ovl.ahg5 = 0;
              ov.dat.ovl.bhg3 = max(0, ov.dat.ovl.bhg3 - ov.dat.ovl.ahg3); ov.dat.ovl.ahg3 = 0;
           }
           // second is b contained (both b hangs can be extended)
           //
           else if (alen >= blen && ov.dat.ovl.bhg5 >= 0 && ov.dat.ovl.bhg3 >= 0 && ov.dat.ovl.ahg5 >= ov.dat.ovl.bhg5 && ov.dat.ovl.ahg3 >= ov.dat.ovl.bhg3 && ((ov.dat.ovl.bhg5 + ov.dat.ovl.bhg3)) < tolerance) {
              ov.dat.ovl.ahg5 = max(0, ov.dat.ovl.ahg5 - ov.dat.ovl.bhg5); ov.dat.ovl.bhg5 = 0;
              ov.dat.ovl.ahg3 = max(0, ov.dat.ovl.ahg3 - ov.dat.ovl.bhg3); ov.dat.ovl.bhg3 = 0;
           }
           // third is 5' dovetal  ---------->
           //                          ---------->
           //                          or
           //                          <---------
           //                         bhg5 here is always first overhang on b read
           //
           else if (ov.dat.ovl.ahg3 <= ov.dat.ovl.bhg3 && (ov.dat.ovl.ahg3 >= 0 && ((double)(ov.dat.ovl.ahg3)) < tolerance) &&
                   (ov.dat.ovl.bhg5 >= 0 && ((double)(ov.dat.ovl.bhg5)) < tolerance)) {
              ov.dat.ovl.ahg5 = max(0, ov.dat.ovl.ahg5 - ov.dat.ovl.bhg5); ov.dat.ovl.bhg5 = 0;
              ov.dat.ovl.bhg3 = max(0, ov.dat.ovl.bhg3 - ov.dat.ovl.ahg3); ov.dat.ovl.ahg3 = 0;
           }
           //
           // fourth is 3' dovetail    ---------->
           //                     ---------->
           //                     or
           //                     <----------
           //                     bhg5 is always first overhang on b read
           else if (ov.dat.ovl.ahg5 <= ov.dat.ovl.bhg5 && (ov.dat.ovl.ahg5 >= 0 && ((double)(ov.dat.ovl.ahg5)) < tolerance) &&
                   (ov.dat.ovl.bhg3 >= 0 && ((double)(ov.dat.ovl.bhg3)) < tolerance)) {
              ov.dat.ovl.bhg5 = max(0, ov.dat.ovl.bhg5 - ov.dat.ovl.ahg5); ov.dat.ovl.ahg5 = 0;
              ov.dat.ovl.ahg3 = max(0, ov.dat.ovl.ahg3 - ov.dat.ovl.bhg3); ov.dat.ovl.bhg3 = 0;
           }
     }

      ov.dat.ovl.forUTG = (partialOverlaps == false) && (ov.overlapIsDovetail() == true);;
      ov.dat.ovl.forOBT = partialOverlaps;
      ov.dat.ovl.forDUP = partialOverlaps;

      // check the length is big enough
      if (ov.a_end() - ov.a_bgn() < minOverlapLength || ov.b_end() - ov.b_bgn() < minOverlapLength) {
         continue;
      }

      //  Overlap looks good, write it!

      of->writeOverlap(&ov);
    }

    arg++;
  }

  delete    of;
  delete [] ovStr;

  gkpStore->gkStore_close();

  exit(0);
}
