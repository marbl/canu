
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
 *    Brian P. Walenz from 2015-MAR-27 to 2015-JUN-25
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Sergey Koren beginning on 2016-FEB-24
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
  bool            asCoords = true;

  char           *outName     = NULL;

  uint32          baseIDhash  = 0;
  uint32          numIDhash   = 0;
  uint32          baseIDquery = 0;

  vector<char *>  files;


  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      outName = argv[++arg];

    } else if (strcmp(argv[arg], "-h") == 0) {
      baseIDhash = atoi(argv[++arg]) - 1;
      numIDhash = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-q") == 0) {
      baseIDquery = atoi(argv[++arg]) - 1;

    } else if (AS_UTL_fileExists(argv[arg])) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((err) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] file.mhap[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  Converts mhap native output to ovb\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o out.ovb     output file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -h id num      base id and number of hash table reads\n");
    fprintf(stderr, "                   (mhap output IDs 1 through 'num')\n");
    fprintf(stderr, "  -q id          base id of query reads\n");
    fprintf(stderr, "                   (mhap output IDs 'num+1' and higher)\n");

    if (files.size() == 0)
      fprintf(stderr, "ERROR:  no overlap files supplied\n");

    exit(1);
  }

  char        *ovStr = new char [1024];

  ovOverlap   ov(NULL);
  ovFile      *of = new ovFile(outName, ovFileFullWrite);

  for (uint32 ff=0; ff<files.size(); ff++) {
    compressedFileReader  *in = new compressedFileReader(files[ff]);

    //  $1    $2   $3       $4  $5  $6  $7   $8   $9  $10 $11  $12
    //  0     1    2        3   4   5   6    7    8   9   10   11
    //  26887 4509 87.05933 301 0   479 2305 4328 1   34  1852 3637
    //  aiid  biid qual     ?   ori bgn end  len  ori bgn end  len

    while (fgets(ovStr, 1024, in->file()) != NULL) {
      splitToWords  W(ovStr);

      ov.a_iid = W(0) + baseIDquery - numIDhash;  //  First ID is the query
      ov.b_iid = W(1) + baseIDhash;               //  Second ID is the hash table

      if (ov.a_iid == ov.b_iid)
        continue;

      assert(W[4][0] == '0');

      ov.dat.ovl.forUTG = true;
      ov.dat.ovl.forOBT = true;
      ov.dat.ovl.forDUP = true;

      ov.dat.ovl.ahg5 = W(5);
      ov.dat.ovl.ahg3 = W(7) - W(6);

      if (W[8][0] == '0') {
        ov.dat.ovl.bhg5 = W(9);
        ov.dat.ovl.bhg3 = W(11) - W(10);
        ov.flipped(false);
      } else {
        ov.dat.ovl.bhg3 = W(9);
        ov.dat.ovl.bhg5 = W(11) - W(10);
        ov.flipped(true);
      }

      ov.erate(atof(W[2]));

      of->writeOverlap(&ov);
    }

    arg++;
  }

  delete    of;
  delete [] ovStr;

  exit(0);
}
