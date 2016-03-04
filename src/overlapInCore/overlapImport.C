
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
 *    Brian P. Walenz from 2015-MAY-14 to 2015-JUN-25
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-17
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"

#include "splitToWords.H"

#include <vector>

using namespace std;

#define  TYPE_NONE    'N'
#define  TYPE_LEGACY  'L'
#define  TYPE_COORDS  'C'
#define  TYPE_HANGS   'H'
#define  TYPE_RAW     'R'



int
main(int argc, char **argv) {
  char                  *gkpStoreName = NULL;
  gkStore               *gkpStore = NULL;

  char                  *ovlFileName = NULL;
  char                  *ovlStoreName = NULL;

  char                   inType = TYPE_NONE;

  vector<char *>         files;


  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      ovlFileName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-legacy") == 0) {
      inType = TYPE_LEGACY;

    } else if (strcmp(argv[arg], "-coords") == 0) {
      fprintf(stderr, "-coords not implemented.\n"), exit(1);
      inType = TYPE_COORDS;

    } else if (strcmp(argv[arg], "-hangs") == 0) {
      fprintf(stderr, "-hangs not implemented.\n"), exit(1);
      inType = TYPE_HANGS;

    } else if (strcmp(argv[arg], "-raw") == 0) {
      inType = TYPE_RAW;

    } else if ((strcmp(argv[arg], "-") == 0) ||
               (AS_UTL_fileExists(argv[arg]))) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if (gkpStoreName == NULL)
    err++;
  if (inType == TYPE_NONE)
    err++;

  if ((err) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] ascii-ovl-file-input.[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -G name.gkpStore   path to valid gatekeeper store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "  -o file.ovb        output file name\n");
    fprintf(stderr, "  -O name.ovlStore   output overlap store");
    fprintf(stderr, "\n");
    fprintf(stderr, "Format options:\n");
    fprintf(stderr, "  -legacy            'CA8 overlapStore -d' format\n");
    fprintf(stderr, "  -coords            'overlapConvert -coords' format (not implemented)\n");
    fprintf(stderr, "  -hangs             'overlapConvert -hangs' format (not implemented)\n");
    fprintf(stderr, "  -raw               'overlapConvert -raw' format\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input file can be stdin ('-') or a gz/bz2/xz compressed file.\n");
    fprintf(stderr, "\n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: need to supply a gkpStore (-G).\n");
    if (inType == TYPE_NONE)
      fprintf(stderr, "ERROR: need to supply a format type (-legacy, -coords, -hangs, -raw).\n");

    exit(1);
  }

  if (gkpStoreName)
    gkpStore = gkStore::gkStore_open(gkpStoreName);

  char         *S     = new char [1024];
  splitToWords  W;
  ovOverlap    ov(gkpStore);

  ovFile       *of    = (ovlFileName  == NULL) ? NULL : new ovFile(ovlFileName, ovFileFullWrite);
  ovStore      *os    = (ovlStoreName == NULL) ? NULL : new ovStore(ovlStoreName, gkpStore, ovStoreWrite);

  for (uint32 ff=0; ff<files.size(); ff++) {
    compressedFileReader   *in = new compressedFileReader(files[ff]);

    fgets(S, 1024, in->file());

    while (!feof(in->file())) {
      W.split(S);

      switch (inType) {
      case TYPE_LEGACY:
        //  Aiid Biid 'I/N' ahang bhang erate erate
        ov.a_iid = W(0);
        ov.b_iid = W(1);

        ov.flipped(W[2][0] == 'I');

        ov.a_hang(W(3));
        ov.b_hang(W(4));

        //  Overlap store reports %error, but we expect fraction error.
        //ov.erate(atof(W[5]);  //  Don't use the original uncorrected error rate
        ov.erate(atof(W[6]) / 100.0);
        break;

      case TYPE_COORDS:
        break;

      case TYPE_HANGS:
        break;

      case TYPE_RAW:
         ov.a_iid = W(0);
         ov.b_iid = W(1);

         ov.flipped(W[2][0] == 'I');

         ov.dat.ovl.span = W(3);

         ov.dat.ovl.ahg5 = W(4);
         ov.dat.ovl.ahg3 = W(5);

         ov.dat.ovl.bhg5 = W(6);
         ov.dat.ovl.bhg3 = W(7);

         ov.erate(atof(W[8]) / 1);

         ov.dat.ovl.forUTG = false;
         ov.dat.ovl.forOBT = false;
         ov.dat.ovl.forDUP = false;

         for (uint32 i = 9; i < W.numWords(); i++) {
           ov.dat.ovl.forUTG |= ((W[i][0] == 'U') && (W[i][1] == 'T') && (W[i][2] == 'G'));  //  Fails if W[i] == "U".
           ov.dat.ovl.forOBT |= ((W[i][0] == 'O') && (W[i][1] == 'B') && (W[i][2] == 'T'));
           ov.dat.ovl.forDUP |= ((W[i][0] == 'D') && (W[i][1] == 'U') && (W[i][2] == 'P'));
         }
        break;

      default:
        break;
      }

      if (of)
        of->writeOverlap(&ov);

      if (os)
        os->writeOverlap(&ov);

      fgets(S, 1024, in->file());
    }

    delete in;
  }

  delete    os;
  delete    of;

  delete [] S;

  gkpStore->gkStore_close();

  exit(0);
}
