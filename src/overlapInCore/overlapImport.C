
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
#include "sqStore.H"
#include "ovStore.H"

#include "AS_UTL_decodeRange.H"

#include "splitToWords.H"
#include "mt19937ar.H"

#include <vector>

using namespace std;

#define  TYPE_NONE    'N'
#define  TYPE_LEGACY  'L'
#define  TYPE_COORDS  'C'
#define  TYPE_HANGS   'H'
#define  TYPE_RAW     'R'
#define  TYPE_OVB     'O'
#define  TYPE_RANDOM  'r'


int
main(int argc, char **argv) {
  char                  *seqStoreName = NULL;
  sqStore               *seqStore = NULL;

  char                  *ovlFileName = NULL;
  char                  *ovlStoreName = NULL;

  char                   inType = TYPE_NONE;

  uint64                 rmin = 0, rmax = 0;
  uint32                 abgn = 1, aend = 0;
  uint32                 bbgn = 1, bend = 0;

  bool                   native = false;

  vector<char *>         files;


  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStoreName = argv[++arg];

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

    } else if (strcmp(argv[arg], "-ovb") == 0) {
      fprintf(stderr, "-ovb not implemented.\n"), exit(1);
      inType = TYPE_OVB;

    } else if (strcmp(argv[arg], "-random") == 0) {
      inType    = TYPE_RANDOM;
      AS_UTL_decodeRange(argv[++arg], rmin, rmax);
      files.push_back(NULL);

    } else if (strcmp(argv[arg], "-a") == 0) {
      AS_UTL_decodeRange(argv[++arg], abgn, aend);

    } else if (strcmp(argv[arg], "-b") == 0) {
      AS_UTL_decodeRange(argv[++arg], bbgn, bend);

    } else if (strcmp(argv[arg], "-native") == 0) {
      native = true;

    } else if ((strcmp(argv[arg], "-") == 0) ||
               (AS_UTL_fileExists(argv[arg]))) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if (seqStoreName == NULL)
    err++;
  if (inType == TYPE_NONE)
    err++;

  if ((err) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] ascii-ovl-file-input.[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -S name.seqStore   path to valid sequence store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output Format:\n");
    fprintf(stderr, "  -o file.ovb        output file name\n");
    fprintf(stderr, "  -O name.ovlStore   output overlap store");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input Format:\n");
    fprintf(stderr, "  -legacy            'CA8 overlapStore -d' format\n");
    fprintf(stderr, "  -coords            'overlapConvert -coords' format (not implemented)\n");
    fprintf(stderr, "  -hangs             'overlapConvert -hangs' format (not implemented)\n");
    fprintf(stderr, "  -raw               'overlapConvert -raw' format\n");
    fprintf(stderr, "  -ovb               'overlapInCore' format (not implemented)\n");
    fprintf(stderr, "  -random N          create N random overlaps, for store testing\n");
    fprintf(stderr, "    -a x-y             A read IDs will be between x and y\n");
    fprintf(stderr, "    -b x-y             B read IDs will be between x and y\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -native            output ovb (-o) files will not be snappy compressed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input file can be stdin ('-') or a gz/bz2/xz compressed file.\n");
    fprintf(stderr, "\n");

    if (seqStoreName == NULL)
      fprintf(stderr, "ERROR: need to supply a seqStore (-S).\n");
    if (inType == TYPE_NONE)
      fprintf(stderr, "ERROR: need to supply a format type (-legacy, -coords, -hangs, -raw).\n");
    if (files.size() == 0)
      fprintf(stderr, "ERROR: need to supply input files.\n");

    exit(1);
  }

  if (seqStoreName)
    seqStore = sqStore::sqStore_open(seqStoreName);

  char          *S     = new char [1024];
  splitToWords   W;
  ovOverlap     ov(seqStore);

  ovFile        *of = (ovlFileName  == NULL) ? NULL : new ovFile(seqStore, ovlFileName, ovFileFullWrite);
  ovStoreWriter *os = (ovlStoreName == NULL) ? NULL : new ovStoreWriter(ovlStoreName, seqStore);

  if ((of) && (native == true))
    of->enableSnappy(false);

  //  Make random inputs first.

  if (inType == TYPE_RANDOM) {
    mtRandom  mt;

    if (aend == 0)   aend = seqStore->sqStore_getNumReads();
    if (bend == 0)   bend = seqStore->sqStore_getNumReads();

    uint64  numRandom = rmin + floor(mt.mtRandomRealOpen() * (rmax - rmin));

    for (uint64 ii=0; ii<numRandom; ii++) {
      uint32   aID      = abgn + floor(mt.mtRandomRealOpen() * (aend - abgn));
      uint32   bID      = bbgn + floor(mt.mtRandomRealOpen() * (bend - bbgn));

      if (aID < bID) {
        uint32 t = aID;
        aID = bID;
        bID = aID;
      }

#if 0
      //  For testing when reads have no overlaps in store building.  Issue #302.
      aID = aID & 0xfffffff0;
      bID = bID & 0xfffffff0;

      if (aID == 0)   aID = 1;
      if (bID == 0)   bID = 1;
#endif

      uint32   aLen     = seqStore->sqStore_getRead(aID)->sqRead_sequenceLength();
      uint32   bLen     = seqStore->sqStore_getRead(bID)->sqRead_sequenceLength();

      bool     olapFlip = mt.mtRandom32() % 2;

      //  We could be fancy and make actual overlaps that make sense, or punt and make overlaps that
      //  are valid but nonsense.

      ov.a_iid = aID;
      ov.b_iid = bID;

      ov.flipped(olapFlip);

      ov.a_hang((int32)(mt.mtRandomRealOpen() * 2 * aLen - aLen));
      ov.b_hang((int32)(mt.mtRandomRealOpen() * 2 * bLen - bLen));

      ov.dat.ovl.forOBT = false;
      ov.dat.ovl.forDUP = false;
      ov.dat.ovl.forUTG = true;

      ov.erate(mt.mtRandomRealOpen() * 0.1);

      if (of)
        of->writeOverlap(&ov);

      if (os)
        os->writeOverlap(&ov);
    }

    files.pop_back();
  }

  //  Now process any files.

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

  seqStore->sqStore_close();

  exit(0);
}
