
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"
#include "sqStore.H"
#include "ovStore.H"

#include "strings.H"
#include "mt19937ar.H"

#include <vector>

using namespace std;


int
main(int argc, char **argv) {
  char const            *seqStoreName = NULL;

  char const            *ovlFileName = NULL;
  char const            *ovlStoreName = NULL;

  bool                   asCoords    = false;    //  Format of input overlaps
  bool                   asHangs     = false;
  bool                   asUnaligned = false;
  bool                   asPAF       = false;
  bool                   asRandom    = false;

  uint64                 rmin = 0, rmax = 0;
  uint32                 abgn = 1, aend = 0;
  uint32                 bbgn = 1, bend = 0;

  vector<char const *>   files;


  vector<char const *>  err;
  int32                 arg = 1;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-S") == 0) {
      seqStoreName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-o") == 0) {
      ovlFileName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-O") == 0) {
      ovlStoreName = argv[++arg];
    }


    else if (strcmp(argv[arg], "-raw") == 0)
      sqRead_setDefaultVersion(sqRead_raw);

    else if (strcmp(argv[arg], "-obt") == 0)
      sqRead_setDefaultVersion(sqRead_corrected);

    else if (strcmp(argv[arg], "-utg") == 0)
      sqRead_setDefaultVersion(sqRead_trimmed);


    else if (strcmp(argv[arg], "-coords") == 0) {
      asCoords    = true;
      asHangs     = false;
      asUnaligned = false;
      asPAF       = false;
      asRandom    = false;
    }

    else if (strcmp(argv[arg], "-hangs") == 0) {
      asCoords    = false;
      asHangs     = true;
      asUnaligned = false;
      asPAF       = false;
      asRandom    = false;
    }

    else if (strcmp(argv[arg], "-unaligned") == 0) {
      asCoords    = false;
      asHangs     = false;
      asUnaligned = true;
      asPAF       = false;
      asRandom    = false;
    }

    else if (strcmp(argv[arg], "-paf") == 0) {
      asCoords    = false;
      asHangs     = false;
      asUnaligned = false;
      asPAF       = true;
      asRandom    = false;
    }

    else if (strcmp(argv[arg], "-random") == 0) {
      asCoords    = false;
      asHangs     = false;
      asUnaligned = false;
      asPAF       = false;
      asRandom    = true;

      decodeRange(argv[++arg], rmin, rmax);
    }

    else if (strcmp(argv[arg], "-a") == 0) {
      decodeRange(argv[++arg], abgn, aend);
    }

    else if (strcmp(argv[arg], "-b") == 0) {
      decodeRange(argv[++arg], bbgn, bend);
    }

    else if ((strcmp(argv[arg], "-") == 0) ||
             (fileExists(argv[arg]))) {
      files.push_back(argv[arg]);
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err.push_back(s);
    }

    arg++;
  }


  if (seqStoreName == NULL)
    err.push_back("ERROR: no input seqStore (-S) supplied.\n");

  if ((files.size() == 0) && (asRandom == false))
    err.push_back("ERROR: no input overlap files supplied.\n");


  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [options] ascii-ovl-file-input.[.gz]\n", argv[0]);
    fprintf(stderr, "  -S name.seqStore    path to valid sequence store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUT FORMAT (pick exactly one):\n");
    fprintf(stderr, "  -o file.ovb         output file name\n");
    fprintf(stderr, "  -O name.ovlStore    output overlap store");
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUT FORMAT:\n");
    fprintf(stderr, "  -coords             as coordinates on each read (default)\n");
    fprintf(stderr, "  -hangs              as dovetail hangs\n");
    fprintf(stderr, "  -unaligned          as unaligned regions on each read\n");
    fprintf(stderr, "  -paf                as miniasm Pairwise mApping Format\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "READ VERSION:\n");
    fprintf(stderr, "  -raw                uncorrected raw reads\n");
    fprintf(stderr, "  -obt                corrected reads\n");
    fprintf(stderr, "  -utg                trimmed reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "RANDOM OPTIONS:\n");
    fprintf(stderr, "  Doesn't read overlaps from an input file, instead generated\n");
    fprintf(stderr, "  random non-sense overlaps.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -random N           create N random overlaps, for store testing\n");
    fprintf(stderr, "  -a x-y              A read IDs will be between x and y\n");
    fprintf(stderr, "  -b x-y              B read IDs will be between x and y\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input file can be stdin ('-') or a gz/bz2/xz compressed file.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  sqStore       *seqStore = new sqStore(seqStoreName);

  char          *S     = new char [1024];
  splitToWords   W;
  ovOverlap     ov;

  ovFile        *of = (ovlFileName  == NULL) ? NULL : new ovFile(seqStore, ovlFileName, ovFileFullWrite);
  ovStoreWriter *os = (ovlStoreName == NULL) ? NULL : new ovStoreWriter(ovlStoreName, seqStore);

  //  Make random inputs first.

  if (asRandom == true) {
    mtRandom  mt;

    if (aend == 0)   aend = seqStore->sqStore_lastReadID();
    if (bend == 0)   bend = seqStore->sqStore_lastReadID();

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

      uint32   aLen     = seqStore->sqStore_getReadLength(aID);
      uint32   bLen     = seqStore->sqStore_getReadLength(bID);

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
  }

  //  Now process any files.

  for (uint32 ff=0; ff<files.size(); ff++) {
    compressedFileReader   *in = new compressedFileReader(files[ff]);

    fgets(S, 1024, in->file());

    while (!feof(in->file())) {
      W.split(S);

      if (asCoords) {
        ov.fromString(W, ovOverlapAsCoords);
      }

      else if (asHangs) {
        ov.fromString(W, ovOverlapAsHangs);
      }

      else if (asUnaligned) {
        ov.fromString(W, ovOverlapAsUnaligned);
      }

      else if (asPAF) {
        ov.fromString(W, ovOverlapAsPaf);
      }

      else {
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

  delete seqStore;

  exit(0);
}
